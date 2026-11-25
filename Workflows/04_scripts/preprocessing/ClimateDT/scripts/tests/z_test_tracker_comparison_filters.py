import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from z_test_comparison_track import (
    add_aligned_wind,
    compare_climatedt_tracks,
    create_aligned_tracks_dataframe,
    plot_comparison_tracks,
    save_comparison_figures,
)
from z_tracker_comparison import align_tracks, compare_all_strategies


class ComparisonFilterTests(unittest.TestCase):
    def setUp(self):
        times = pd.date_range("2020-01-01", periods=5, freq="h")
        self.reference = pd.DataFrame({
            "time": times,
            "lat": [-39.0, -40.0, -41.0, -38.0, -39.0],
            "lon": [40.0] * 5,
            "min_pres": [1000.0] * 5,
            "usa_sshs": [0.0, 1.0, 3.0, np.nan, 2.0],
        })
        self.tracked = pd.DataFrame({
            "time": times,
            "lat": [-39.0, -40.0, -41.0, -38.0, -39.0],
            "lon": [40.0] * 5,
            "msl": [1000.0, 1000.0, 900.0, 1000.0, 1000.0],
        })

    def test_create_aligned_tracks_dataframe_for_inspection(self):
        tracks = {
            ("strategy", "hist", "1", "storm"): self.tracked,
            ("cont", "2", "storm"): self.tracked,
        }
        reference_before = self.reference.copy(deep=True)
        tracked_before = self.tracked.copy(deep=True)

        aligned = create_aligned_tracks_dataframe(
            {"storm": self.reference}, tracks, tracking_strategy="single",
        )

        self.assertEqual(len(aligned), 10)
        self.assertEqual(aligned["strategy"].unique().tolist(), ["strategy", "single"])
        self.assertEqual(aligned["should_exclude"].tolist(), [False, False, True, False, False] * 2)
        pd.testing.assert_series_equal(
            aligned["reference_time_ref"].reset_index(drop=True),
            aligned["time"].reset_index(drop=True),
            check_names=False,
        )
        pd.testing.assert_frame_equal(self.reference, reference_before)
        pd.testing.assert_frame_equal(self.tracked, tracked_before)

    def test_create_aligned_tracks_preserves_matched_reference_time(self):
        track = self.tracked.iloc[:3].copy()
        track["time"] += pd.Timedelta(minutes=10)
        aligned = create_aligned_tracks_dataframe(
            {"storm": self.reference}, {("strategy", "hist", "1", "storm"): track},
        )
        self.assertEqual(aligned["reference_time_ref"].tolist(), self.reference["time"].iloc[:3].tolist())
        self.assertEqual(aligned["time"].tolist(), track["time"].tolist())

    def test_add_aligned_wind_matches_track_key_and_time(self):
        key = ("all_ibtracs", "hist", "1", "storm")
        other_key = ("other", "hist", "1", "storm")
        aligned = create_aligned_tracks_dataframe(
            {"storm": self.reference}, {key: self.tracked, other_key: self.tracked},
        )
        wind = self.tracked.iloc[[2, 0]].copy()
        wind["max_wind_model"] = [20.0, 10.0]
        result = add_aligned_wind(aligned, {key: wind}, "200km")
        result = add_aligned_wind(result, {key: wind.assign(max_wind_model=[8.0, 4.0])}, "rmw")
        self.assertEqual(result["max_wind_model_200km"].iloc[0], 10.0)
        self.assertEqual(result["max_wind_model_200km"].iloc[2], 20.0)
        self.assertEqual(result["max_wind_model_rmw"].iloc[2], 8.0)
        self.assertTrue(result["max_wind_model_200km"].iloc[5:].isna().all())
        self.assertTrue(pd.isna(result["max_wind_model_200km"].iloc[1]))
        self.assertNotIn("max_wind_model_200km", aligned.columns)

    def test_add_aligned_wind_rejects_duplicate_or_mismatched_points(self):
        key = ("all_ibtracs", "hist", "1", "storm")
        aligned = create_aligned_tracks_dataframe({"storm": self.reference}, {key: self.tracked})
        wind = self.tracked.assign(max_wind_model=10.0)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            add_aligned_wind(aligned, {key: pd.concat([wind, wind.iloc[[0]]])}, "200km")
        wind.loc[0, "lat"] += 1
        with self.assertRaisesRegex(ValueError, "coordinates"):
            add_aligned_wind(aligned, {key: wind}, "200km")

    def test_plot_comparison_tracks_filters_reference_points(self):
        key = ("strategy", "hist", "1", "storm")
        results = pd.DataFrame({"storm": ["storm"], "corr_pres": [0.5]})
        original = self.tracked.copy(deep=True)
        with tempfile.TemporaryDirectory() as directory:
            with patch("test_comparison_track.plot_tracks_all") as plot:
                plot_comparison_tracks(
                    results, {"storm": self.reference}, {key: self.tracked}, directory,
                    plot_mode="all", min_ref_lat=-40, min_ref_sshs=1,
                )
                plotted_track = plot.call_args.args[1][key]
                self.assertEqual(plotted_track["time"].tolist(), self.tracked["time"].iloc[[1, 4]].tolist())
                self.assertEqual(Path(plot.call_args.args[2]), Path(directory) / "figures" / "tracking")
            pd.testing.assert_frame_equal(self.tracked, original)

    def test_plot_comparison_tracks_low_correlation_and_none(self):
        results = pd.DataFrame({"storm": ["storm", "other"], "corr_pres": [0.1, 0.9]})
        tracks = {
            ("strategy", "hist", "1", "storm"): self.tracked,
            ("strategy", "hist", "1", "other"): self.tracked,
        }
        with tempfile.TemporaryDirectory() as directory:
            with patch("test_comparison_track.plot_tracks_all") as plot:
                plot_comparison_tracks(
                    results, {"storm": self.reference, "other": self.reference}, tracks,
                    directory, plot_mode="none",
                )
                plot.assert_not_called()
                plot_comparison_tracks(
                    results, {"storm": self.reference, "other": self.reference}, tracks,
                    directory, plot_mode="low_correlation", min_ref_lat=None,
                )
                self.assertEqual(list(plot.call_args.args[1]), [("strategy", "hist", "1", "storm")])
            with self.assertRaises(ValueError):
                plot_comparison_tracks(results, {}, {}, directory, plot_mode="invalid")

    def test_legacy_plot_mode_delegates_to_separate_function(self):
        tracks = {("strategy", "hist", "1", "storm"): self.tracked}
        comparison = pd.DataFrame({"storm": ["storm"], "strategy": ["strategy"]})
        with tempfile.TemporaryDirectory() as directory:
            with patch("test_comparison_track.compare_all_strategies", return_value=comparison):
                with patch("test_comparison_track.plot_comparison_tracks") as plot:
                    results, _ = compare_climatedt_tracks(
                        {"storm": self.reference}, tracks, results_dir=directory,
                    )
                    plot.assert_not_called()
                    compare_climatedt_tracks(
                        {"storm": self.reference}, tracks, results_dir=directory,
                        plot_mode="all",
                    )
                    self.assertEqual(plot.call_args.args[4], "all")
                    self.assertEqual(results.iloc[0]["storm"], "storm")

    def test_filters_matched_reference_points_inclusive_without_rematching(self):
        latitude_only = align_tracks(self.tracked, self.reference, min_ref_lat=-40)
        self.assertEqual(latitude_only["time"].tolist(), self.tracked["time"].iloc[[0, 1, 3, 4]].tolist())
        category_only = align_tracks(self.tracked, self.reference, min_ref_sshs=1)
        self.assertEqual(category_only["time"].tolist(), self.tracked["time"].iloc[[1, 2, 4]].tolist())
        aligned = align_tracks(self.tracked, self.reference, min_ref_lat=-40, min_ref_sshs=1)
        self.assertEqual(aligned["time"].tolist(), self.tracked["time"].iloc[[1, 4]].tolist())
        self.assertEqual(aligned["lat_ref"].tolist(), [-40.0, -39.0])

        result = compare_all_strategies(
            "storm", {"strategy": self.tracked}, self.reference,
            min_ref_lat=-40, min_ref_sshs=1,
        ).iloc[0]
        self.assertEqual(result["n_points_compared"], 2)
        self.assertEqual(result["RMSE_pres_hPa"], 0)
        self.assertEqual(result["track_duration_hrs"], 3)

    def test_disabled_filters_preserve_unfiltered_metrics(self):
        result = compare_all_strategies("storm", {"strategy": self.tracked}, self.reference).iloc[0]
        self.assertEqual(result["n_points_compared"], 5)
        self.assertAlmostEqual(result["RMSE_pres_hPa"], np.sqrt(2000))
        self.assertEqual(result["track_duration_hrs"], 4)

    def test_nearest_reference_outside_tracked_time_range_is_used(self):
        reference = pd.DataFrame({
            "time": pd.to_datetime(["2020-01-01 00:00", "2020-01-01 03:00", "2020-01-01 06:00"]),
            "lat": [-41.0, -39.0, -41.0],
            "lon": [40.0] * 3,
            "min_pres": [1000.0] * 3,
        })
        tracked = pd.DataFrame({
            "time": pd.to_datetime(["2020-01-01 01:00", "2020-01-01 03:00", "2020-01-01 05:00"]),
            "lat": [-39.0] * 3,
            "lon": [40.0] * 3,
            "msl": [1000.0] * 3,
        })
        aligned = align_tracks(tracked, reference, min_ref_lat=-40)
        self.assertEqual(aligned["time"].tolist(), [tracked["time"].iloc[1]])

    def test_empty_filter_warns_and_skips(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = compare_all_strategies(
                "storm", {"strategy": self.tracked}, self.reference, min_ref_lat=0,
            )
        self.assertTrue(result.empty)
        self.assertTrue(any("no aligned points" in str(w.message) for w in caught))

    def test_single_overlapping_point_is_compared(self):
        result = compare_all_strategies(
            "storm", {"strategy": self.tracked.iloc[[1]]}, self.reference,
            min_ref_lat=-40,
        ).iloc[0]
        self.assertEqual(result["n_points_compared"], 1)
        self.assertEqual(result["track_duration_hrs"], 0)

    def test_comparison_default_and_optional_plotting(self):
        tracks = {("strategy", "hist", "1", "storm"): self.tracked}
        with tempfile.TemporaryDirectory() as directory:
            with patch("test_comparison_track.plot_tracks_all") as plot:
                results, _ = compare_climatedt_tracks(
                    {"storm": self.reference}, tracks, results_dir=directory,
                )
                plot.assert_not_called()
                self.assertEqual(results.iloc[0]["n_points_compared"], 4)
                filtered_results, _ = compare_climatedt_tracks(
                    {"storm": self.reference}, tracks, results_dir=directory,
                    plot_mode="all", min_ref_lat=None, min_ref_sshs=1,
                )
                plot.assert_called_once()
                self.assertEqual(filtered_results.iloc[0]["n_points_compared"], 3)
                plotted_track = plot.call_args.args[1][next(iter(tracks))]
                self.assertEqual(plotted_track["time"].tolist(), self.tracked["time"].iloc[[1, 2, 4]].tolist())

            with self.assertRaises(ValueError):
                compare_climatedt_tracks(
                    {"storm": self.reference}, tracks, results_dir=directory,
                    plot_mode="invalid",
                )

    def test_saves_per_strategy_figures(self):
        with tempfile.TemporaryDirectory() as directory:
            figures = {"first": plt.figure(), "second": plt.figure()}
            save_comparison_figures(figures, Path(directory), "metrics")
            self.assertTrue((Path(directory) / "metrics_first.png").is_file())
            self.assertTrue((Path(directory) / "metrics_second.png").is_file())


if __name__ == "__main__":
    unittest.main()
