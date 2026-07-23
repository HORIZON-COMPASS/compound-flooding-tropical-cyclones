"""
Snakemake lint tests — syntax validation for all .smk files.

`snakemake --lint` checks for syntax errors, deprecated features, and common
best-practice violations without requiring input data or a resolved DAG.

Requires snakemake on PATH. Run with:
    pixi run -e compass-v1 pytest tests/snakemake/ -m snakemake -v

Or within any environment that has snakemake installed:
    pytest tests/snakemake/ -m snakemake -v
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SMK_DIR     = Path(__file__).parents[2] / "Workflows" / "02_workflow_rules"
SMK_FILES   = sorted(SMK_DIR.glob("*.smk"))
LINT_CONFIG = Path(__file__).parent / "config_lint.yml"


def _snakemake_binary() -> str | None:
    """Locate the snakemake binary: check PATH first, then the current Python's bin dir."""
    found = shutil.which("snakemake")
    if found:
        return found
    # When running inside a pixi env, Python lives in .pixi/envs/<name>/bin/ alongside snakemake
    candidate = Path(sys.executable).parent / "snakemake"
    return str(candidate) if candidate.exists() else None


@pytest.mark.snakemake
@pytest.mark.parametrize("smk", SMK_FILES, ids=lambda p: p.name)
def test_snakemake_lint(smk: Path):
    smk_bin = _snakemake_binary()
    if smk_bin is None:
        pytest.skip("snakemake not found — activate a pixi env that includes it")

    result = subprocess.run(
        [
            smk_bin, "--lint",
            "--snakefile", str(smk),
            "--configfile", str(LINT_CONFIG),
        ],
        capture_output=True,
        text=True,
        cwd=str(smk.parents[1]),  # run from Workflows/ so include: paths resolve
    )
    # snakemake --lint exits non-zero for both parse errors AND style warnings.
    # Only fail the test on actual parse/runtime errors (Python Traceback in stderr).
    # Style warnings (no log directive, absolute paths, no conda env) are expected
    # in a scientific workflow and are not treated as failures here.
    if result.returncode != 0 and "Traceback" in result.stderr:
        pytest.fail(
            f"Parse/runtime error in {smk.name}\n"
            f"--- stderr ---\n{result.stderr}"
        )
    # Emit style warnings as pytest warnings so they're visible but non-blocking
    if result.returncode != 0 and result.stderr.strip():
        import warnings
        warnings.warn(f"Lint style warnings in {smk.name} (not a failure):\n{result.stderr[:500]}")


def test_smk_files_exist():
    """Sanity check: the workflow rules directory contains at least 10 .smk files."""
    assert len(SMK_FILES) >= 10, (
        f"Expected at least 10 .smk files in {SMK_DIR}, found {len(SMK_FILES)}"
    )


def test_smk_files_are_readable():
    """All .smk files can be opened and read (not empty, no permissions issues)."""
    for smk in SMK_FILES:
        content = smk.read_text(encoding="utf-8")
        assert len(content) > 0, f"{smk.name} is empty"
