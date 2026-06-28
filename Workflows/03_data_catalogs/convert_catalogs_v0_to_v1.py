"""Convert v0 HydroMT data catalogs to the v1 format (branch hydromt_v1_test).

Wraps `hydromt check -d <cat> --format v0 --upgrade`, which performs the bulk of the
v0 -> v1 transformation (path->uri, driver+driver_kwargs->driver{name,options},
rename/unit_add/unit_mult->data_adapter, meta/crs/nodata->metadata). The upgrader is strict
and aborts on the first invalid entry, so this script first sanitizes known issues that exist
in the COMPASS catalogs:

  * invalid metadata URL fields (e.g. ``source_url: none``) -> dropped
  * empty metadata values -> dropped
  * entries hand-edited to v1 (``metadata:`` key) -> renamed back to ``meta:`` so the v0
    upgrader accepts them
  * ``crs``/``nodata``/``attrs`` nested inside ``meta`` -> hoisted to the source top level
    (the upgrader passes them explicitly and also spreads ``meta``, which otherwise collides)
  * Windows backslashes in paths -> forward slashes
  * sources with no usable path/uri (placeholders) and variant-based sources the upgrader
    cannot convert (only ``pcr_globwb``) -> dropped (reported)

It also restores the top-level ``root`` (the upgrader drops it to ``roots: null``) and fixes
the committed local-path root in the CF_forcing catalog.

Usage (from this directory, inside the compass-v1 env):
    pixi run -e compass-v1 python convert_catalogs_v0_to_v1.py

Writes ``<name>_v1___linux.yml`` next to each input and validates it with `hydromt check`.
"""
import subprocess
import pathlib
import sys
import yaml

CAT_DIR = pathlib.Path(__file__).resolve().parent

# in-scope catalogs for the SFINCS + Wflow migration (FIAT deferred). Adjust for Windows variants.
CATALOGS = [
    "datacatalog_general___linux.yml",
    "datacatalog_CF_forcing___linux.yml",
    "datacatalog_SFINCS_coastal_coupling___linux.yml",
    "datacatalog_SFINCS_obspoints___linux.yml",
]

URL_KEYS = {"source_url"}
HOIST_KEYS = ("nodata", "crs", "attrs")


def is_valid_url(v):
    return isinstance(v, str) and v.strip().lower().startswith(("http://", "https://"))


def sanitize(doc):
    n, dropped = 0, []
    for name in list(doc.keys()):
        src = doc[name]
        if name == "meta" or not isinstance(src, dict):
            continue
        if not src.get("path") and not src.get("uri") and not src.get("variants"):
            dropped.append(name); del doc[name]; continue
        if "variants" in src and not src.get("path") and not src.get("uri"):
            dropped.append(name + " [variants]"); del doc[name]; continue
        if "metadata" in src:
            md = src.pop("metadata")
            if isinstance(md, dict):
                src.setdefault("meta", {}).update(md)
            n += 1
        if isinstance(src.get("path"), str) and "\\" in src["path"]:
            src["path"] = src["path"].replace("\\", "/"); n += 1
        m = src.get("meta")
        if isinstance(m, dict):
            for k in HOIST_KEYS:
                if k in m:
                    src.setdefault(k, m.pop(k)); n += 1
            for k in list(m.keys()):
                if k in URL_KEYS and not is_valid_url(m[k]):
                    del m[k]; n += 1
                elif m[k] in (None, "", "none", "None"):
                    del m[k]; n += 1
    return n, dropped


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def main():
    rc = 0
    for cat in CATALOGS:
        src_path = CAT_DIR / cat
        if not src_path.exists():
            print(f"!! {cat} not found, skipping"); continue
        doc = yaml.safe_load(src_path.read_text())
        root = (doc.get("meta") or {}).get("root")
        if root and "compound-flooding-tropical-cyclones" in str(root):
            root = "/p/"  # CF_forcing has a committed local repo path; data uris are absolute
        n, dropped = sanitize(doc)
        tmp_in = CAT_DIR / (cat + ".sanitized.yml")
        tmp_in.write_text(yaml.safe_dump(doc, sort_keys=False))
        print(f"\n=== {cat}  (sanitized {n} fields, dropped {len(dropped)}, root={root}) ===")
        if dropped:
            print("  dropped:", ", ".join(dropped))
        r = run(["hydromt", "check", "-d", str(tmp_in), "--format", "v0", "--upgrade"])
        v1_tmp = tmp_in.with_name(tmp_in.stem + "_v1.yml")
        if not v1_tmp.exists():
            rc = 1
            print("  !! upgrade produced no output:")
            for l in (r.stdout + r.stderr).splitlines():
                if "error" in l.lower() or "Input should" in l:
                    print("    ", l)
            tmp_in.unlink(missing_ok=True)
            continue
        v1doc = yaml.safe_load(v1_tmp.read_text())
        if root and isinstance(v1doc.get("meta"), dict):
            v1doc["meta"]["root"] = root
            v1doc["meta"].pop("roots", None)
        out_path = CAT_DIR / cat.replace("___linux.yml", "_v1___linux.yml")
        out_path.write_text(yaml.safe_dump(v1doc, sort_keys=False))
        tmp_in.unlink(missing_ok=True)
        v1_tmp.unlink(missing_ok=True)
        rv = run(["hydromt", "check", "-d", str(out_path), "--format", "v1"])
        ok = not any("error" in l.lower() for l in (rv.stdout + rv.stderr).splitlines()
                     if "ERROR" in l or "validation error" in l)
        print(f"  -> {out_path.name}: {len(v1doc) - 1} sources, validate v1 {'OK' if ok else 'FAILED'}")
        if not ok:
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
