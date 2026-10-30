"""Execute the whole Source Discovery Tool notebook against a SPHEREx release that is not public yet.

Usage (from the repo root, IRSA network or VPN; needs the `sex` command from SExtractor --
`brew install sextractor` works, no conda env needed):
    python run_sdt_release_test.py <path to spherex_source_discovery_tool_demo.md> \
        --release dr2 [--host https://irsatest.ipac.caltech.edu] [--out sdt_dr2_test_run]

Runs every cell (secs. 4-12, including SExtractor, bokeh and Firefly) on a copy of the notebook's
folder in --out; the tracked notebook is not changed. Test-only changes to the copy:
- SIA -> --host. IRSA must have loaded the release's SIA collection there (check_release.sh).
- s3://nasa-irsa-spherex/<release>/... reads -> the IBE mirror of the same path on --host
  (public S3 is empty before release).
- FireflyClient.make_lab_client() -> make_client() on --host's irsaviewer. A headless run has no
  JupyterLab, and production Firefly can't load pre-release URLs. Browser tabs open; look at them.
At the end it prints a summary: which detectors sec. 6 picked, the SExtractor catalog, the
hard-coded sec. 12.1 source, and its spectrum.
"""
import argparse
import os
import shutil
import sys

import jupytext
import nbformat
from nbclient import NotebookClient

_ap = argparse.ArgumentParser()
_ap.add_argument("notebook")
_ap.add_argument("--release", required=True, help="e.g. dr2")
_ap.add_argument("--host", default="https://irsatest.ipac.caltech.edu")
_ap.add_argument("--out", help="run folder (default: sdt_<release>_test_run in the current folder)")
ARGS = _ap.parse_args()
REL, HOST = ARGS.release, ARGS.host.rstrip("/")
OUT = os.path.abspath(ARGS.out or f"sdt_{REL}_test_run")
NB_NAME = os.path.basename(ARGS.notebook)

if not shutil.which("sex"):
    sys.exit("SExtractor's `sex` command is not on PATH (try `brew install sextractor`).")

OVERRIDE = f'''
# TEST HARNESS ONLY (run_sdt_release_test.py)
Irsa.sia_url = "{HOST}/SIA"
_S3, _IBE = "s3://nasa-irsa-spherex/{REL}/", "{HOST}/ibe/data/spherex/{REL}/"
_real_open = fits.open
def _open_via_ibe(name, *args, **kwargs):
    if isinstance(name, str) and name.startswith(_S3):
        kwargs.pop("fsspec_kwargs", None)
        name = _IBE + name[len(_S3):]
    return _real_open(name, *args, **kwargs)
fits.open = _open_via_ibe
'''

SUMMARY = '''
# TEST HARNESS ONLY: summary
print("=== SUMMARY ===")
print("sec. 6 images:", d2_result["energy_bandpassname"], "and", d4_result["energy_bandpassname"],
      "(expect SPHEREx-D2 and SPHEREx-D4)")
_by_flux = extracted.copy()
_by_flux.sort("FLUX_AUTO", reverse=True)
_nums = list(_by_flux["NUMBER"])
for _sid in selected_ids:
    _rank = _nums.index(_sid) + 1 if _sid in _nums else None
    print(f"SExtractor: {len(extracted)} sources; selected_ids {_sid} is",
          f"#{_rank} by FLUX_AUTO" if _rank else "NOT in the catalog")
_w, _f, _v = np.asarray(wvl, float), np.asarray(flx, float), np.asarray(var, float)
_nb = np.asarray(n_bad, float)
print(f"spectrum: {len(_w)} points, {np.nanmin(_w):.2f}-{np.nanmax(_w):.2f} um,",
      f"median S/N {np.nanmedian(_f / np.sqrt(_v)):.1f},",
      f"{int((_nb > 0).sum())} points with a bad_flags pixel in the aperture")
'''

if os.path.exists(OUT):
    shutil.rmtree(OUT)
shutil.copytree(os.path.dirname(os.path.abspath(ARGS.notebook)), OUT)

nb = jupytext.read(os.path.join(OUT, NB_NAME))
code = [i for i, c in enumerate(nb.cells) if c.cell_type == "code"]
anchor = [i for i in code if "output_notebook(" in nb.cells[i].source]
assert len(anchor) == 1, "can't find the cell to insert the test overrides after"
n_ff = 0
for i in code:
    if "FireflyClient.make_lab_client()" in nb.cells[i].source:
        n_ff += 1
        nb.cells[i].source = nb.cells[i].source.replace(
            "FireflyClient.make_lab_client()",
            f'FireflyClient.make_client(url="{HOST}/irsaviewer", channel_override="sdt_test_{n_ff}")')
print(f"Firefly cells redirected: {n_ff}")
nb.cells.insert(anchor[0] + 1, nbformat.v4.new_code_cell(OVERRIDE))
nb.cells.append(nbformat.v4.new_code_cell(SUMMARY))

print(f"Executing {NB_NAME} in {OUT} (several minutes) ...", flush=True)
NotebookClient(nb, timeout=1800, kernel_name="python3", allow_errors=True,
               resources={"metadata": {"path": OUT}}).execute()
nbformat.write(nb, os.path.join(OUT, NB_NAME.replace(".md", ".ipynb")))

n_err = 0
for c in nb.cells:
    for o in c.get("outputs", []) if c.cell_type == "code" else []:
        if o.get("output_type") == "error":
            n_err += 1
            print("=" * 70, "\n" + c.source[:300], "\n--> ERROR", o["ename"], o["evalue"][:500])
summary = "".join(o.get("text", "") for o in nb.cells[-1].get("outputs", []))
print("\n" + (summary or "(summary cell did not run)"))
print(f"{len(code)} notebook code cells, {n_err} errors. Executed copy: "
      f"{os.path.join(OUT, NB_NAME.replace('.md', '.ipynb'))}")
