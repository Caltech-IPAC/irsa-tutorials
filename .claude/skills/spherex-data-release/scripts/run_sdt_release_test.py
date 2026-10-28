"""Run the Source Discovery Tool notebook's own sec. 5, 7, 8 and 12 cells against a SPHEREx
release that is not yet in SIA or public S3.

Usage (from the repo root, IRSA network or VPN):
    python run_sdt_release_test.py <path to spherex_source_discovery_tool_demo.md> \
        --release dr2 [--host https://irsatest.ipac.caltech.edu]

Stand-ins:
- sec. 6 SIA results -> TAP rows on --host, with a `cloud_access` JSON built from the uri
- s3://nasa-irsa-spherex/<release>/... reads -> the IBE mirror of the same path on --host
- secs. 10-11 (SExtractor + interactive selection) -> one 2MASS star near the field center
Secs. 9-11 (SExtractor, Firefly, bokeh) are not exercised.
"""
import argparse
import os
import sys
import json
import traceback

import jupytext
import astropy.io.fits as fits_mod

_ap = argparse.ArgumentParser()
_ap.add_argument("notebook")
_ap.add_argument("--release", required=True, help="e.g. dr2")
_ap.add_argument("--host", default="https://irsatest.ipac.caltech.edu")
ARGS = _ap.parse_args()
REL, HOST = ARGS.release, ARGS.host.rstrip("/")

NB = ARGS.notebook
SDT_DIR = os.path.dirname(os.path.abspath(NB))
sys.path.insert(0, SDT_DIR)

S3_PREFIX = f"s3://nasa-irsa-spherex/{REL}/"
IBE_PREFIX = f"{HOST}/ibe/data/spherex/{REL}/"
_real_open = fits_mod.open


def open_via_ibe(name, *args, **kwargs):
    if isinstance(name, str) and name.startswith(S3_PREFIX):
        kwargs.pop("fsspec_kwargs", None)
        name = IBE_PREFIX + name[len(S3_PREFIX):]
    return _real_open(name, *args, **kwargs)


fits_mod.open = open_via_ibe

cells = [c["source"] for c in jupytext.read(NB).cells if c["cell_type"] == "code"]


def cell(marker):
    hits = [c for c in cells if marker in c]
    assert len(hits) == 1, (marker, len(hits))
    return hits[0]


ns = {"HOST": HOST, "REL": REL}


def run(label, src):
    print(f"\n===== {label} =====", flush=True)
    try:
        exec(src, ns)
        return True
    except Exception:
        traceback.print_exc()
        return False


run("imports", cell("import concurrent.futures"))
run("sec 5 params", cell('data_release = "'))
print("sapm_s3[1] =", ns["sapm_s3"][1])

# ---- sec. 6 stand-in: first D2 and D4 images of the release at the SDT position, from TAP
run("sec 6 stand-in", '''
coord = SkyCoord(150.7817877, 3.3502204, unit='deg')
_svc = pyvo.dal.TAPService(HOST + "/TAP")
_q = """SELECT a.uri, p.energy_bandpassname, p.time_bounds_lower
FROM spherex.artifact a JOIN spherex.plane p ON a.planeid = p.planeid
WHERE 1 = CONTAINS(POINT('ICRS', 150.7817877, 3.3502204), p.poly)
  AND a.uri LIKE 'ibe/data/spherex/%s/%%'
ORDER BY p.time_bounds_lower""" % REL
_t = _svc.search(_q).to_table()
print(len(_t), REL, "images at SDT position;", sorted(set(_t["energy_bandpassname"])))
def _fake_row(det):
    r = _t[_t["energy_bandpassname"] == det][0]
    key = r["uri"][len("ibe/data/spherex/"):]
    return {"cloud_access": json.dumps({"aws": {"bucket_name": "nasa-irsa-spherex", "key": key,
                                                "region": "us-east-1"}})}
d2_result = _fake_row("SPHEREx-D2")
d4_result = _fake_row("SPHEREx-D4")
''')

ok7 = run("sec 7 get_s3_fpath", cell("def get_s3_fpath")) and \
      run("sec 7 load images", cell("s3_fpath_d2 = get_s3_fpath"))
if ok7:
    print("D2 path:", ns["s3_fpath_d2"])
ok8 = ok7 and all([
    run("sec 8 SAPM", cell("s3_fpath_sapm2 = sapm_s3[1]")),
    run("sec 8 units", cell('print(f"Spectral image BUNIT')),
    run("sec 8 background", cell("bkg2 = sep.Background")),
    run("sec 8 masks", cell("mask2 = np.zeros")),
    run("sec 8 reproject", cell("img_data2_reproj, footprint = reproject_exact")),
])
if ok8:
    import numpy as np
    print("D2 BUNIT now:", ns["img_hdr2"]["BUNIT"], "| reprojected finite px:",
          int(np.isfinite(ns["img_data2_reproj"]).sum()),
          "| good-pixel fraction D2 mask:", float(np.nanmean(ns["mask2"] == 1)))

# ---- secs. 10-11 stand-in: one 2MASS star near the field center
run("matched_stars stand-in", '''
import pandas as pd
# 2MASS from --host's TAP (production IRSA TAP can return 502s)
import pyvo as _pyvo
_stars = _pyvo.dal.TAPService(HOST + "/TAP").run_sync(
    f"SELECT ra, dec, k_m FROM fp_psc WHERE CONTAINS(POINT('ICRS', ra, dec), "
    f"CIRCLE('ICRS', {coord.ra.deg}, {coord.dec.deg}, {50/3600})) = 1").to_table()
_stars = _stars[(_stars["k_m"] > 11) & (_stars["k_m"] < 14)]
_stars.sort("k_m")
matched_stars = pd.DataFrame({"NUMBER": [65], "ra": [float(_stars["ra"][0])], "dec": [float(_stars["dec"][0])]})
print(matched_stars, "Kmag", _stars["k_m"][0])
''')

ok12 = all([
    run("sec 12.1 TAP query", cell("service = pyvo.dal.TAPService")),
    run("sec 12.1 process_cutout def", cell("def process_cutout")),
    run("sec 12.1 results_table", cell('results_table["cutout_index"]')),
    run("sec 12.1 download cutouts", cell("Process cutouts in parallel")),
])
os.makedirs("data", exist_ok=True)
if ok12:
    run("sec 12.2 build_phot_table", cell("build_phot_table(results_table"))
    run("sec 12.2 load pickle", cell("with open(phot_path, 'rb')"))
    if run("sec 12.2 grab_star", cell("bad_flags = [")):
        print("n points:", len(ns["wvl"]), "| wl range:", ns["wvl"].min(), ns["wvl"].max())

if "wvl" in ns:
    import numpy as np
    w, f, nb = np.asarray(ns["wvl"]), np.asarray(ns["flx"]), np.asarray(ns["n_bad"])
    k = np.abs(w - 2.16) < 0.1
    print("flux near 2.16 um (mJy, median):", np.nanmedian(f[k]), "| expected for Kmag",
          "~", 666.7e3 * 10 ** (-11.946 / 2.5))
    print("points with >=1 flagged pixel in aperture:", int((nb > 0).sum()), "of", len(nb))
    entry = ns["all_stars_band_phot"]["obs_log"][0][0]
    print("flags recorded:", sorted(entry["dqf_sum"].keys()))
