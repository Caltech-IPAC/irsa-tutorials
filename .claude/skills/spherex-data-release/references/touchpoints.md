# Release-dependent spots in the SPHEREx notebooks

Values shown are what DR1 used (branch `spherex-dr1`, Nov 2026). Line numbers drift; search for the
quoted text instead. Paths are relative to `tutorials/spherex/`.

## spherex_cutouts.md (TAP + IBE cutouts)

- `IRSA_HOST` test block (`# TESTING ONLY`): sets `IRSA_HOST` and `Irsa.tap_url`. Note
  `Irsa.tap_url` must be set directly; `astroquery.ipac.irsa.conf.tap_url` is read only at import.
- TAP release filter: `AND a.uri LIKE 'ibe/data/spherex/dr1/%'`. `spherex.plane.datarelease` is a
  date, not a release name, so filter on the URI path. Without a filter the query mixes releases
  and `results[...][0]` picks the oldest.
- Cutout URL prefix built from `{IRSA_HOST}/` || `a.uri`.
- Supplement link (`SPHEREx_Expsupp.pdf`).
- `**Updated:**` date at the bottom.

## spherex_psf.md (TAP + IBE cutouts + S3)

- Same `IRSA_HOST` block, release filter and cutout URL prefix as cutouts.
- Sec. 5.1 reads the image's own `EPSF` extension (image version >= 7.0). Sec. 5.2 is a fallback
  for older images that reads ePSF calibration files from S3 (`<rel>/epsf/cal-epsf-v*-*`). It must
  stay guarded so it does not overwrite a new image's own ePSF (a DR1 bug).
- "primary HDU plus 6 extensions (IMAGE, FLAGS, VARIANCE, ZODI, EPSF, WCS-WAVE)" -- check against a
  new image.
- Supplement link, `**Updated:**`.

## spherex_intro.md (SIA)

- Collection list and the `Irsa.query_sia(..., collection='spherex_dr1')` call; bullet list of
  "two collections" (wide `spherex_dr1`, deep `spherex_dr1_deep`). Show only the newest release.
- "This notebook uses data from SPHEREx Data Release 1 (DR1)."
- Sec. 6 extension list and sec. 12 describe the `EPSF` bintable (DR1: 441 rows = 21x21 regions,
  33x33 ePSF, 5x oversampled via `OVSMPX`/`OVSMPY`; columns XCENTER/YCENTER/CWAVE/CBAND).
- Sec. 8d: prints `spectral_image_header['L2 N_*']` and explains it. DR1 had 26 keys: one per flag
  for most flags (STREAK has none), `_NB` = not also flagged bad, `_NS` = not also flagged SOURCE,
  `N_SAT_HALO` = SATELLITE_HALO. Cites the supplement's flag table (Table 16 in supplement v2.0).
- Sec. 9: loops over FLAGS header `MSKNnnnn` (name) / `MSKDnnnn` (description). QR2 used `MP_*`.
  DR1 defined 22 flags.
- Supplement link, `**Updated:**`.

## spherex_source_discovery/spherex_source_discovery_tool_demo.md (SIA + S3 + TAP)

- Sec. 5 parameters: `data_release = "dr1"`, `sia_collection = f"spherex_{data_release}"`,
  `sapm_collection = "cal-sapm-v3-2026-191"`. `data_release` also builds the S3 path
  `nasa-irsa-spherex/{data_release}/solid_angle_pixel_map/{sapm_collection}/...` (folder is
  `solid_angle_pixel_map/`, not `sapm/`). If only the SIA name changes, set `sia_collection`
  directly. Remove any "valid as of <date>" note that no longer holds.
- "Note that we have used the collection `spherex_dr1`" sentence after the SIA query.
- Sec. 7 `get_s3_fpath(cloud_access)`: needs `cloud_access` filled in SIA results.
- Sec. 8 masks: D2 mask from `flags_data2`, D4 mask from `flags_data4` (DR1 fixed a typo here).
- Sec. 12.1 TAP query: `IRSA_HOST` test block (uses `pyvo.dal.TAPService(f"{IRSA_HOST}/TAP")`),
  release filter, cutout URL prefix.
- Sec. 12.2 `bad_flags`: the supplement's CONTAINS_BAD_PIXEL list (see decisions.md). Any listed
  flag in the aperture drops the point from the plotted spectrum (`bokeh_viz.plot_spectrum`).
- Acknowledgements: DOI sentence. Cite only the DOIs of releases the notebook actually loads.
- Supplement link, `**Updated:**`.
- Package `spherex_source_discovery_tool/aperture_photometry.py`:
  `build_phot_table` reads flag names from `MSKN*` (falls back to `MP_*`); `grab_star` skips
  `bad_flags` entries an image does not define. Check both against new FLAGS headers.

## spherex_mosaic.md

- Supplement link. The static `M101_irsa.fits` demo cube was not regenerated for DR1.

## spherex.md (landing page)

- Wording about Quick Releases and Data Releases.

## techniques-and-tools/MOST_queries.md

- No release-specific code (`catalog="spherex"`). Re-run once the release is queryable in MOST to
  refresh its outputs.

## Repo files

- `deploy_to_fornax_manifest.in` and `notebook_metadata.yml` list the SPHEREx files; update them if
  files are added, renamed or removed (a stale manifest entry fails CI).
