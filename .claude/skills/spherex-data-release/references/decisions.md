# Decisions and gotchas from DR1

## Rules (settled; don't re-ask unless something changed)

- **Show only the newest release.** Older collections stay searchable, but the tutorials
  demonstrate only the newest one. TAP notebooks filter on `a.uri LIKE 'ibe/data/spherex/<rel>/%'`.
- **DOIs:** cite only the DOI of each release the notebook actually loads (both if it loads two).
- **`bad_flags` in the Source Discovery Tool:** use the flag list that sets the supplement's
  CONTAINS_BAD_PIXEL warning (supplement v2.0, p. 43). This came from the SSDC developer who wrote
  the SDT. Do not use the shorter "flags for masking pixels for photometry calculations" list on
  the same page. For DR1 the list was: TRANSIENT, OVERFLOW, SUR_ERROR, NONFUNC, DICHROIC,
  MISSING_DATA, HOT, COLD, NONLINEAR, CROSSTALK, GHOST, GHOST_FPA, GHOST_EXT, STREAK, PERSIST,
  BLOOM, SNOWBALL, HALO, SATELLITE_HALO. FULLSAMPLE is informational and was removed.
- **Supplement:** link `SPHEREx_Expsupp.pdf` (the `_QR` suffix was dropped for DR1). For flag
  definitions, the online supplement is the reference -- not plans for a future release. Cite
  the flag table by number (Table 16 in v2.0) and recheck the number when a new version appears.
- **Ownership:** IRSA updates the Source Discovery Tool notebook, not SSDC.
- **Test hosts** go in a `# TESTING ONLY` block with an `IRSA_HOST` variable, so release day is a
  mechanical revert. Push only to a personal fork until release.
- **`cloud_access`:** assume the new SIA collection fills it. If a pre-release collection has it
  empty, build the S3 path from the URI for local testing only, and don't commit that.
- **Demo mosaic cube:** not regenerated per release.

## Gotchas

- irsatest/irsadev are reachable only on the IRSA network or VPN. Cloud notebook servers get 403 or
  cannot resolve them.
- IRSA production TAP/SIA can return 502 for a while. For non-SPHEREx lookups in tests (like
  2MASS), irsatest works too.
- `spherex.plane.datarelease` is a date. Filter by URI path or `obs_publisher_did`.
- Each release has a `_cal` SIA collection listing calibration files (ePSF, SAPM) per detector --
  the quickest way to learn calibration version names.
- IBE mirrors the S3 paths, so `s3://nasa-irsa-spherex/<rel>/...` reads can be pointed at
  `https://<host>/ibe/data/spherex/<rel>/...` for testing.
- Pre-release S3 data sat in a staging bucket (`nasa-irsa-spherex-dr1-staging`) that moved to the
  public `nasa-irsa-spherex` on release day.
- FLAGS keywords changed between releases (`MP_*` in QR2, `MSKN*`/`MSKD*` from QR3). New flags
  can appear in the supplement before they appear in image headers (SATELLITE_HALO, bit 29, was in
  supplement v2.0 but not in DR1 v7.0.5 headers) -- code must tolerate missing flag names.
- To test SIA notebooks against irsadev, set `Irsa.sia_url = "https://irsadev.ipac.caltech.edu/SIA"`
  in a temporary copy (don't commit it). `Irsa.list_collections()` still lists production
  collections, so check irsadev collections with a TAP query on `ivoa.obscore` instead. irsadev SIA
  `access_url`s point at a dev IBE host. IRSA may load only part of the sky, so confirm the notebook
  positions are covered.
- Select detectors by name, never by position. The SDT once took D4 as
  `results.group_by('energy_bandpassname').groups[3]`, which is D4 only if all six detectors have
  images at the position. DR1 at the demo position had D1, D2, D4 and D5, so it silently used D5
  with every label still saying D4. Fixed for DR1 by matching `groups.keys` to `'SPHEREx-D4'`.
- The SDT's sec. 12.1 measures a hard-coded SExtractor source number (`selected_ids = [65]`). The
  numbering changes whenever the images change, so check the number on every release, including
  release day. For DR1 it was left as 65 (a real source, mid-brightness) with a release-day check.
- Production Firefly (irsa.ipac.caltech.edu/irsaviewer) can't load images from pre-release hosts:
  it shows "failed". Use the same host's viewer (`https://irsatest.ipac.caltech.edu/irsaviewer`)
  when testing, and check production Firefly on release day.
- MOST on the pre-release hosts doesn't index the new release, so `MOST_queries.md` can't be
  tested end to end before release. Its SPHEREx cutout cell can be: take the exposures production
  MOST returns, find the new-release file with the same observation ID
  (`level2_<week>_<n>_<det>`) on the pre-release IBE, and run the cell on it.
- Expect overlapping cleanup PRs. Merge them into the release branch before opening the release
  PR, or the same changes show up twice.

## Release-day checklist (DR1's final shape)

1. Confirm the release is live: `scripts/check_release.sh <rel>`; verify the guessed SIA collection
   names and fix intro + SDT if they differ; confirm production SIA rows have `cloud_access`.
2. Revert test hosts: delete each `# TESTING ONLY` block and `IRSA_HOST`, inline
   `https://irsa.ipac.caltech.edu`, and check that
   `grep -rn "irsatest\|irsadev\|IRSA_HOST" tutorials/` finds nothing.
3. Put in the new DOI (SDT acknowledgements).
4. Bump `**Updated:**` in every notebook that changed.
5. Run each notebook against production with `jupytext --to notebook --execute` (SDT needs
   SExtractor: its conda env, `conda-spherex_sdt.yml`, or `brew install sextractor`). Check the
   Firefly tabs show images, and that the SDT's hard-coded `selected_ids` still gives a real
   spectrum.
6. `techniques-and-tools/MOST_queries.md` belongs to other authors and isn't in the release PR.
   For DR1 a separate issue asked them to re-run it once MOST serves the release and check which
   release MOST returns.
7. Recheck the supplement's flag table number.
8. Mark the draft PR ready for review (it closes the release issue).
9. Update this skill with what was learned.
