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
5. Run each notebook against production with `jupytext --to notebook --execute` (SDT needs its
   conda env, `conda-spherex_sdt.yml`).
6. Re-run `techniques-and-tools/MOST_queries.md` once the release is in MOST.
7. Recheck the supplement's flag table number.
8. Mark the draft PR ready for review (it closes the release issue).
9. Update this skill with what was learned.
