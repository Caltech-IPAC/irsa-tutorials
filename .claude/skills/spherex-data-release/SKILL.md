---
name: spherex-data-release
description: Update the SPHEREx tutorial notebooks in irsa-tutorials for a new SPHEREx data release (DR2, DR3, a new Quick Release, a reprocessing). Use this whenever someone mentions a new SPHEREx release, new SPHEREx collections or DOIs, a new SPHEREx Explanatory Supplement, new FLAGS bits, or asks to "update the SPHEREx notebooks", "test SPHEREx against irsatest", or prepare the SPHEREx release-day PR -- even if they only ask about one of those pieces.
---

# SPHEREx data-release update

This skill repeats the work done for DR1 (issue #392, released 2026-11-12): find every place in the
SPHEREx notebooks that depends on the release, test against pre-release IRSA hosts, collect the
answers that only other people can give, and land everything in one PR on release day.

Most of the time goes into waiting on other people, not editing code. So the main job is to keep a
written plan current, ask for missing information early, and do everything that does not depend on
those answers right away.

## Files in this skill

- `references/touchpoints.md` -- every release-dependent spot in each notebook, with what it looked
  like for DR1. Read this first; it is the work list.
- `references/decisions.md` -- rules settled during DR1 (which release to show, which DOI to cite,
  which flags count as bad, who owns which notebook) and gotchas that cost time. Read before
  asking the user anything; most questions were already answered once.
- `scripts/check_release.sh <rel>` -- what is live for release `<rel>` (e.g. `dr2`) on the
  pre-release and production hosts, S3 buckets and docs.
- `scripts/run_sdt_release_test.py` -- executes the whole Source Discovery Tool notebook against
  pre-release data and prints a summary of what to check.

## Workflow

### 1. Start a plan file

Create `SPHEREX_<REL>_PLAN.md` at the repo root and add it (and any local test folder) to
`.git/info/exclude` so it never gets committed. Put a short "RESUME HERE" block at the top: branch,
last commit, what is blocked and on whom, next actions. Update it at the end of every work chunk.
The work spans weeks and many sessions, and this file is the only thing that carries over.

Work on a branch named `spherex-<rel>` off `main`. Push it to the user's fork only until release;
the test-host settings must not reach `main` early.

### 2. Find out what exists

Run `scripts/check_release.sh <rel>` (needs IRSA VPN for the pre-release hosts). It reports, for
TAP, SIA, IBE, S3 and the supplement PDF, whether the new release is visible. Record the results as
a "data products needed" table in the plan (product, which notebook uses it, test location,
available yes/no, date checked). Re-run it whenever something might have changed.

### 3. Ask early for what only people can provide

Send these questions as soon as the release date is known. Draft each as a paste-ready message
in a fenced code block. Keep each one short and specific.

- **IRSA** (archive/backend staff): Will the new SIA collections be loaded on the pre-release host
  (irsadev) before release, so the SIA notebooks can be tested? What will the collections be
  named? Will `cloud_access` be filled in those rows? (The Source Discovery Tool reads its S3
  paths from it. Check it in SIA results, e.g. with `scripts/check_release.sh`: TAP's
  `ivoa.obscore` leaves it blank even when SIA fills it.) Will the data be on the staging
  S3 bucket, and when does it move to the public one?
- **IRSA** (DOIs): the DOI for the new release. It is usually not available until close to release.
- **SSDC**: a copy of the new Explanatory Supplement, or confirmation that the online one is
  current. Ask about any new FLAGS bits and which ones should count as bad for photometry.

Then decide with the user what happens if the pre-release SIA load does not come through
(for DR1 the fallback was letting the tutorials go out a few days after release). Set a date to
follow up and write it in the plan.

### 4. Make the changes that do not depend on answers

Go through `references/touchpoints.md` notebook by notebook. For each one, make the change, then
test it against the pre-release host. In order of how early they can be done:

1. TAP notebooks (`spherex_cutouts.md`, `spherex_psf.md`, SDT sec. 12): point at the pre-release
   host through an `IRSA_HOST` variable inside a clearly marked `# TESTING ONLY` block, and add the
   release filter on `a.uri`. Run them end to end with
   `jupytext --to notebook --execute <file>`.
2. Header-driven code and text: open one new-release image and compare its extensions, FLAGS
   `MSKN*` keywords and `L2 N_*` counts with what the notebooks say and parse (see
   touchpoints). New flags or renamed keywords break code silently -- the DR1 photometry recorded
   no flags at all until the parser was fixed.
3. SIA notebooks (`spherex_intro.md`, SDT secs. 6-11): switch to the guessed collection names now
   (pattern below), so only a check is left for release day.
4. Narrative: supplement link and table number, DOIs, release wording on the landing page.

Collection names follow the pattern of the previous release: for DR1 they were guessed as
`spherex_dr1`, `spherex_dr1_deep`, `spherex_dr1_cal` from `spherex_qr3`, `spherex_qr3_deep`,
`spherex_qr3_cal`. Put a "verify the guessed names" step in the release-day checklist with the
exact lines to change if they differ.

Commit in small, single-purpose commits and push to the fork. Ask before committing, and say in
the plan what each commit did.

### 5. Test the Source Discovery Tool

Run the whole notebook, not selected cells. A DR1 test that ran only some cells and stood in a
2MASS star for secs. 9-11 passed while sec. 6 was silently picking a D5 image instead of D4
(see decisions.md). `scripts/run_sdt_release_test.py` executes every cell on a copy, with the
pre-release host's SIA, its IBE mirror in place of S3, and its Firefly viewer. It needs the
release's SIA collection on that host and SExtractor's `sex` command (`brew install sextractor`
is enough; the notebook's conda env is not needed). Run it from the repo root:

```bash
.venv/bin/python .claude/skills/spherex-data-release/scripts/run_sdt_release_test.py \
  tutorials/spherex/spherex_source_discovery/spherex_source_discovery_tool_demo.md --release dr2 \
  --out <your excluded test folder>/sdt_run
```

Check its summary: sec. 6 must pick D2 and D4, the hard-coded `selected_ids` source must exist,
and its spectrum should look like a real source. Ask the user to look at the two Firefly tabs.
For DR1: 33 notebook code cells, 0 errors; 79 SExtractor sources, #65 was 40th by flux; 276
points, 0.75-4.99 um, median S/N about 11, 50 points with a `bad_flags` pixel in the aperture.

### 6. Release day

Keep a release-day checklist in the plan from the start and refine it as answers arrive. Its
final shape for DR1 is in `references/decisions.md` under "Release-day checklist". Open a draft PR
well before release so it can be reviewed, with the checklist in the PR description. Merge any
pending cleanup PRs into the release branch first so the diff shows only release work.

As the last checklist item, update this skill with what was learned (see below).

## Keeping this skill current

After each release, update `references/touchpoints.md` with the new values and line numbers, add
new decisions and gotchas to `references/decisions.md`, and fix anything in this file that was
wrong. Leave out people's names: say which group to ask (IRSA, SSDC), not who.
