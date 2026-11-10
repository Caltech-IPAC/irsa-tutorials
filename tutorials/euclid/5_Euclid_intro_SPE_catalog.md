---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  name: python3
  display_name: python3
  language: python
authors:
  - name: Tiffany Meshkat
  - name: Anahita Alavi
  - name: Anastasia Laity
  - name: Andreas Faisst
  - name: Brigitta Sipőcz
  - name: Dan Masters
  - name: Harry Teplitz
  - name: Jaladh Singhal
  - name: Keunho Kim
  - name: Shoubaneh Hemmati
  - name: Vandana Desai
  - name: Troy Raen
  - name: Jessica Krick
---

# Euclid DR1: SPE Catalogs

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- know which Euclid DR1 SPE catalog tables are available at IRSA, and how to view their columns;
- know how to join the SPE redshift, quality, and lines tables with the MER and PHZ catalogs in ADQL;
- know how to filter the SPE tables on their redshift, quality, and line columns, with the example of selecting galaxies with a strong Hα line, and see how many galaxies each criterion removes;
- understand the line quality flags, and how to pick the combined-grism line measurements;
- know how to find and plot the spectra of galaxies, and mark their measured lines;
- know the limitations of the DR1 line measurements for line ratios.

+++

## Introduction

+++
Euclid launched in July 2023 as a European Space Agency (ESA) mission with contributions from NASA.
The primary science goals of Euclid are to better understand the composition and evolution of the dark Universe.
The Euclid mission provides space-based imaging and spectroscopy, complemented by ground-based imaging, to achieve these goals.
These data will be archived by multiple global repositories, including IRSA, where they will support a wide range of astrophysical research.

Euclid Data Release 1 (DR1) will include imaging, spectroscopy, and catalogs covering approximately 1,900 square degrees of the Euclid Wide Survey, along with observations of the Euclid Deep Fields. Further information about Euclid DR1 can be found at: https://euclid.caltech.edu/page/euclid-dr1

Among the DR1 data products are multiple spectroscopic catalogs produced by the Spectroscopic Analysis (SPE) Processing Function, including spectroscopic redshifts, emission-line measurements, and other derived spectral properties.

This notebook provides an introduction to the Euclid DR1 SPE catalogs.
If you have questions about this notebook, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## Imports

```{important}
We rely on ``astroquery`` features that have been recently added, so please make sure you have version v0.4.10 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed
# %pip install matplotlib astropy 'astroquery>=0.4.11'
```

```{code-cell} ipython3
import warnings

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from astropy.table import QTable, Table, join, vstack
from astropy.utils.data import conf
from astropy.visualization import quantity_support
from astropy import units as u

from astroquery.ipac.irsa import Irsa

# Suppress warnings about deprecated units in the spectrum files, and about the astropy cache directory
warnings.filterwarnings('ignore', message="The unit 'Angstrom' has been deprecated", category=u.UnitsWarning)
warnings.filterwarnings('ignore', message="The unit 'erg' has been deprecated", category=u.UnitsWarning)
warnings.filterwarnings('ignore', message='XDG_CACHE_HOME is set')

# Let matplotlib pick up the axis units from the quantities we plot
quantity_support()

# The Euclid spectrum files can take a while to read,
# which can exceed astropy's default timeout limit. Increase it.
conf.remote_timeout = 120
```

<!-- TODO: remove this cell when DR1 is public on the production IRSA services. -->
```{code-cell} ipython3
# DR1 is not public yet, so we query IRSA's development server.
Irsa.tap_url = 'https://irsadev.ipac.caltech.edu/TAP'
```

## 1. Find the SPE catalog tables at IRSA

The SPE catalogs are large, so we never download them: we query them on IRSA's server with ADQL through the Table Access Protocol (TAP) service, and get back only the rows we ask for.
In this section we only look at which tables are available and what their columns are.

Like the MER catalog (see the [MER catalogs notebook](2_Euclid_intro_MER_catalog.md)), the DR1 SPE catalogs are split by survey, with a `_wide_survey` and a `_deep` version of each table:
* `euclid_dr1_spectro_z_spe_galaxy_candidates_*` -- the spectroscopic redshift candidates of each source fitted with galaxy templates, ranked by `spe_rank` (0 is the best candidate), with their probability `spe_z_prob` and reliability `spe_z_rel`.
* `euclid_dr1_spectro_line_features_spe_line_features_*` -- the line measurements for each redshift candidate (matched on `object_id` and `spe_rank`): the lines are measured at their expected wavelengths at that redshift, with their wavelength, flux, signal-to-noise ratio, and quality, from a Gaussian fit (`_gf`) and from direct integration (`_di`). The 11 main lines (Lyα, [OII], Hβ, [OIII], Hα, [NII]6585, [SII]6718, [SIII]9068, [SIII]9530, Paγ, and Paβ) are measured for every candidate, whether or not they are detected, so many rows are non-detections (for example `spe_line_snr_gf = -99`); other lines are only stored when they are detected. A line can have several rows, one for each grism, with different `spe_line_qual_gf` values.
* `euclid_dr1_spectro_z_spe_quality_*` -- one row per source processed by SPE, including sources without a redshift candidate, describing its spectrum and the processing: grism(s), wavelength range, number of valid pixels, number of dithers, and error and warning flags.
* `euclid_dr1_spectro_z_spe_classification_*` -- the classification of each source as a galaxy, star, or QSO from its spectrum (not used in this notebook).

The SPE tables have no coordinates; they are joined to the MER catalog of the same survey on `object_id`.
In this notebook we work in the Euclid Deep Field North, so we use the deep tables.

The DR1 catalogs at IRSA can be listed with the `list_catalogs` method:

```{note}
In the example below, we use ``include_metadata_tables=True``, as IRSA lists the PHZ and SPE tables among the metadata tables rather than the catalogs.
This parameter was introduced in astroquery 0.4.11. For earlier versions, please omit this parameter to receive the combined list of available catalogs and metadata tables.
```

```{code-cell} ipython3
Irsa.list_catalogs(filter='euclid_dr1', include_metadata_tables=True)
```

We will use the following tables, including the MER catalog for the coordinates and morphology of the sources, and the PHZ catalog for their photometric redshifts:

```{code-cell} ipython3
table_mer = 'euclid_dr1_mer_catalogue_deep_survey'
table_galaxy_candidates = 'euclid_dr1_spectro_z_spe_galaxy_candidates_deep'
table_quality = 'euclid_dr1_spectro_z_spe_quality_deep'
table_lines = 'euclid_dr1_spectro_line_features_spe_line_features_deep'
table_photo_z = 'euclid_dr1_phz_photo_z_deep_survey'
```

### Explore the columns of the SPE tables

Before writing a query, it helps to know which columns each table has and what they mean.
The `list_columns` method returns the name and description of every column of a table, without querying any data.
We look at the two tables we use most below: the lines table, where each line measurement has a set of columns from the Gaussian fit (ending in `_gf`) and the same quantities from direct integration (ending in `_di`), and the redshift table, where each row is one redshift candidate of a source:

```{code-cell} ipython3
columns_lines = Irsa.list_columns(catalog=table_lines)
columns_galaxy_candidates = Irsa.list_columns(catalog=table_galaxy_candidates)
print(f'{len(columns_lines)} columns in the lines table, {len(columns_galaxy_candidates)} in the redshift table')
```

```{tip}
A few highlights of the columns we use below:
- `object_id`, the MER object ID used to join the SPE tables with each other and with the MER and PHZ catalogs;
- `spe_rank`, the rank of the redshift candidate, 0 for the best one;
- `spe_z`, `spe_z_prob`, and `spe_z_rel`, the spectroscopic redshift, its probability, and its reliability;
- `spe_line_name`, the name of the measured line, for example `'Halpha'`;
- `spe_line_central_wl_gf`, `spe_line_flux_gf`, and `spe_line_snr_gf`, the central wavelength (in Angstrom), flux (in erg s<sup>-1</sup> cm<sup>-2</sup>), and signal-to-noise ratio of the line from the Gaussian fit;
- `spe_line_qual_gf`, the quality of the Gaussian fit measurements;
- we use the Gaussian fit measurements throughout: in DR1, the direct integration fluxes (`_di`) are about 15% larger than the Gaussian fit ones, because of a too aggressive continuum subtraction (Le Brun et al. 2026);
- `spe_npix` and `spe_n_dith_med` in the quality table, the number of valid pixels and the median number of dithers of the spectrum.
```

```{code-cell} ipython3
# Full list of columns of the lines table and their description
columns_lines
```

```{code-cell} ipython3
# Full list of columns of the redshift table and their description
columns_galaxy_candidates
```

## 2. Find objects within a specified sky region defined by RA and Dec ranges
Euclid DR1 covers a large sky area of approximately 1,900 square degrees. To reduce query execution time and avoid potential timeout issues, we recommend restricting searches to smaller sky regions.

As an example, let's work in a 0.1 square-degree region of the Euclid Deep Field North (EDF-N), and first count the MER sources in it:

```{code-cell} ipython3
ra_min, ra_max = 269.34391, 270.12192
dec_min, dec_max = 65.85958, 66.17581

query_mer = f"""
SELECT COUNT(*) AS n_sources
FROM {table_mer}
WHERE ra BETWEEN {ra_min} AND {ra_max}
  AND dec BETWEEN {dec_min} AND {dec_max}
"""
n_mer = Irsa.query_tap(query_mer).to_table()['n_sources'][0]
print(f'Found {n_mer} MER sources in the region')
```

Even this small region contains tens of thousands of MER sources.
Only a fraction of them have an SPE redshift, and fewer still have a redshift and line measurements we can trust; finding those is the subject of the next section.

## 3. Example: filter the SPE tables for galaxies with a strong H-alpha line

The SPE tables can be filtered on any of their columns, and the right filters depend on the science case.
In this section, we show how to combine the SPE redshift, quality, and lines tables with the MER and PHZ catalogs, using one common case as an example: selecting galaxies with a strong, well-measured Hα line, the main target of the Euclid spectroscopic survey.

The criteria below are one recipe for this case, adapted from the ones used by Le Brun et al. (2026, Sect. 4.2) to build a reference sample of Hα emitters in the EDF-N, which agrees with DESI redshifts for about 99% of the galaxies.
This is the redshift range and the line that the Euclid spectroscopic survey targets, and for which the SPE galaxy model is optimized.
The recipe is deliberately strict, and is not a definition of a "good" spectrum: a galaxy at z > 1.85, for example, fails criterion 2, but can have an excellent spectrum and redshift from its [OIII] or [OII] lines.
For another science case, drop or change the criteria that do not apply, for example the redshift range and the line in criteria 2 and 6.

| | Criterion | Columns and thresholds | Why |
|---|---|---|---|
| 1 | High-confidence redshift | `spe_rank = 0`, `spe_z_prob > 0.99`, `spe_z_rel > 0.05` | The best redshift candidate has a high probability (the fraction of the redshift probability distribution in its peak) and is not flagged as unreliable |
| 2 | Hα in the red grism spectrum | `0.85 < spe_z < 1.85` | Hα falls in the red grism spectrum (1.21 to 1.89 micron), away from its edges |
| 3 | Both grisms, full wavelength coverage and depth | `spe_grism = 'BGS+RGS'`, `spe_npix >= 800`, `spe_n_dith_med >= 15` (quality table) | The source has a blue and a red grism spectrum, with almost none of the wavelength range masked and close to the full depth of the survey |
| 4 | Extended, non-spurious source | `point_like_prob < 0.05`, `spurious_prob < 0.5`, `kron_radius < 50` (MER catalog) | Removes stars (and also compact galaxies and quasars), spurious detections, and very large galaxies, which are most likely nearby |
| 5 | Consistent with the photometric redshift | `phz_95_int_1 < spe_z < phz_95_int_2` (PHZ catalog) | A check against a wrong line identification |
| 6 | Strong, well-measured Hα line | Gaussian fit good in at least one grism; combined signal-to-noise ratio > 10 and flux > 1.5×10<sup>-16</sup> erg s<sup>-1</sup> cm<sup>-2</sup> (lines table, see section 3.2) | The line is clearly detected, as in the visually inspected sample of the paper |

```{note}
How the criteria compare with Le Brun et al. (2026), and when to change them:
- Criterion 1: `spe_z_prob > 0.99` is the threshold of Sect. 4.1 of the paper, while Sect. 4.2 uses 0.94. The reliability `spe_z_rel` comes from a machine-learning classifier, which in DR1 is trained on simulations, so only a low threshold is recommended.
- Criterion 2: at z < 1.08, Hα is also in the blue grism spectrum (up to 1.37 micron).
- Criterion 3: the combined blue and red grism spectrum has at most about 850 valid pixels, and most sources in this region have about 20 dithered exposures. Drop the `BGS+RGS` requirement for the Euclid Wide Survey, which only has red grism spectra.
- Criterion 4: these are the criteria of the Euclid Galaxy Clustering Science Working Group.
- Criterion 5: the paper finds that with this criterion, 98% of the redshifts are correct at 0.84 < z < 1.6 in the EDF-N.
```

### 3.1 Get the redshifts, spectrum quality, and morphology of the sources

Criteria 1 to 5 only involve one row per source, from the redshift, quality, MER, and PHZ tables, so we get all of them with a single query, which joins the tables on IRSA's server.
We use a left outer join for the PHZ catalog, to keep the sources without a photometric redshift:

```{code-cell} ipython3
query_candidates = f"""
SELECT mer.object_id, mer.ra, mer.dec,
       mer.point_like_prob, mer.spurious_prob, mer.kron_radius,
       spe.spe_z, spe.spe_z_prob, spe.spe_z_rel,
       quality.spe_grism, quality.spe_npix, quality.spe_n_dith_med,
       phz.phz_median, phz.phz_95_int_1, phz.phz_95_int_2
FROM {table_mer} AS mer
JOIN {table_galaxy_candidates} AS spe ON spe.object_id = mer.object_id
JOIN {table_quality} AS quality ON quality.object_id = mer.object_id
LEFT OUTER JOIN {table_photo_z} AS phz ON phz.object_id = mer.object_id
WHERE mer.ra BETWEEN {ra_min} AND {ra_max}
  AND mer.dec BETWEEN {dec_min} AND {dec_max}
  AND spe.spe_rank = 0
"""
candidates = Irsa.query_tap(query_candidates, async_job=True).to_table()
print(f'{len(candidates)} of the {n_mer} MER sources have an SPE redshift')
```

We now apply the criteria one at a time and count how many galaxies are left after each of them.
Missing values (for example a missing photometric redshift) are masked in the table; we fill them with NaN, so that they fail every comparison.

The helper functions of this notebook, like the two below, are in collapsed cells: click on a cell to expand it.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def as_float(table, column):
    """Get a column of a table as floats, with NaN for the missing values.

    Parameters
    ----------
    table : astropy.table.Table
        Table to read the column from.
    column : str
        Name of the column.

    Returns
    -------
    numpy.ndarray
        Values of the column as floats, with the masked (missing) values set to NaN.
    """
    return np.ma.filled(np.ma.asarray(table[column], dtype=float), np.nan)


def cut_flow(cuts):
    """Count the galaxies left after applying each criterion in turn.

    Parameters
    ----------
    cuts : dict
        Boolean array of each criterion (True for the galaxies that pass it), keyed by the name of the criterion.

    Returns
    -------
    astropy.table.Table
        One row per criterion, with the number of galaxies that pass it and all of the criteria before it.
    """
    remaining = np.ones(len(next(iter(cuts.values()))), dtype=bool)
    steps, counts = ['SPE redshift'], [len(remaining)]
    for name, cut in cuts.items():
        remaining &= cut
        steps.append(name)
        counts.append(remaining.sum())
    return Table({'criterion': steps, 'n_galaxies': counts})
```

```{code-cell} ipython3
def criteria(table):
    """Apply criteria 1 to 5 to a table of sources.

    Parameters
    ----------
    table : astropy.table.Table
        Sources, with the columns of the redshift, quality, MER, and PHZ tables used by the criteria.

    Returns
    -------
    dict
        Boolean array of each criterion (True for the sources that pass it), keyed by the name of the criterion.
    """
    z = as_float(table, 'spe_z')
    return {
        '1. High-confidence redshift': (as_float(table, 'spe_z_prob') > 0.99) & (as_float(table, 'spe_z_rel') > 0.05),
        '2. Hα in the red grism (0.85 < z < 1.85)': (z > 0.85) & (z < 1.85),
        '3. Both grisms, full coverage and depth': ((table['spe_grism'] == 'BGS+RGS')
                                                    & (as_float(table, 'spe_npix') >= 800)
                                                    & (as_float(table, 'spe_n_dith_med') >= 15)),
        '4. Extended, non-spurious source': ((as_float(table, 'point_like_prob') < 0.05)
                                             & (as_float(table, 'spurious_prob') < 0.5)
                                             & (as_float(table, 'kron_radius') < 50)),
        '5. Consistent with the photo-z': (z > as_float(table, 'phz_95_int_1')) & (z < as_float(table, 'phz_95_int_2')),
    }
```

```{code-cell} ipython3
cut_flow(criteria(candidates))
```

### 3.2 Add the H-alpha measurements

As we saw in section 1, the lines table has several rows for each line of a source: one for each grism in which the line is measured.
The quality flag `spe_line_qual_gf` has three digits, described in the [DR1 documentation of the SPE catalogs](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmdr1/spedpd/dpcards/spe_spepfoutputcatalog.html): the first digit is the measurement method (2 for the Gaussian fit), the last digit tells which spectrum the measurement comes from, and the middle digit tells how good it is:

| Digit | Value | Meaning |
|---|---|---|
| last (grism) | 1, 2, 3 | blue grism, red grism, blue and red grisms combined |
| middle (quality) | 0 | not measured |
| | 1 | good measurement |
| | 2 | low signal-to-noise ratio (below 2) |
| | 3 | positive equivalent width (probably a negative continuum) |
| | 4 | combined measurement only: the blue and red grism measurements are not both good |
| | 5 | acceptable, but up to 15% of the pixels of the line are missing |
| | 7 | combined measurement only: error when combining the blue and red grism measurements |
| | 8 | unreliable, more than 15% of the pixels of the line are missing |
| | 9 | the line is outside of the wavelength range (with a negative sign, for example -92) |

For each galaxy, we get all of the Hα measurements in the wavelength range (`spe_line_qual_gf > 0`), and use them in two ways:
- the quality of the single-grism measurements (last digit 1 or 2) tells us whether Hα is well measured in at least one of the grisms: we require a quality digit of 1 (good) or 5 (acceptable), that is a quality flag of 211, 212, 251, or 252. Hα does not need to be well measured in both grisms: at z > 1.08 it is only in the red grism, and near the start of the red grism spectrum it is often only well measured in the blue grism;
- the combined measurement (last digit 3), which uses both grisms where the line falls in both, gives the flux and signal-to-noise ratio of the line.

```{note}
The combined measurement of most lines has a quality digit of 4 (flag 243), because the line is not well measured in both grisms, most often because it is outside of the wavelength range of one of them; the combined measurement then has the same flux and signal-to-noise ratio as the single-grism one.
Only lines well measured in both grisms have a combined flag of 213.
This is why we take the quality from the single-grism measurements, and only require Hα to be well measured in one of them.
```

The IRSA TAP service accepts at most 1000 values in an `IN (...)` list, so we query the Hα measurements in batches of 1000 object IDs.
We only query the galaxies that pass the first two criteria (the others are rejected anyway), which still includes the ones failing the later criteria, so that we can look at some of them in section 5:

```{code-cell} ipython3
def query_halpha(object_ids, batch_size=1000):
    """Get the Hα measurements of a list of sources from the lines table.

    Parameters
    ----------
    object_ids : list of int
        MER object IDs of the sources.
    batch_size : int, optional
        Number of object IDs in each query (the IRSA TAP service accepts at most 1000).

    Returns
    -------
    astropy.table.Table
        One row per Hα measurement in the wavelength range, for the best (rank 0) redshift of each source:
        one for each grism the line falls in, and one for the two grisms combined.
    """
    tables = []
    for start in range(0, len(object_ids), batch_size):
        ids = ', '.join(str(object_id) for object_id in object_ids[start:start + batch_size])
        query = f"""
            SELECT object_id, spe_line_qual_gf, spe_line_flux_gf, spe_line_snr_gf, spe_line_central_wl_gf
            FROM {table_lines}
            WHERE object_id IN ({ids})
              AND spe_rank = 0
              AND spe_line_name = 'Halpha'
              AND spe_line_qual_gf > 0
        """
        tables.append(Irsa.query_tap(query, async_job=True).to_table())
    return vstack(tables)
```

```{code-cell} ipython3
cuts = criteria(candidates)
in_range = cuts['1. High-confidence redshift'] & cuts['2. Hα in the red grism (0.85 < z < 1.85)']
halpha_rows = query_halpha(list(candidates['object_id'][in_range]))

grism_digit = halpha_rows['spe_line_qual_gf'] % 10
quality_digit = (halpha_rows['spe_line_qual_gf'] // 10) % 10
good_single_grism = halpha_rows[(grism_digit < 3) & np.isin(quality_digit, [1, 5])]

# Keep the combined measurement, and flag whether Hα is good in at least one grism
halpha = halpha_rows[grism_digit == 3]
halpha['halpha_good_in_a_grism'] = np.isin(halpha['object_id'], good_single_grism['object_id'])
for column in ['qual_gf', 'flux_gf', 'snr_gf', 'central_wl_gf']:
    halpha.rename_column(f'spe_line_{column}', f'halpha_{column}')

galaxies = join(candidates, halpha, keys='object_id', join_type='left')
print(f"Got the Hα measurements of {len(halpha)} of the {in_range.sum()} galaxies; "
      f"Hα is good in at least one grism for {halpha['halpha_good_in_a_grism'].sum()} of them")
```

We add the Hα criterion to the others, and look at the number of galaxies left after each of them:

```{code-cell} ipython3
def all_criteria(table):
    """Apply criteria 1 to 6 to a table of sources.

    Parameters
    ----------
    table : astropy.table.Table
        Sources, with the columns used by `criteria` and their Hα measurements.

    Returns
    -------
    dict
        Boolean array of each criterion (True for the sources that pass it), keyed by the name of the criterion.
    """
    cuts = criteria(table)
    cuts['6. Strong, well-measured Hα'] = (np.ma.filled(table['halpha_good_in_a_grism'], False)
                                           & (as_float(table, 'halpha_snr_gf') > 10)
                                           & (as_float(table, 'halpha_flux_gf') > 1.5e-16))
    return cuts
```

```{code-cell} ipython3
cuts = all_criteria(galaxies)
passes_all = np.logical_and.reduce(list(cuts.values()))
cut_flow(cuts)
```

Out of the tens of thousands of MER sources in the region, about 130 galaxies pass all of the criteria.
The Hα criterion removes most of the galaxies left after the first five criteria, mainly because of the flux cut: most of the Hα lines are fainter than 1.5×10<sup>-16</sup> erg s<sup>-1</sup> cm<sup>-2</sup>, close to the flux limit of the Euclid Wide Survey (2×10<sup>-16</sup> erg s<sup>-1</sup> cm<sup>-2</sup>), which the deep fields go below (see section 4).
Let's sort the galaxies that pass all of the criteria by the signal-to-noise ratio of their Hα line:

```{code-cell} ipython3
strong_halpha = galaxies[passes_all]
if len(strong_halpha) == 0:
    raise ValueError('No galaxy passes all of the criteria; try a larger region or relax some of the criteria.')
strong_halpha.sort('halpha_snr_gf', reverse=True)
strong_halpha['object_id', 'spe_z', 'phz_median', 'halpha_flux_gf', 'halpha_snr_gf'][:10]
```

## 4. Redshifts and Hα fluxes

Let's compare the redshifts of the galaxies that pass all of the criteria with those of all the galaxies with a high-confidence redshift in the same range, and look at their Hα fluxes as a function of redshift:

```{code-cell} ipython3
high_confidence = cuts['1. High-confidence redshift'] & cuts['2. Hα in the red grism (0.85 < z < 1.85)']
has_halpha = as_float(galaxies, 'halpha_flux_gf') > 0

fig, (ax_z, ax_flux) = plt.subplots(1, 2, figsize=(13, 4.5))

bins = np.arange(0.85, 1.86, 0.05)
ax_z.hist(galaxies['spe_z'][high_confidence], bins=bins, color='lightgrey', label='High-confidence redshift')
ax_z.hist(galaxies['spe_z'][passes_all], bins=bins, color='tab:red', label='All criteria')
ax_z.set_xlabel('Spectroscopic redshift')
ax_z.set_ylabel('Number of galaxies')
ax_z.legend()

ax_flux.scatter(galaxies['spe_z'][high_confidence & has_halpha], as_float(galaxies, 'halpha_flux_gf')[high_confidence & has_halpha],
                s=5, color='lightgrey', label='High-confidence redshift')
points = ax_flux.scatter(galaxies['spe_z'][passes_all], as_float(galaxies, 'halpha_flux_gf')[passes_all],
                         s=12, c=as_float(galaxies, 'halpha_snr_gf')[passes_all], cmap='viridis', vmin=10, vmax=50,
                         label='All criteria')
ax_flux.axhline(1.5e-16, color='grey', linestyle='--', linewidth=0.8)
ax_flux.set_yscale('log')
ax_flux.set_ylim(1e-17, 3e-15)
ax_flux.set_xlabel('Spectroscopic redshift')
ax_flux.set_ylabel('Hα flux (erg s$^{-1}$ cm$^{-2}$)')
ax_flux.legend(loc='upper right')
fig.colorbar(points, ax=ax_flux, label='Hα signal-to-noise ratio')
plt.show()
```

## 5. Look at the spectra

To see what the criteria mean, let's look at the spectra of some of the galaxies that pass all of them, and of some galaxies that fail one of them.
We find and read in the spectra with the functions of the [SIR 1D spectra notebook](3_Euclid_intro_1D_spectra.md), which explains them step by step:
the table `euclid.objectid_spectrafile_association_dr1` lists the file and HDU of each spectrum of a source, and IRSA's spectrumDM service returns a given HDU as a VOTable, at a URL built from the file path and the HDU number.
A galaxy in the Euclid Deep Fields can have spectra from both the red and the blue grism, and several spectra of the same grism, so we keep the spectrum with the most usable flux bins for each grism.
If a spectrum cannot be read (for example, because of a network error or a timeout), we skip it and use the other spectra of the galaxy.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
SPECTRUMDM_URL = 'https://irsacloud.ipac.caltech.edu/api/spectrumdm/convert/'
SPECTRUM_DATASET_ID = 'euclid_combspec_dr1'


def spectrum_url(path, hdu):
    """Build the spectrumDM URL of one spectrum.

    Parameters
    ----------
    path : str
        Path of the spectrum file, from the ``path`` column of the association table.
    hdu : int
        Number of the HDU of the spectrum in the file.

    Returns
    -------
    str
        URL that returns the spectrum as a VOTable.
    """
    return f'{SPECTRUMDM_URL}{path}?dataset_id={SPECTRUM_DATASET_ID}&hdu={hdu}'


def find_spectra(object_id):
    """Find the spectra of an object in the spectrum association table.

    Parameters
    ----------
    object_id : int
        MER object ID of the source.

    Returns
    -------
    astropy.table.Table
        One row per spectrum, with the file path and HDU number of the spectrum, and its spectrumDM URL.
    """
    query = f"""
        SELECT objectid, tileid, path, hdu
        FROM euclid.objectid_spectrafile_association_dr1
        WHERE objectid = {object_id}
    """
    spectrum_files = Irsa.query_tap(query).to_table()
    spectrum_files['url'] = [spectrum_url(path, hdu) for path, hdu in zip(spectrum_files['path'], spectrum_files['hdu'])]
    return spectrum_files


def bad_flux_bins(spectrum):
    """Find the flux bins of a spectrum that should not be used.

    Parameters
    ----------
    spectrum : astropy.table.QTable
        Spectrum, with a MASK column.

    Returns
    -------
    numpy.ndarray
        Boolean array, True for the flux bins with an odd MASK value or a MASK value of 64 or more.
    """
    mask = spectrum['MASK'].value
    return (mask % 2 == 1) | (mask >= 64)


def grism_name(spectrum):
    """Tell which kind of grism a spectrum comes from.

    Parameters
    ----------
    spectrum : astropy.table.QTable
        Spectrum, with a WAVELENGTH column.

    Returns
    -------
    str
        ``'Blue grism'`` if the spectrum starts below 1.1 micron, ``'Red grism'`` otherwise.
    """
    return 'Blue grism' if spectrum['WAVELENGTH'].min() < 1.1 * u.micron else 'Red grism'


def best_spectra(spectra):
    """Pick the spectrum with the most usable flux bins for each kind of grism.

    Parameters
    ----------
    spectra : list of astropy.table.QTable
        Spectra of one object.

    Returns
    -------
    dict
        Best spectrum of each kind of grism, keyed by the name of the grism.
    """
    best = {}
    for spectrum in sorted(spectra, key=lambda spectrum: np.sum(~bad_flux_bins(spectrum)), reverse=True):
        best.setdefault(grism_name(spectrum), spectrum)
    return best


def get_best_spectra(object_id):
    """Read the spectra of an object and pick the best one of each kind of grism.

    Spectra that cannot be read are skipped, with a message.

    Parameters
    ----------
    object_id : int
        MER object ID of the source.

    Returns
    -------
    dict
        Best spectrum of each kind of grism, keyed by the name of the grism; empty if no spectrum could be read.
    """
    spectra = []
    for url in find_spectra(object_id)['url']:
        try:
            spectra.append(QTable.read(url))
        except OSError as error:
            print(f'Skipping a spectrum that could not be read ({error}): {url}')
    return best_spectra(spectra)
```

We plot the usable flux bins as solid lines and the flux bins not to use as dotted lines, set the y-axis range from the usable flux bins (ignoring the few highest ones, which are often spikes at the edges of the spectra), and mark the expected wavelengths of some of the strongest emission lines at the redshift of the galaxy:

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
grism_colors = {'Red grism': 'tab:red', 'Blue grism': 'tab:blue'}
emission_lines = {'[OII]': 3728.5 * u.AA, 'Hβ': 4862.7 * u.AA, '[OIII]': 5008.2 * u.AA,
                  'Hα': 6564.6 * u.AA, '[SII]': 6725.5 * u.AA}


def plot_spectra(ax, spectra, title, line_wavelengths):
    """Plot the spectra of an object on an axis, and mark its emission lines.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to plot on.
    spectra : dict
        Spectra to plot, keyed by the name of their grism, as returned by `best_spectra`.
    title : str
        Title of the plot.
    line_wavelengths : dict
        Observed wavelength of each line to mark, keyed by the label of the line.
        Only the lines in the wavelength range of the spectra are marked.
    """
    for name, spectrum in spectra.items():
        wavelength = spectrum['WAVELENGTH'].to(u.micron)
        bad = bad_flux_bins(spectrum)
        ax.plot(wavelength, np.ma.masked_where(bad, spectrum['SIGNAL']), color=grism_colors[name],
                alpha=0.8, linewidth=0.8)
        ax.plot(wavelength, np.ma.masked_where(~bad, spectrum['SIGNAL']), color=grism_colors[name],
                alpha=0.5, linewidth=0.8, linestyle=':')

    good_flux = np.concatenate([spectrum['SIGNAL'][~bad_flux_bins(spectrum)].filled(np.nan)
                                for spectrum in spectra.values()])
    low, high = np.nanpercentile(good_flux, [1, 99.5])
    ax.set_ylim(low - 0.1 * (high - low), high + 0.3 * (high - low))

    # Mark the lines in the wavelength range, with a single label for lines closer than 100 Angstrom
    wave_min = min(spectrum['WAVELENGTH'].min() for spectrum in spectra.values())
    wave_max = max(spectrum['WAVELENGTH'].max() for spectrum in spectra.values())
    labels = []
    for line, wavelength in sorted(line_wavelengths.items(), key=lambda item: item[1]):
        if wave_min < wavelength < wave_max:
            ax.axvline(wavelength.to(u.micron), color='grey', linestyle='--', linewidth=0.6)
            if labels and wavelength - labels[-1][1] < 100 * u.AA:
                labels[-1][0] += f'+{line}'
            else:
                labels.append([line, wavelength])
    for label, wavelength in labels:
        ax.text(wavelength.to_value(u.micron), 0.98, f' {label}', transform=ax.get_xaxis_transform(),
                rotation=90, va='top', ha='right', fontsize=9)
    ax.set_title(title, fontsize='medium')


def plot_galaxy(ax, galaxy, title=None):
    """Read and plot the best spectra of a galaxy, with the expected wavelengths of its emission lines.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to plot on.
    galaxy : astropy.table.Row
        Row of the galaxy, with its object ID, redshift, and Hα signal-to-noise ratio.
    title : str, optional
        First line of the title of the plot, above the description of the galaxy.
    """
    z = galaxy['spe_z']
    description = f"Object {galaxy['object_id']}, z = {z:.3f}, Hα S/N = {galaxy['halpha_snr_gf']:.0f}"
    spectra = get_best_spectra(galaxy['object_id'])
    if not spectra:
        ax.set_title(f'{description}\nNo spectrum could be read', fontsize='medium')
        return
    lines = {line: rest_wavelength * (1 + z) for line, rest_wavelength in emission_lines.items()}
    plot_spectra(ax, spectra, f'{title}\n{description}' if title else description, lines)


def add_legend(fig):
    """Add a legend of the grism colors and line styles to a figure, outside of its axes.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to add the legend to; it should use the constrained layout.
    """
    handles = [Line2D([], [], color=color, label=name) for name, color in grism_colors.items()]
    handles.append(Line2D([], [], color='grey', linestyle=':', label='Do not use'))
    handles.append(Line2D([], [], color='grey', linestyle='--', linewidth=0.6, label='Emission line'))
    fig.legend(handles=handles, loc='outside right upper', fontsize='small')


def fails_only(cuts, criterion):
    """Find the galaxies that fail one criterion and pass all of the others.

    Parameters
    ----------
    cuts : dict
        Boolean array of each criterion, keyed by the name of the criterion, as returned by `all_criteria`.
    criterion : str
        Name of the criterion to fail.

    Returns
    -------
    numpy.ndarray
        Boolean array, True for the galaxies that fail only this criterion.
    """
    others = np.logical_and.reduce([cut for name, cut in cuts.items() if name != criterion])
    return others & ~cuts[criterion]
```

### 5.1 Galaxies passing all of the criteria

Here are six of the galaxies that pass all of the criteria: five spread evenly from the strongest to the weakest Hα signal-to-noise ratio, and the galaxy with the strongest Hα line among those at z < 1.08, where Hα falls in both grisms:

```{code-cell} ipython3
examples = strong_halpha[np.linspace(0, len(strong_halpha) - 1, 5).astype(int)]
in_both_grisms = strong_halpha[(strong_halpha['spe_z'] < 1.08) & ~np.isin(strong_halpha['object_id'], examples['object_id'])]
examples = vstack([examples, in_both_grisms[:1]])

fig, axes = plt.subplots(3, 2, figsize=(14, 11), sharex=True, layout='constrained')
for ax, example in zip(axes.flat, examples):
    plot_galaxy(ax, example)
# Only label the outer axes, and hide the unused ones
for ax in axes[:, 1:].flat:
    ax.set_ylabel('')
for ax in axes[:-1, :].flat:
    ax.set_xlabel('')
for ax in axes.flat[len(examples):]:
    ax.set_visible(False)
add_legend(fig)
plt.show()
```

Hα is detected in all of them, down to a signal-to-noise ratio of about 10.
Where the Hα line falls in both grisms (1.21 to 1.37 micron), it is seen at the same wavelength in the red and the blue grism spectra, and in some of the galaxies other lines, such as [OIII] or [OII] in the blue grism, confirm the redshift.

```{note}
The blue and red grism spectra of a galaxy can have different continuum levels where they overlap (1.21 to 1.37 micron), sometimes by a factor of two or more, and the blue grism spectra can have spikes and other artifacts that are not flagged in the MASK column.
This is a property of the DR1 spectra, not a problem with how they are read in here: see Scodeggio et al. (2026) and the [Euclid DR1 documentation](https://euclid.caltech.edu/page/euclid-dr1) for the details.
```

### 5.2 Galaxies failing one of the criteria

Now let's look at galaxies that pass all of the criteria but one: for two of the criteria, we take the galaxy with the highest Hα signal-to-noise ratio among those that fail only this criterion.

```{code-cell} ipython3
fig, axes = plt.subplots(1, 2, figsize=(14, 4.5), sharex=True, layout='constrained')
for ax, criterion in zip(axes, ['3. Both grisms, full coverage and depth', '5. Consistent with the photo-z']):
    rejected = galaxies[fails_only(cuts, criterion)]
    if len(rejected) == 0:
        ax.set_title(f'No galaxy fails only criterion {criterion}', fontsize='medium')
        continue
    rejected.sort('halpha_snr_gf', reverse=True)
    plot_galaxy(ax, rejected[0], f'Fails criterion {criterion}')
axes[1].set_ylabel('')
add_legend(fig)
plt.show()
```

The galaxy on the left fails the grism, coverage, and depth criterion because its spectrum has 780 valid pixels and a median of 14 dithered exposures, just below our thresholds of 800 and 15: its blue grism spectrum is noisy, but Hα is detected in the red grism spectrum (S/N = 15), and its redshift is consistent with its photometric redshift. Thresholds like these also reject some usable spectra.
The galaxy on the right has a strong line at 1.23 micron, which the SPE processing identified as Hα, but its redshift is outside of the 95% interval of its photometric redshift.
The line is only seen in the blue grism spectrum, although the red grism spectrum also covers this wavelength: it is most likely not a line of this galaxy, but a contamination of the blue grism spectrum by a neighboring source, and the photometric redshift criterion rightly rejects this galaxy.

## 6. A closer look at the line measurements of one galaxy

Hα is only one of the lines the SPE processing measures.
Let's get all of the combined-grism line measurements of the galaxy with the strongest Hα line among those that pass all of the criteria, and keep those with a signal-to-noise ratio greater than 3:

```{code-cell} ipython3
strongest_galaxy = strong_halpha[0]

query_lines = f"""
SELECT spe_line_name, spe_line_central_wl_gf, spe_line_flux_gf, spe_line_flux_err_gf,
       spe_line_snr_gf, spe_line_ew_gf, spe_line_fwhm_gf, spe_line_qual_gf
FROM {table_lines}
WHERE object_id = {strongest_galaxy['object_id']}
  AND spe_rank = 0
  AND MOD(spe_line_qual_gf, 10) = 3
  AND spe_line_snr_gf > 3
ORDER BY spe_line_central_wl_gf
"""
galaxy_lines = Irsa.query_tap(query_lines).to_table()
galaxy_lines
```

Besides Hα, the [SII] doublet is well detected, with a signal-to-noise ratio of about 11.
[OIII] and [SIII]9068 also pass our cut, but with a signal-to-noise ratio between 4 and 5: these are marginal measurements, and the lines cannot be seen by eye in the spectrum below.
Note that the Hα and [NII]6585 lines have the same signal-to-noise ratio and FWHM: they are fitted together, as they are blended at the resolution of the spectra.
The table gives an [NII]/Hα flux ratio of about 0.6 for this galaxy, but this ratio comes from deblending the two lines and is poorly constrained (see the warning below).

We plot the spectra of this galaxy again, now marking the lines at the central wavelengths measured by the Gaussian fits, with the same labels as above for the main lines of the catalog:

```{code-cell} ipython3
line_labels = {'LyA': 'Lyα', 'OII3727': '[OII]', 'Hbeta': 'Hβ', 'OIII5008': '[OIII]', 'Halpha': 'Hα',
               'NII6585': '[NII]', 'SII6718': '[SII]', 'SIII9068': '[SIII]9068', 'SIII9530': '[SIII]9530',
               'PaGamma': 'Paγ', 'PaBeta': 'Paβ'}
measured_lines = {line_labels.get(name, name): wavelength for name, wavelength
                  in zip(galaxy_lines['spe_line_name'], u.Quantity(galaxy_lines['spe_line_central_wl_gf']))}

fig, ax = plt.subplots(figsize=(13, 4.5), layout='constrained')
plot_spectra(ax, get_best_spectra(strongest_galaxy['object_id']),
             f"Object {strongest_galaxy['object_id']}, z = {strongest_galaxy['spe_z']:.3f}: lines with S/N > 3", measured_lines)
add_legend(fig)
plt.show()
```

```{warning}
Be careful when using line ratios from the DR1 SPE lines catalog.
At the resolution of the NISP spectra, Hα is blended with the [NII] doublet: they are fitted together, with the same width and fixed relative positions, so the [NII]/Hα ratio comes from deblending the two and is poorly constrained.
For doublets fitted together, such as [OIII] and [NII], the flux is the total flux of the doublet: multiply it by 2/3 to get the flux of [OIII]5008 or [NII]6585 alone (Le Brun et al. 2026).
Which grism a line falls in depends on the redshift. In the Euclid Deep Fields, Hβ is in the blue grism spectrum at about 0.89 < z < 1.82 and [OIII] at about 0.84 < z < 1.74, while Hα is only in the red grism spectrum at z > 1.08. For the galaxies selected here (0.85 < z < 1.85) that are at z > 1.08, ratios such as Hα/Hβ therefore compare a line in the red grism with one in the blue grism, and also depend on the relative calibration of the two grisms.
```

## References

- Euclid Collaboration: Le Brun, V., Moresco, M., Bethermin, M., et al. 2026, *Euclid Data Release 1 (DR1): Performances of the redshift measurement and spectroscopic analysis*
<!-- TODO: add the link to the DR1 SPE paper once it is public. -->
<!-- TODO: update the Scodeggio et al. reference (journal, volume, arXiv ID) once the paper is public. -->
- Euclid Collaboration: Scodeggio, M., Copin, Y., et al. 2026, *Euclid Data Release 1 (DR1): From spectrograms to spectra: The SIR spectroscopic processing function*, A&A, in prep.
- Euclid Collaboration: Le Brun, V., Bethermin, M., Moresco, M., et al. 2026, A&A, 711, A7, *Euclid Quick Data Release (Q1): Characteristics and limitations of the spectroscopic measurements*, [arXiv:2503.15308](https://arxiv.org/abs/2503.15308)
- [Euclid DR1 SPE output catalog data product description](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmdr1/spedpd/dpcards/spe_spepfoutputcatalog.html), which describes the columns of the SPE catalogs and the line quality flags

+++

## About this Notebook

**Updated**: 2026-11-12

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.

**Runtime:** As of the date above, this notebook takes about 2 to 3 minutes, depending on the load of the IRSA services, to run to completion on a machine with 8GB RAM and 4 CPU.
