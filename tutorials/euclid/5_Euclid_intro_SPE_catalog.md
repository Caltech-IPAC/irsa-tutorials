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
from urllib.error import HTTPError

import matplotlib.pyplot as plt
import numpy as np

from astropy.table import QTable, Table, join, vstack
from astropy import units as u
from astropy.utils.data import conf
from astropy.visualization import quantity_support

from astroquery.ipac.irsa import Irsa

# Let matplotlib pick up the axis units from the quantities we plot
quantity_support()

# The Euclid spectrum files can take a while to read,
# which can exceed astropy's default timeout limit. Increase it.
conf.remote_timeout = 120

#suppress warnings about deprecated units and cache
import warnings
warnings.filterwarnings(
    "ignore",
    message="The unit 'Angstrom' has been deprecated",
    category=u.UnitsWarning,
)

warnings.filterwarnings(
    "ignore",
    message="The unit 'erg' has been deprecated",
    category=u.UnitsWarning,
)
warnings.filterwarnings(
    "ignore",
    message="XDG_CACHE_HOME is set",
)
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
* `euclid_dr1_spectro_z_spe_classification_*` -- the classification of each source as a galaxy, star, or QSO from its spectrum. The galaxy model only fits emission lines, so galaxies with a strong continuum and weak lines are often best fitted by the QSO model, which in DR1 also includes templates of quiescent, star-forming, and type 2 AGN galaxies: not all of the sources classified as QSO are quasars.

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
The criteria are deliberately strict, and are not a definition of a "good" spectrum: a galaxy at z > 1.85, for example, has its Hα line beyond the end of the red grism spectrum and fails criterion 2, but can have an excellent spectrum and redshift from its [OIII] or [OII] lines, and a galaxy at z < 0.85 can have its Hα line well measured in the blue grism spectrum of the Euclid Deep Fields.
For another science case, drop or change the criteria that do not apply, for example the redshift range and the line in criteria 2 and 6.

1. **High-confidence redshift**: the best redshift candidate (`spe_rank = 0`) has a high probability (`spe_z_prob > 0.99`, as in Sect. 4.1 of the paper, while Sect. 4.2 uses 0.94; this is the fraction of the redshift probability distribution in its peak) and is not flagged as unreliable (`spe_z_rel > 0.05`, a reliability from a machine-learning classifier, which in DR1 is trained on simulations, so only a low threshold is recommended).
2. **Hα in the red grism spectrum**: the redshift is between 0.85 and 1.85, the range in which Hα falls in the red grism spectrum (1.21 to 1.89 micron), away from its edges. At z < 1.08, Hα is also in the blue grism spectrum (up to 1.37 micron).
3. **Both grisms, full wavelength coverage and depth**: in the quality table, the source has both a blue and a red grism spectrum (`spe_grism = 'BGS+RGS'`; drop this requirement for the Euclid Wide Survey, which only has red grism spectra), its combined spectrum has at least 800 valid pixels (`spe_npix >= 800`, out of at most about 850 for the combined blue and red grism spectrum), so almost none of its wavelength range is masked, and the median number of dithered exposures per pixel is at least 15 (`spe_n_dith_med >= 15`, about 20 for most sources in this region), so it reaches close to the full depth of the survey.
4. **Extended, non-spurious source**: in the MER catalog, the source is not point-like (`point_like_prob < 0.05`, which removes most stars, but also compact galaxies and quasars), is unlikely to be spurious (`spurious_prob < 0.5`), and is not excessively large (`kron_radius < 50`): a large galaxy is most likely nearby, and a redshift in the Hα range is then most likely wrong. These are the criteria of the Euclid Galaxy Clustering Science Working Group.
5. **Consistent with the photometric redshift**: the spectroscopic redshift is within the 95% interval of the photometric redshift of the PHZ catalog (`phz_95_int_1` to `phz_95_int_2`), a check against a wrong line identification: Le Brun et al. (2026) find that with this criterion, 98% of the redshifts are correct at 0.84 < z < 1.6 in the EDF-N.
6. **Strong, well-measured Hα line**: the Gaussian fit of the Hα line is good in at least one of the two grisms, and the combined measurement of the two grisms has a signal-to-noise ratio greater than 10 and a flux greater than 1.5×10<sup>-16</sup> erg s<sup>-1</sup> cm<sup>-2</sup>, as in the visually inspected sample of the paper.

### 3.1 Get the redshifts, spectrum quality, and morphology of the sources

Criteria 1 to 5 only involve one row per source, from the redshift, quality, MER, and PHZ tables, so we get all of them with a single query, which joins the tables on IRSA's server.
We use a left outer join for the PHZ catalog, to keep the sources without a photometric redshift:

```{code-cell} ipython3
query_candidates = f"""
SELECT mer.object_id, mer.ra, mer.dec,
       mer.point_like_prob, mer.spurious_prob, mer.kron_radius,
       galaxy.spe_z, galaxy.spe_z_prob, galaxy.spe_z_rel,
       quality.spe_grism, quality.spe_npix, quality.spe_n_dith_med,
       phz.phz_median, phz.phz_95_int_1, phz.phz_95_int_2
FROM {table_mer} AS mer
JOIN {table_galaxy_candidates} AS galaxy ON galaxy.object_id = mer.object_id
JOIN {table_quality} AS quality ON quality.object_id = mer.object_id
LEFT OUTER JOIN {table_photo_z} AS phz ON phz.object_id = mer.object_id
WHERE mer.ra BETWEEN {ra_min} AND {ra_max}
  AND mer.dec BETWEEN {dec_min} AND {dec_max}
  AND galaxy.spe_rank = 0
"""
candidates = Irsa.query_tap(query_candidates, async_job=True).to_table()
print(f'{len(candidates)} of the {n_mer} MER sources have an SPE redshift')
```

We now apply the criteria one at a time and count how many galaxies are left after each of them.
Missing values (for example a missing photometric redshift) are masked in the table; we fill them with NaN, so that they fail every comparison:

```{code-cell} ipython3
def values(table, column):
    return np.ma.filled(np.ma.asarray(table[column], dtype=float), np.nan)


def criteria(table):
    z = values(table, 'spe_z')
    return {
        '1. High-confidence redshift': (values(table, 'spe_z_prob') > 0.99) & (values(table, 'spe_z_rel') > 0.05),
        '2. Hα in the red grism (0.85 < z < 1.85)': (z > 0.85) & (z < 1.85),
        '3. Both grisms, full coverage and depth': ((table['spe_grism'] == 'BGS+RGS')
                                                    & (values(table, 'spe_npix') >= 800)
                                                    & (values(table, 'spe_n_dith_med') >= 15)),
        '4. Extended, non-spurious source': ((values(table, 'point_like_prob') < 0.05)
                                     & (values(table, 'spurious_prob') < 0.5)
                                     & (values(table, 'kron_radius') < 50)),
        '5. Consistent with the photo-z': (z > values(table, 'phz_95_int_1')) & (z < values(table, 'phz_95_int_2')),
    }


def cut_flow(cuts):
    remaining = np.ones(len(next(iter(cuts.values()))), dtype=bool)
    steps, counts = ['SPE redshift'], [len(remaining)]
    for name, cut in cuts.items():
        remaining &= cut
        steps.append(name)
        counts.append(remaining.sum())
    return Table({'criterion': steps, 'n_galaxies': counts})


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
    cuts = criteria(table)
    cuts['6. Strong, well-measured Hα'] = (np.ma.filled(table['halpha_good_in_a_grism'], False)
                                           & (values(table, 'halpha_snr_gf') > 10)
                                           & (values(table, 'halpha_flux_gf') > 1.5e-16))
    return cuts


cuts = all_criteria(galaxies)
passes_all = np.logical_and.reduce(list(cuts.values()))
cut_flow(cuts)
```

Out of the tens of thousands of MER sources in the region, about 130 galaxies pass all the criteria of this strict recipe.
The Hα criterion removes most of the galaxies left after the first five criteria, mainly because of the flux cut: most of the Hα lines are fainter than 1.5×10<sup>-16</sup> erg s<sup>-1</sup> cm<sup>-2</sup>, close to the flux limit of the Euclid Wide Survey (2×10<sup>-16</sup> erg s<sup>-1</sup> cm<sup>-2</sup>), which the deep fields go below (see section 4).
Let's sort the galaxies that pass all of the criteria by the signal-to-noise ratio of their Hα line:

```{code-cell} ipython3
strong_halpha = galaxies[passes_all]
strong_halpha.sort('halpha_snr_gf', reverse=True)
strong_halpha['object_id', 'spe_z', 'phz_median', 'halpha_flux_gf', 'halpha_snr_gf'][:10]
```

## 4. Redshifts and Hα fluxes

Let's compare the redshifts of the galaxies that pass all of the criteria with those of all the galaxies with a high-confidence redshift in the same range, and look at their Hα fluxes as a function of redshift:

```{code-cell} ipython3
high_confidence = cuts['1. High-confidence redshift'] & cuts['2. Hα in the red grism (0.85 < z < 1.85)']
has_halpha = values(galaxies, 'halpha_flux_gf') > 0

fig, (ax_z, ax_flux) = plt.subplots(1, 2, figsize=(13, 4.5))

bins = np.arange(0.85, 1.86, 0.05)
ax_z.hist(galaxies['spe_z'][high_confidence], bins=bins, color='lightgrey', label='High-confidence redshift')
ax_z.hist(galaxies['spe_z'][passes_all], bins=bins, color='tab:red', label='All criteria')
ax_z.set_xlabel('Spectroscopic redshift')
ax_z.set_ylabel('Number of galaxies')
ax_z.legend()

ax_flux.scatter(galaxies['spe_z'][high_confidence & has_halpha], values(galaxies, 'halpha_flux_gf')[high_confidence & has_halpha],
                s=5, color='lightgrey', label='High-confidence redshift')
points = ax_flux.scatter(galaxies['spe_z'][passes_all], values(galaxies, 'halpha_flux_gf')[passes_all],
                         s=12, c=values(galaxies, 'halpha_snr_gf')[passes_all], cmap='viridis', vmin=10, vmax=50,
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
We find and read in the spectra as in the [SIR 1D spectra notebook](3_Euclid_intro_1D_spectra.md):
the table `euclid.objectid_spectrafile_association_dr1` lists the file and HDU of each spectrum of a source, and IRSA's spectrumDM service converts a given HDU into a VOTable, whose URL we build in the query.
A galaxy in the Euclid Deep Fields can have spectra from both the red and the blue grism, and several spectra of the same grism, so we keep the spectrum with the most usable flux bins for each grism:

```{code-cell} ipython3
def read_spectra(object_id):
    query = f"""
        SELECT 'https://irsacloud.ipac.caltech.edu/api/spectrumdm/convert/' || path
               || '?dataset_id=euclid_combspec_dr1&hdu=' || hdu AS url
        FROM euclid.objectid_spectrafile_association_dr1
        WHERE objectid = {object_id}
    """
    urls = Irsa.query_tap(query).to_table()['url']
    spectra = []
    for url in urls:
        # TODO: remove the try/except once IRSA's spectrumDM service converts all of the DR1 HDUs
        # (a few of them currently return an HTTP 422 error)
        try:
            spectra.append(QTable.read(url))
        except HTTPError as error:
            print(f'Could not read a spectrum of object {object_id}: {error}')
    return spectra


def bad_flux_bins(spectrum):
    # Flux bins with an odd MASK value or a MASK value of 64 or more should not be used
    mask = spectrum['MASK'].value
    return (mask % 2 == 1) | (mask >= 64)


def grism_name(spectrum):
    # Blue grism spectra start below 1.1 micron
    return 'Blue grism' if spectrum['WAVELENGTH'].min() < 1.1 * u.micron else 'Red grism'


def best_spectra(spectra):
    best = {}
    for spectrum in sorted(spectra, key=lambda spectrum: np.sum(~bad_flux_bins(spectrum)), reverse=True):
        best.setdefault(grism_name(spectrum), spectrum)
    return best
```

We plot the usable flux bins as solid lines and the flux bins not to use as dotted lines, set the y-axis range from the usable flux bins (ignoring the few highest ones, which are often spikes at the edges of the spectra), and mark the expected wavelengths of some of the strongest emission lines at the redshift of the galaxy:

```{code-cell} ipython3
grism_colors = {'Red grism': 'tab:red', 'Blue grism': 'tab:blue'}
emission_lines = {'[OII]': 3728.5 * u.AA, 'Hβ': 4862.7 * u.AA, '[OIII]': 5008.2 * u.AA,
                  'Hα': 6564.6 * u.AA, '[SII]': 6725.5 * u.AA}


def plot_spectra(ax, spectra, title, line_wavelengths):
    for name, spectrum in spectra.items():
        wavelength = spectrum['WAVELENGTH'].to(u.micron)
        bad = bad_flux_bins(spectrum)
        ax.plot(wavelength, np.ma.masked_where(bad, spectrum['SIGNAL']), color=grism_colors[name],
                alpha=0.8, linewidth=0.8, label=name)
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
                rotation=90, va='top', ha='right', fontsize=8)
    ax.set_title(title, fontsize='small')


def plot_galaxy(ax, galaxy, title=None, best=None):
    z = galaxy['spe_z']
    lines = {line: rest_wavelength * (1 + z) for line, rest_wavelength in emission_lines.items()}
    if best is None:
        best = best_spectra(read_spectra(galaxy['object_id']))
    description = f"Object {galaxy['object_id']}, z = {z:.3f}, Hα S/N = {galaxy['halpha_snr_gf']:.0f}"
    plot_spectra(ax, best, f'{title}\n{description}' if title else description, lines)
```

### 5.1 Galaxies passing all of the criteria

Here are 20 of the galaxies that pass all of the criteria, spread evenly from the strongest to the weakest Hα signal-to-noise ratio (reading the spectra takes a minute or two):

```{code-cell} ipython3
n_rows, n_columns = 4, 5
positions = np.linspace(0, len(strong_halpha) - 1, n_rows * n_columns).astype(int)

fig, axes = plt.subplots(n_rows, n_columns, figsize=(25, 14), sharex=True)
for ax, position in zip(axes.flat, positions):
    # If one of the spectra of a galaxy cannot be read, plot the next galaxy of the list instead
    for galaxy in strong_halpha[position:]:
        best = best_spectra(read_spectra(galaxy['object_id']))
        if len(best) == 2:
            break
    plot_galaxy(ax, galaxy, best=best)
# Only label the outer axes
for ax in axes[:, 1:].flat:
    ax.set_ylabel('')
for ax in axes[:-1, :].flat:
    ax.set_xlabel('')
axes[0, 0].legend(loc='lower right', fontsize='small')
plt.tight_layout()
plt.show()
```

Hα is clearly detected in all of them, down to a signal-to-noise ratio of about 10.
Where the Hα line falls in both grisms (1.2 to 1.37 micron), it is seen at the same wavelength in the red and the blue grism spectra, and in many of the galaxies other lines, such as [OIII] in the blue grism, confirm the redshift.

### 5.2 Galaxies failing one of the criteria

Now let's look at galaxies that pass all of the criteria but one: for two of the criteria, we take the galaxy with the highest Hα signal-to-noise ratio among those that fail only this criterion.

```{code-cell} ipython3
def fails_only(cuts, criterion):
    others = np.logical_and.reduce([cut for name, cut in cuts.items() if name != criterion])
    return others & ~cuts[criterion]


fig, axes = plt.subplots(1, 2, figsize=(14, 4.5), sharex=True)
for ax, criterion in zip(axes, ['3. Both grisms, full coverage and depth', '5. Consistent with the photo-z']):
    rejected = galaxies[fails_only(cuts, criterion)]
    rejected.sort('halpha_snr_gf', reverse=True)
    plot_galaxy(ax, rejected[0], f'Fails criterion {criterion}')
axes[0].legend(loc='lower right', fontsize='small')
plt.tight_layout()
plt.show()
```

The galaxy on the left fails the grism, coverage, and depth criterion because its spectrum has 780 valid pixels and a median of 14 dithered exposures, just below our thresholds of 800 and 15: its blue grism spectrum is noisy, but its Hα line is clearly seen in the red grism spectrum, and its redshift is consistent with its photometric redshift. This is a reminder that the thresholds of a strict recipe like this one also reject some usable spectra.
The galaxy on the right has a strong line at 1.23 micron, which the SPE processing identified as Hα, but its redshift is outside of the 95% interval of its photometric redshift.
The line is only seen in the blue grism spectrum, although the red grism spectrum also covers this wavelength: it is most likely not a line of this galaxy, but a contamination of the blue grism spectrum by a neighboring source, and the photometric redshift criterion rightly rejects this galaxy.

## 6. A closer look at the line measurements of one galaxy

Hα is only one of the lines the SPE processing measures.
Let's get all of the combined-grism line measurements of the galaxy with the strongest Hα line among those that pass all of the criteria, and keep those with a signal-to-noise ratio greater than 3:

```{code-cell} ipython3
galaxy = strong_halpha[0]

query_lines = f"""
SELECT spe_line_name, spe_line_central_wl_gf, spe_line_flux_gf, spe_line_flux_err_gf,
       spe_line_snr_gf, spe_line_ew_gf, spe_line_fwhm_gf, spe_line_qual_gf
FROM {table_lines}
WHERE object_id = {galaxy['object_id']}
  AND spe_rank = 0
  AND MOD(spe_line_qual_gf, 10) = 3
  AND spe_line_snr_gf > 3
ORDER BY spe_line_central_wl_gf
"""
galaxy_lines = Irsa.query_tap(query_lines).to_table()
galaxy_lines
```

Besides Hα, lines such as [OIII] (in the blue grism) and [SII] are detected.
Note that the Hα and [NII]6585 lines have the same signal-to-noise ratio and FWHM: they are fitted together, as they are blended at the resolution of the spectra (see the warning below).

We plot the spectra of this galaxy again, now marking the lines at the central wavelengths measured by the Gaussian fits:

```{code-cell} ipython3
best = best_spectra(read_spectra(galaxy['object_id']))
measured_lines = {row['spe_line_name']: u.Quantity(galaxy_lines['spe_line_central_wl_gf'])[i]
                  for i, row in enumerate(galaxy_lines)}

fig, ax = plt.subplots(figsize=(13, 4.5))
plot_spectra(ax, best, f"Object {galaxy['object_id']}, z = {galaxy['spe_z']:.3f}: lines with S/N > 3", measured_lines)
ax.legend(loc='lower right', fontsize='small')
plt.show()
```

```{warning}
Be careful when using line ratios from the DR1 SPE lines catalog.
At the resolution of the NISP spectra, Hα is blended with the [NII] doublet: they are fitted together, with the same width and fixed relative positions, so the [NII]/Hα ratio comes from deblending the two and is poorly constrained.
For doublets fitted together, such as [OIII] and [NII], the flux is the total flux of the doublet: multiply it by 2/3 to get the flux of [OIII]5008 or [NII]6585 alone (Le Brun et al. 2026).
Hβ and [OIII] fall in the blue grism for most of the galaxies in the Euclid Deep Fields, while Hα falls in the red grism, so ratios such as Hα/Hβ also depend on the relative calibration of the two grisms.
```

## References

- Euclid Collaboration: Le Brun, V., Moresco, M., Bethermin, M., et al. 2026, *Euclid Data Release 1 (DR1): Performances of the redshift measurement and spectroscopic analysis*
<!-- TODO: add the link to the DR1 SPE paper once it is public. -->
- Euclid Collaboration: Le Brun, V., Bethermin, M., Moresco, M., et al. 2026, A&A, 711, A7, *Euclid Quick Data Release (Q1): Characteristics and limitations of the spectroscopic measurements*, [arXiv:2503.15308](https://arxiv.org/abs/2503.15308)
- [Euclid DR1 SPE output catalog data product description](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmdr1/spedpd/dpcards/spe_spepfoutputcatalog.html), which describes the columns of the SPE catalogs and the line quality flags

+++

## About this Notebook

**Updated**: 2026-11-09

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.

**Runtime:** As of the date above, this notebook takes about 2 to 5 minutes, depending on the load of the IRSA services, to run to completion on a machine with 8GB RAM and 4 CPU.
