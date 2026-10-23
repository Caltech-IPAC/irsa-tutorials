---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
authors:
  - name: Tiffany Meshkat
  - name: Anahita Alavi
  - name: Anastasia Laity
  - name: Andreas Faisst
  - name: Brigitta Sipőcz
  - name: Dan Masters
  - name: Harry Teplitz
  - name: Jaladh Singhal
  - name: Shoubaneh Hemmati
  - name: Vandana Desai
---

# Euclid DR1: MER Catalogs

```{danger} TODO before the DR1 release (remove this box when done)
**1. Remove the cell that points `Irsa.tap_url` to irsadev, once DR1 is public on the production IRSA services, then re-run the notebook against production and check the outputs and the plot.**

**2. The column names used in the query were checked against the irsadev TAP schema (2026-10-23). Re-check once the production tables are final.**

**3. Confirm with IRSA which table holds the deep visits (single visits of the deep fields processed like wide survey data), and where the EMDS / HSC fields (HSC3, HSC4) ended up. Update the table descriptions below accordingly.**

**4. Ask IRSA whether DR1 morphology and cutouts tables (as `euclid_q1_mer_morphology` / `euclid_q1_mer_cutouts` in Q1) will be added, and mention them in the introduction if so.**

**5. Update the Kümmel et al. reference (journal, volume, arXiv ID) and verify the DR1 DOI.**
```

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- Understand the basic characteristics of Euclid DR1 MER catalogs.
- Know which of the four DR1 MER catalog tables to use.
- What columns are available in the MER catalog.
- How to query with ADQL in the MER catalog.
- How to make a simple color-magnitude diagram with the data.

+++

## Introduction

+++

Euclid launched in July 2023 as a European Space Agency (ESA) mission with involvement by NASA.
The primary science goals of Euclid are to better understand the composition and evolution of the dark Universe.
The Euclid mission is providing space-based imaging and spectroscopy as well as supporting ground-based imaging to achieve these primary goals.
These data will be archived by multiple global repositories, including IRSA, where they will support transformational work in many areas of astrophysics.

Euclid Data Release 1 (DR1) comprises about 1900 sq deg of imaging, spectroscopy, and catalogs from the Euclid Wide Survey (EWS), which will eventually cover about 14,000 sq deg, and from the Euclid Deep Survey (EDS), which covers about 63 sq deg in the three Euclid Deep Fields (North, South, and Fornax).
DR1 also includes several photometric redshift and instrument calibration fields, such as COSMOS, CDFS, and AEGIS.
It follows Quick Release 1 (Q1), which covered four much smaller fields; all of the DR1 data, including the areas observed in Q1, have been reprocessed with updated pipelines and calibrations.

Among the data products included in DR1 are the MER catalogs.
The MER processing detects sources on the VIS and NISP mosaics (see the MER mosaics notebook) and measures their photometry in the Euclid VIS, Y, J, H bands and in the accompanying ground-based u, g, r, i, z bands, along with other basic measurements, like size and shape.
This notebook provides an introduction to the MER final catalog, in which each entry is a single source.
The MER processing for DR1 is described in detail in Euclid Collaboration: Kümmel et al. (2026), see the [References](#references) section below.
If you have questions about this notebook, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

<!-- TODO: if IRSA adds DR1 morphology and cutouts tables, mention them here (Q1 had euclid_q1_mer_morphology and euclid_q1_mer_cutouts). -->

+++

## The four DR1 MER catalog tables

Unlike Q1, which had a single MER final catalog table, DR1 splits the MER catalog into four tables, according to how the data were observed and which MER tiling was used.
MER processes the sky in tiles: **wide** tiles (tile IDs starting with `102`, about 32 x 32 arcmin, the same tiling as in Q1) and smaller **deep** tiles (tile IDs starting with `101`, about 17 x 17 arcmin), used for fields too deep to process as wide tiles.

* `euclid_dr1_mer_catalogue_wide_survey` -- **the Euclid Wide Survey**, the standard-depth survey intended for cosmology (wide tiles).
  **If you are looking for the standard Euclid survey, this is the table you want**, and it is the one used in this notebook.
* `euclid_dr1_mer_catalogue_deep_survey` -- **the Euclid Deep Survey**: the deep stacks of the Euclid Deep Fields North, South, and Fornax (deep tiles).
* `euclid_dr1_mer_catalogue_deep_mode` -- other fields deep enough to need the deep tiling, but not part of the Deep Survey: the Euclid auxiliary fields COSMOS, CDFS, and AEGIS.
* `euclid_dr1_mer_catalogue_wide_mode` -- other fields processed with the wide tiling that are neither part of the Wide Survey nor visits of the deep fields, for example some of the fields observed during the thermal stabilization period of the PSF calibration campaign, just before the main survey started (January/February 2024).

<!-- TODO: confirm with IRSA where the deep visits are (each visit of a deep field processed as though it were wide survey data,
identified by its sky patch; the ESA archive help puts them in wide mode, but Euclid MER says they should not be there),
and where the EMDS / HSC fields (HSC3, HSC4: wide tiling but with EXT bands that need the MER deep catalog) ended up.
Do not quote sky patch numbers until confirmed: the MER paper's "patch" numbers (44, 64, 71, 126) differ from the SkyPatch
values given for Q1 (49 EDF-N, 58 EDF-F, 62 EDF-S). Once confirmed, add a bullet for the deep visits here.
irsadev (2026-10-23) also has euclid_dr1_mer_catalogue_deep_survey_photometry and euclid_dr1_mer_catalogue_deep_mode_photometry
(716 columns each); ask IRSA what they hold (EXT bands for the deep tiles?) and describe them here. -->


The deep fields are covered by both the Deep Survey (deep tiles) and the wide tiling, so the same object can appear in more than one table, measured on data of different depth.
Table 5 of Kümmel et al. (2026) lists all of the DR1 fields with their tiling and exposure times.

+++

## Things to know about the DR1 MER photometry

Kümmel et al. (2026) describe several changes and known issues compared to Q1.
The most important ones for using the catalog are:

* **The PSFs are no longer normalized to unit flux**, which corrects the large offsets seen in Q1 for bands with extended PSFs (such as Pan-STARRS i), but means that DR1 fluxes differ systematically from the Q1 ones. Colors measured with aperture (A-PHOT) or template-fitting (T-PHOT) photometry are largely unaffected.
* **The Euclid NIR bands have offsets of a few to about 13 percent** with respect to Gaia, depending on the photometric method. The offsets are similar in the three NIR bands, so NIR colors are less affected than absolute magnitudes.
* **Some tiles are only partially covered**, at the survey borders or in individual bands, for example after solar proton events. Requiring a detection in all bands therefore removes more sources than in Q1.
* **Spurious sources**: `spurious_flag` is set for `spurious_prob > 0.5`, which has proven too aggressive; a threshold of 0.7 is suggested instead. The spurious probability is optimized for wide data.

<!-- TODO: the MER paper points to a filling_factor column for partial coverage, but it is not in the irsadev tables
as of 2026-10-23; mention it here if IRSA adds it. -->

+++

## Imports

```{important}
We rely on ``astroquery`` features that have been recently added, so please make sure you have version v0.4.10 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed
# !pip install numpy matplotlib 'astroquery>=0.4.10'
```

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

from astroquery.ipac.irsa import Irsa
```

<!-- TODO: remove this cell when DR1 is public on the production IRSA services. -->

```{code-cell} ipython3
# DR1 is not public yet, so we query IRSA's development server.
Irsa.tap_url = 'https://irsadev.ipac.caltech.edu/TAP'
```

## 1. Download MER catalog from IRSA directly to this notebook

+++

First, have a look at what Euclid catalogs are available. With the ``list_catalogs`` functionality, we'll receive a list of the name of the catalogs as well as their brief desciption.

```{code-cell} ipython3
Irsa.list_catalogs(filter='euclid')
```

### Choose the Euclid MER table

We use the Wide Survey table. To use another one of the four DR1 MER tables described above, change the table name here.

```{code-cell} ipython3
table_mer = 'euclid_dr1_mer_catalogue_wide_survey'
```

### Learn some information about the MER catalog:
- How many columns are there?
- List the column names

```{code-cell} ipython3
columns_info = Irsa.list_columns(catalog=table_mer)
print(len(columns_info))
```

```{code-cell} ipython3
columns_info
```

<!-- TODO: flux_radius, sersic_visnir_flags, sersic_ext_flags (new DR1 columns in the MER paper) are not in the
irsadev tables as of 2026-10-23; ask IRSA whether they will be added, and highlight them here if so. -->

```{tip}
Each DR1 MER catalog table contains 480 columns, below are a few highlights:

- object_id
- flux_detection_total: the total (Kron) flux in the detection band
- flux_y_templfit, flux_j_templfit, flux_h_templfit: template-fitting fluxes
- flux_vis_to_y_templfit, flux_vis_to_j_templfit, flux_vis_to_h_templfit: template-fitting fluxes measured on the VIS image degraded to the resolution of each NIR band
- fwhm
- point_like_flag, point_like_prob, spurious_prob
```

+++

### Define the following ADQL query to find the first 10k stars in the MER catalog

Since we are just using the MER catalog alone, it does not have a column for classification.
We can use `point_like_flag = 1` for stars, which selects VIS-detected objects with `point_like_prob > 0.96`, the threshold that reached a 99% star purity in simulations.
A lower `point_like_prob` threshold gives a more complete but less pure star sample.

To get total fluxes in the NIR bands, we follow the recipe of the [MER Photometry Cookbook](https://euclid.esac.esa.int/dr/q1/dpdd/merdpd/merphotometrycookbook.html).
The total flux of an object is measured with a Kron aperture in the detection band, which is VIS for VIS-detected objects (`vis_det = 1`), and is given in `flux_detection_total`.
The total flux in another band X is then obtained by scaling the detection-band total flux with the color measured by template fitting:

$$F_{X,\mathrm{tot}} = \mathtt{flux\_detection\_total} \times \frac{\mathtt{flux\_X\_templfit}}{\mathtt{flux\_vis\_to\_X\_templfit}}$$

where `flux_vis_to_X_templfit` is the template-fitting flux measured on the VIS image degraded to the resolution of band X.
The cookbook gives an equivalent recipe with the aperture fluxes (`flux_X_nfwhm_aper / flux_vis_nfwhm_aper`).

We select VIS-detected objects and require all of these fluxes to be greater than 0, so the object is detected in all four Euclid MER mosaic images.

<!-- TODO: check whether there is a DR1 version of the MER Photometry Cookbook and update the link, and confirm with MER whether
they recommend the template-fitting or the aperture version of the recipe for DR1. -->

```{code-cell} ipython3
adql_stars = ("SELECT TOP 10000 mer.object_id, mer.ra, mer.dec, "
              "mer.flux_detection_total, mer.fluxerr_detection_total, "
              "mer.flux_y_templfit,mer.fluxerr_y_templfit, mer.flux_j_templfit, mer.fluxerr_j_templfit, "
              "mer.flux_h_templfit, mer.fluxerr_h_templfit, "
              "mer.flux_vis_to_y_templfit, mer.flux_vis_to_j_templfit, mer.flux_vis_to_h_templfit, "
              "mer.point_like_prob, mer.extended_prob "
              f"FROM {table_mer} AS mer "
              "WHERE mer.vis_det = 1 "
              "AND mer.flux_detection_total > 0 "
              "AND mer.flux_y_templfit > 0 "
              "AND mer.flux_j_templfit > 0 "
              "AND mer.flux_h_templfit > 0 "
              "AND mer.flux_vis_to_y_templfit > 0 "
              "AND mer.flux_vis_to_j_templfit > 0 "
              "AND mer.flux_vis_to_h_templfit > 0 "
              "AND mer.point_like_flag = 1 ")
```

We can run the query with the TAP service, and then look at some of the results.

```{code-cell} ipython3
result_stars = Irsa.query_tap(adql_stars).to_table()
result_stars[:5]
```

## 2. Make a color-magnitude diagram using the catalogs pulled from IRSA

- Compute the total Y and H fluxes with the recipe above
- Convert from flux in uJy to magnitudes using the zero point correction
- Convert the error bars to magnitudes as well
- Plot the color-magnitude diagram

```{code-cell} ipython3
def total_flux(band):
    return (result_stars["flux_detection_total"] * result_stars[f"flux_{band}_templfit"]
            / result_stars[f"flux_vis_to_{band}_templfit"])


mag_y = -2.5 * np.log10(total_flux("y")) + 23.9
mag_h = -2.5 * np.log10(total_flux("h")) + 23.9

x = mag_y - mag_h  # Y - H
y = mag_y

# For simplicity, the error bars only include the template-fitting flux errors
xerr = (2.5 / np.log(10) * np.sqrt((result_stars["fluxerr_y_templfit"] / result_stars["flux_y_templfit"])**2
                                   + (result_stars["fluxerr_h_templfit"] / result_stars["flux_h_templfit"])**2))
yerr = (2.5 / np.log(10) * (result_stars["fluxerr_y_templfit"] / result_stars["flux_y_templfit"]))

plt.errorbar(x, y, xerr=xerr, yerr=yerr,
             fmt='o', markersize=1.5, ecolor='lightgrey', elinewidth=0.5, capsize=2)

plt.xlabel('Y-H')
plt.ylabel('Y')
# Note that these limits exclude a handful of points with large error bars.
plt.xlim(-2, 2)
plt.ylim(24, 16)
plt.title('10k Stars in DR1 MER Wide Survey catalog -- IRSA')
```

<!-- TODO: after running against DR1, check whether the plot limits (Y-H from -2 to 2, Y from 16 to 24) still suit the data. -->

## References

<!-- TODO: update the Kümmel et al. reference (journal, volume, arXiv ID) and verify the DR1 DOI once the paper and DR1 are public. -->

- Euclid Collaboration: Kümmel, M., et al. 2026, *Euclid Data Release (DR1): The MERge processing*, A&A, in prep.
- Euclid Collaboration: Romelli, E., Kümmel, M., Dole, H., et al. 2026, A&A, 711, A4 (the Q1 MER processing)
- Euclid Data Release 1, 2026, [doi:10.57780/esa-d079c51](https://doi.org/10.57780/esa-d079c51)

## About this Notebook

**Updated**: 2026-10-15

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
