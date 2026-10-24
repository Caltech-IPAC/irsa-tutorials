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
**1. Remove the cell that points `Irsa.tap_url` and `Irsa.sia_url` to irsadev, once DR1 is public on the production IRSA services, then re-run the notebook against production and check the outputs and the plots.**

**2. Update the Kümmel et al. reference (journal, volume, arXiv ID) and verify the DR1 DOI.**
```

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- Understand the basic characteristics of Euclid DR1 MER catalogs.
- Know which of the four DR1 MER catalog tables to use.
- What columns are available in the MER catalog.
- How to search the MER catalog around a position on the sky.
- How to use the MER quality columns to clean a catalog.
- How to select stars and compute their total fluxes in the Euclid and ground-based bands.
- How to look at the MER VIS images of catalog sources.
- How to query the MER catalog with ADQL and make a color-color diagram of stars.

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
In this notebook we follow a typical path through the MER catalog: we search for the sources around a position on the sky, clean the sample with the quality columns, select the stars, look at the spectral energy distributions and VIS images of a few of them, and finally make a color-color diagram of many stars.
The MER processing for DR1 is described in detail in Euclid Collaboration: Kümmel et al. (2026), see the [References](#references) section below.
If you have questions about this notebook, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

<!-- TODO: if IRSA adds DR1 morphology and cutouts tables, mention them here (Q1 had euclid_q1_mer_morphology and euclid_q1_mer_cutouts). -->

+++

## The MER catalog tables

The DR1 MER catalog is split into four tables, according to how the data were observed and which MER tiling was used.
MER processes the sky in tiles: **wide** tiles (tile IDs starting with `102`, about 32 x 32 arcmin) and smaller **deep** tiles (tile IDs starting with `101`, about 17 x 17 arcmin), used for fields too deep to process as wide tiles.

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

The MER catalog from the earlier Quick Release 1 (Q1) is also available at IRSA, as the single table `euclid_q1_mer_catalogue`.
All of the Q1 areas have been reprocessed for DR1, so Q1 and DR1 fluxes should not be mixed (see the [notes on the DR1 MER photometry](#notes-on-the-dr1-mer-photometry) below).

+++

## Imports

```{important}
We rely on ``astroquery`` features that have been recently added, so please make sure you have version v0.4.10 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed
# !pip install numpy matplotlib astropy 'astroquery>=0.4.10' fsspec s3fs
```

```{code-cell} ipython3
import json

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa
```

<!-- TODO: remove this cell when DR1 is public on the production IRSA services. -->

```{code-cell} ipython3
# DR1 is not public yet, so we query IRSA's development server.
Irsa.tap_url = 'https://irsadev.ipac.caltech.edu/TAP'
Irsa.sia_url = 'https://irsadev.ipac.caltech.edu/SIA'
```

## 1. Find the MER catalog tables at IRSA

+++

The DR1 MER catalog is far too large to download, so we use IRSA's services to run queries on the server and retrieve only the rows and columns we need.

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

## 2. Search the MER catalog around a position

We start with all of the sources within 0.5 deg of a position in the Wide Survey.
This area contains a few hundred thousand sources, so rather than retrieving all 480 columns, we write an ADQL query that returns only the columns we need for this notebook:
the positions and tile IDs, the quality columns, the star-galaxy classification, and the template-fitting fluxes in the Euclid and ground-based bands (see section 5).
The `CONTAINS` and `CIRCLE` functions select the sources in a circle on the sky.

<!-- TODO: check that this position is in the production DR1 Wide Survey table. -->

```{code-cell} ipython3
coord = SkyCoord(ra=90.7299, dec=-45.3529, unit='deg')
search_radius = 0.5 * u.deg

# Template-fitting fluxes in the Euclid NIR bands and in all of the ground-based bands
flux_cols = ([f'flux_{band}_templfit' for band in 'yjh'] + [f'flux_vis_to_{band}_templfit' for band in 'yjh']
             + [col for col in columns_info
                if col.startswith('flux_') and '_ext_' in col and col.endswith('_templfit')])
columns = (['object_id', 'ra', 'dec', 'tileid', 'vis_det', 'flux_detection_total',
            'point_like_flag', 'point_like_prob', 'spurious_prob', 'det_quality_flag',
            'flag_vis', 'flag_y', 'flag_j', 'flag_h'] + flux_cols)

adql = (f"SELECT {', '.join(columns)} "
        f"FROM {table_mer} "
        "WHERE CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {coord.ra.deg}, {coord.dec.deg}, {search_radius.to_value(u.deg)})) = 1")

# Alternatively, uncomment the query below to get 10,000 stars from anywhere in the table instead of all of the
# sources around a position. The sky map below assumes a search around a position, so skip it in that case.
# adql = f"SELECT TOP 10000 {', '.join(columns)} FROM {table_mer} WHERE point_like_flag = 1"
```

We run the query with the TAP service, which takes about a minute.

```{code-cell} ipython3
result = Irsa.query_tap(adql).to_table()

# Keep only the sources with a positive total flux in the detection band, so that we can compute magnitudes
result = result[result['flux_detection_total'] > 0]
print(len(result))
```

```{tip}
For a small search area, the `query_region` method of `astroquery` is a simpler alternative to writing ADQL, and it can return all of the columns, for example:

    Irsa.query_region(coordinates=coord, spatial='Cone', catalog=table_mer, radius=3 * u.arcmin, columns='*')
```

Let's look at the density of sources on the sky, as offsets from the center of the search.
The MER catalog is built tile by tile, and the map shows the tile structure, areas without sources (for example around bright stars), and the edge of the DR1 coverage.

```{code-cell} ipython3
source_coords = SkyCoord(ra=result['ra'], dec=result['dec'], unit='deg')
offset_ra, offset_dec = coord.spherical_offsets_to(source_coords)
offset_ra = offset_ra.to_value(u.arcmin)
offset_dec = offset_dec.to_value(u.arcmin)

limit = 1.05 * search_radius.to_value(u.arcmin)
fig, ax = plt.subplots(figsize=(7, 6))
_, _, _, image = ax.hist2d(offset_ra, offset_dec, bins=np.linspace(-limit, limit, 120), cmap='Greys')
fig.colorbar(image, ax=ax, label='Number of sources per bin')
ax.add_patch(plt.Circle((0, 0), search_radius.to_value(u.arcmin), fill=False, linestyle='--', color='tab:blue'))
ax.set_xlim(limit, -limit)  # RA increases to the left
ax.set_ylim(-limit, limit)
ax.set_aspect('equal')
ax.set_xlabel('RA offset (arcmin)')
ax.set_ylabel('Dec offset (arcmin)')
```

## 3. Clean the sample with the quality columns

Before using the sources, we check the MER columns that identify sources whose detection or photometry may be unreliable:

- `spurious_prob` is the probability that a source is a spurious detection, for example an artifact around a bright object.
  `spurious_flag` is set for `spurious_prob > 0.5`, but this has proven too aggressive, and a threshold of 0.7 is suggested instead (see the [notes on the DR1 MER photometry](#notes-on-the-dr1-mer-photometry) below).
- `flag_vis`, `flag_y`, `flag_j`, and `flag_h` keep track of the flagged pixels (for example bad or saturated pixels) in the flag image of each band.
- `det_quality_flag` flags problems in the detection step that could corrupt the magnitude used for the star-galaxy separation.

Nonzero flags do not necessarily mean that a source is unusable, but they are worth checking for your science case.

<!-- TODO: describe the bits of det_quality_flag and of the flag_<band> columns once the DR1 data product description is available. -->

```{code-cell} ipython3
is_spurious = result['spurious_prob'] > 0.7
print(f"spurious_prob > 0.5: {np.sum(result['spurious_prob'] > 0.5)} of {len(result)} sources")
print(f"spurious_prob > 0.7: {np.sum(is_spurious)} of {len(result)} sources")
print(f"flag_vis != 0: {np.sum(result['flag_vis'] != 0)} of {len(result)} sources")
print(f"det_quality_flag != 0: {np.sum(result['det_quality_flag'] != 0)} of {len(result)} sources")
```

We remove the likely spurious sources from our sample:

```{code-cell} ipython3
clean_sources = result[~is_spurious]
print(len(clean_sources))
```

## 4. Select the stars

The MER catalog does not have a column for classification, but it provides `point_like_prob`, the probability that a VIS-detected source is point-like.
`point_like_flag = 1` selects sources with `point_like_prob > 0.96`, the threshold that reached a 99% star purity in simulations.
A lower `point_like_prob` threshold gives a more complete but less pure star sample.

```{code-cell} ipython3
stars = clean_sources[clean_sources['point_like_flag'].filled(0) == 1]
print(len(stars))
```

## 5. Look at the SEDs and VIS images of three stars

### Total fluxes in the Euclid and ground-based bands

Besides the Euclid VIS, Y, J, and H bands, the MER catalog contains photometry measured on ground-based images in the u, g, r, i, and z bands.
Which ground-based surveys are available depends on the position on the sky, for example DECam in the south, or a combination of MegaCam, HSC, and Pan-STARRS in the north.
The catalog therefore has columns for all of the possible survey and band combinations, named `flux_<band>_ext_<survey>_...`, and only some of them are filled at a given position.
Let's check which ground-based bands have measurements for our stars.

```{code-cell} ipython3
ext_flux_cols = [col for col in stars.colnames
                 if col.startswith('flux_') and '_ext_' in col and col.endswith('_templfit')
                 and not col.startswith('flux_vis_to_')]

for col in ext_flux_cols:
    n_measured = np.sum(~np.ma.getmaskarray(stars[col]))
    print(f"{col}: {n_measured} of {len(stars)} stars")
```

At this position, the ground-based photometry comes from DECam, in the g, r, i, and z bands.

To get total fluxes, we follow the recipe of the [MER Photometry Cookbook](https://euclid.esac.esa.int/dr/q1/dpdd/merdpd/merphotometrycookbook.html).
The total flux of an object is measured with a Kron aperture in the detection band, which is VIS for VIS-detected objects (`vis_det = 1`), and is given in `flux_detection_total`.
The total flux in another band X, from Euclid or from the ground, is then obtained by scaling the detection-band total flux with the color measured by template fitting:

$$F_{X,\mathrm{tot}} = \mathtt{flux\_detection\_total} \times \frac{\mathtt{flux\_X\_templfit}}{\mathtt{flux\_vis\_to\_X\_templfit}}$$

where `flux_vis_to_X_templfit` is the template-fitting flux measured on the VIS image degraded to the resolution of band X.
The cookbook gives an equivalent recipe with the aperture fluxes (`flux_X_nfwhm_aper / flux_vis_nfwhm_aper`).
The fluxes are in uJy, which we convert to AB magnitudes.

<!-- TODO: check whether there is a DR1 version of the MER Photometry Cookbook and update the link, and confirm with MER whether
they recommend the template-fitting or the aperture version of the recipe for DR1. -->

```{code-cell} ipython3
def total_mag(sources, band):
    """Total AB magnitude of VIS-detected sources in a band, following the MER Photometry Cookbook."""
    flux = sources['flux_detection_total']
    if band != 'vis':
        flux = flux * sources[f'flux_{band}_templfit'] / sources[f'flux_vis_to_{band}_templfit']
    return -2.5 * np.log10(flux) + 23.9
```

### Choose three stars

We pick three stars with positive fluxes in all bands, no flagged pixels in the Euclid bands, no detection problems, and of intermediate brightness.
To span the range of stellar colors, we take the stars closest to three points along the stellar locus in the VIS - J vs. J - H color-color diagram that we make below: a blue, an intermediate, and a red star.

```{code-cell} ipython3
# Approximate central wavelengths of the bands, in micron
band_wavelengths = {'g_ext_decam': 0.48, 'r_ext_decam': 0.64, 'vis': 0.72, 'i_ext_decam': 0.78,
                    'z_ext_decam': 0.92, 'y': 1.08, 'j': 1.37, 'h': 1.77}

has_all_bands = np.all([np.ma.filled(stars[f'flux_{band}_templfit'], 0) > 0
                        for band in band_wavelengths if band != 'vis'], axis=0)
mag_vis = total_mag(stars, 'vis')
no_flags = ((stars['flag_vis'] == 0) & (stars['flag_y'] == 0) & (stars['flag_j'] == 0) & (stars['flag_h'] == 0)
            & (stars['det_quality_flag'] == 0))
candidates = stars[has_all_bands & no_flags & (mag_vis > 18.5) & (mag_vis < 21)]

color_vis_j = total_mag(candidates, 'vis') - total_mag(candidates, 'j')
color_j_h = total_mag(candidates, 'j') - total_mag(candidates, 'h')
# (VIS - J, J - H) of a blue, an intermediate, and a red star on the stellar locus
locus_points = [(0.0, -0.2), (0.5, -0.03), (1.5, 0.03)]
three_stars = candidates[[np.argmin(np.hypot(color_vis_j - x, color_j_h - y)) for x, y in locus_points]]
three_stars['object_id', 'ra', 'dec', 'tileid']
```

The MER `object_id` is built from the position of the source: its digits are the RA and the absolute value of the Dec, both in units of 1e-7 deg, and its sign is the sign of the Dec.
All of our sources are in the southern sky, so they have negative IDs.

### Get VIS cutouts of the three stars

The MER mosaics are available through IRSA's Simple Image Access (SIA) service and are hosted in the cloud (see the MER mosaics notebook for more details).
For each star, we search for the wide VIS mosaics at its position, and read a small cutout directly from the cloud, without downloading the full mosaic.
Each source in the MER catalog belongs to one tile, so we use the mosaic of the tile given by its `tileid`.

```{code-cell} ipython3
def get_s3_fpath(cloud_access):
    cloud_info = json.loads(cloud_access)  # converts str to dict
    return f"{cloud_info['aws']['bucket_name']}/{cloud_info['aws']['key']}"


def vis_cutout(source, size=3 * u.arcsec):
    position = SkyCoord(ra=source['ra'], dec=source['dec'], unit='deg')
    sia_results = Irsa.query_sia(pos=(position, 1 * u.arcsec), collection='euclid_dr1_wide_DpdMerBksMosaic')
    # The tile ID is the first part of the obs_id of the mosaic
    is_vis_mosaic = ((sia_results['dataproduct_subtype'] == 'science') & (sia_results['energy_bandpassname'] == 'VIS')
                     & [obs_id.startswith(str(source['tileid'])) for obs_id in sia_results['obs_id']])
    mosaic = sia_results[is_vis_mosaic][0]
    with fits.open(f"s3://{get_s3_fpath(mosaic['cloud_access'])}", fsspec_kwargs={'anon': True}) as hdu:
        return Cutout2D(hdu[0].section, position=position, size=size, wcs=WCS(hdu[0].header)).data


cutouts = [vis_cutout(star) for star in three_stars]
```

### Plot the SEDs and the images

```{code-cell} ipython3
star_colors = ['tab:blue', 'tab:green', 'tab:red']

# The SEDs on the left, and the three VIS cutouts stacked in a narrow column on the right
fig = plt.figure(figsize=(12, 6))
grid = fig.add_gridspec(3, 2, width_ratios=[4, 1])
ax_sed = fig.add_subplot(grid[:, 0])
for star, color in zip(three_stars, star_colors):
    mags = [total_mag(star, band) for band in band_wavelengths]
    ax_sed.plot(list(band_wavelengths.values()), mags, 'o-', color=color, label=f"object_id {star['object_id']}")
ax_sed.invert_yaxis()
ax_sed.set_xscale('log')
ax_sed.set_xticks(list(band_wavelengths.values()), ['g', 'r', 'VIS', 'i', 'z', 'Y', 'J', 'H'])
ax_sed.minorticks_off()
ax_sed.set_xlabel('Band (approximate central wavelength on a log scale)')
ax_sed.set_ylabel('Total AB magnitude')
ax_sed.legend()

for i, (star, cutout, color) in enumerate(zip(three_stars, cutouts, star_colors)):
    ax = fig.add_subplot(grid[i, 1])
    norm = ImageNormalize(cutout, interval=PercentileInterval(99.9), stretch=AsinhStretch())
    ax.imshow(cutout, origin='lower', cmap='gray', norm=norm)
    ax.set_title(f"VIS = {total_mag(star, 'vis'):.1f}", color=color, fontsize='medium')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor(color)
        spine.set_linewidth(3)
fig.tight_layout()
```

All three stars are point sources in the VIS images, but their SEDs have very different shapes.
The bluest star is hot, with a nearly flat SED in AB magnitudes, while the cooler stars are much fainter in the blue bands than in the near-infrared.
A color-color diagram captures these differences for many stars at once.

## 6. Make a color-color diagram of the stars

Finally, we make a color-color diagram of all of the stars in our search area.
Some sources have no measurement in a band, for example where a tile is only partially covered (see the [notes on the DR1 MER photometry](#notes-on-the-dr1-mer-photometry) below).
These fluxes are masked in the table, and faint sources can also have negative fluxes, so we keep only the stars with positive fluxes in all of the bands we need.

```{code-cell} ipython3
cc_bands = ['j', 'h', 'g_ext_decam', 'r_ext_decam']
# Masked fluxes are filled with 0, so they are removed as well
has_all_fluxes = np.all([np.ma.filled(stars[f'flux_{band}_templfit'], 0) > 0 for band in cc_bands], axis=0)
many_stars = stars[has_all_fluxes]
print(f"{len(many_stars)} of {len(stars)} stars have VIS, J, H, g, and r fluxes")

color_x = total_mag(many_stars, 'vis') - total_mag(many_stars, 'j')
color_y = total_mag(many_stars, 'j') - total_mag(many_stars, 'h')
color_g_r = total_mag(many_stars, 'g_ext_decam') - total_mag(many_stars, 'r_ext_decam')
```

We plot the stars as points colored by their ground-based g - r color, add contours of the density of stars, and mark the three stars from above.

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(7, 6))
image = ax.scatter(color_x, color_y, c=color_g_r, s=5, alpha=0.6, cmap='coolwarm', vmin=0, vmax=1.6, edgecolor='none')
fig.colorbar(image, ax=ax, label='g - r (DECam)')

# Contours of the star counts in coarser bins, which smooths them
counts, xedges, yedges = np.histogram2d(color_x, color_y,
                                        bins=[np.arange(-0.75, 2.6, 0.12), np.arange(-0.6, 0.65, 0.05)])
ax.contour((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2, counts.T,
           levels=[5, 15, 40], colors='black', linewidths=0.8)

for star, color in zip(three_stars, star_colors):
    ax.scatter(total_mag(star, 'vis') - total_mag(star, 'j'), total_mag(star, 'j') - total_mag(star, 'h'),
               s=250, marker='*', color=color, edgecolor='black', linewidth=1.5, zorder=3)

ax.set_xlim(-0.75, 2.5)
ax.set_ylim(-0.6, 0.6)
ax.set_xlabel('VIS - J')
ax.set_ylabel('J - H')
ax.set_title(f'{len(many_stars)} stars in the DR1 MER Wide Survey catalog -- IRSA')
```

The stars follow the stellar locus, a narrow sequence in color-color space set mainly by their temperature.
From hot to cool stars, both colors first become redder together, then J - H stays nearly constant while VIS - J keeps increasing: for the coolest stars, the M dwarfs, molecular absorption shapes the near-infrared spectrum and J - H no longer tracks temperature.
The longer wavelength baseline of VIS - J keeps separating these stars, and our three stars lie at different positions along the locus.
The ground-based g - r color, which also traces temperature, changes smoothly along the locus, from the blue, hot stars to the red, cool ones.

<!-- TODO: re-check the plot limits after the final DR1 run on production. -->

## Notes on the DR1 MER photometry

Kümmel et al. (2026) describe the DR1 MER photometry and its known limitations in detail.
The most important points for using the catalog are:

* **NIR calibration offsets**: the Euclid Y, J, and H fluxes have offsets of a few to about 13 percent with respect to Gaia, depending on the photometric method. The offsets are similar in the three bands, so NIR colors are more reliable than absolute NIR magnitudes.
* **Partial coverage**: some tiles are only partially covered, at the survey borders or in individual bands, for example after solar proton events. Requiring a detection in all bands therefore removes sources in these areas.
* **Spurious sources**: `spurious_flag` is set for `spurious_prob > 0.5`, which has proven too aggressive; a threshold of 0.7 is suggested instead. The spurious probability is optimized for wide data.

```{note}
If you have used the Euclid Quick Release 1 (Q1) MER catalog: the DR1 PSFs are no longer normalized to unit flux.
This corrects the large offsets seen in Q1 for bands with extended PSFs (such as Pan-STARRS i), but means that DR1 fluxes differ systematically from the Q1 ones, so the two releases should not be mixed.
Colors measured with aperture (A-PHOT) or template-fitting (T-PHOT) photometry are largely unaffected.
```

<!-- TODO: the MER paper points to a filling_factor column for partial coverage, but it is not in the irsadev tables
as of 2026-10-23; mention it here if IRSA adds it. -->

+++

## References

<!-- TODO: update the Kümmel et al. reference (journal, volume, arXiv ID) and verify the DR1 DOI once the paper and DR1 are public. -->

- Euclid Collaboration: Kümmel, M., et al. 2026, *Euclid Data Release (DR1): The MERge processing*, A&A, in prep.
- Euclid Collaboration: Romelli, E., Kümmel, M., Dole, H., et al. 2026, A&A, 711, A4 (the Q1 MER processing)
- Euclid Data Release 1, 2026, [doi:10.57780/esa-d079c51](https://doi.org/10.57780/esa-d079c51)

## About this Notebook

**Updated**: 2026-10-23

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
