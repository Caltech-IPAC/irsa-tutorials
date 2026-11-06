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
  - name: Shoubaneh Hemmati
  - name: Vandana Desai
  - name: Troy Raen
  - name: Jessica Krick
---

# Euclid DR1: PHZ Catalogs

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- know which Euclid DR1 PHZ (photometric redshift) catalog tables are available at IRSA, and how to view their columns;
- know how to join the PHZ and MER catalogs with ADQL to select galaxies with good photometric redshifts in a redshift range;
- know how to compare the photometric redshifts with the spectroscopic redshifts of the DR1 SPE catalog;
- know how to make a cutout of a MER mosaic and overplot the selected galaxies on it;
- know how to retrieve and plot the spectrum of one of the galaxies, and compare its photometric and spectroscopic redshifts;
- know how to load the image and catalog into Firefly to inspect individual sources interactively.

+++

## Introduction

+++

Euclid launched in July 2023 as a European Space Agency (ESA) mission with involvement by NASA.
The primary science goals of Euclid are to better understand the composition and evolution of the dark Universe.
The Euclid mission is providing space-based imaging and spectroscopy as well as supporting ground-based imaging to achieve these primary goals.
These data are archived by multiple global repositories, including IRSA, where they will support transformational work in many areas of astrophysics.

Euclid Data Release 1 (DR1) comprises about 1900 sq deg of imaging, spectroscopy, and catalogs from the Euclid Wide Survey (EWS), which will eventually cover about 14,000 sq deg, and from the Euclid Deep Survey (EDS), which covers about 53 sq deg in the three Euclid Deep Fields (North, South, and Fornax).
DR1 also includes several photometric redshift and instrument calibration fields, such as COSMOS, CDFS, and AEGIS.
It follows Quick Release 1 (Q1), which covered about 63 sq deg with a single pass over each of the three Euclid Deep Fields and the dark cloud LDN 1641; all of the DR1 data, including the areas observed in Q1, have been reprocessed with updated pipelines and calibrations.

Among the data products included in DR1 are the catalogs created by the PHZ processing function, which estimates the photometric redshift and the classification (star, galaxy, QSO) of every source in the MER catalog from its photometry.
This notebook provides an introduction to the main PHZ catalog, which contains the photometric redshift probability distribution, the unified fluxes, and the classification for each source.
We follow a typical path through the catalog: we select galaxies with good photometric redshifts in a redshift range, look at them on a MER mosaic, and compare the photometric redshift of one of them with its spectrum and spectroscopic redshift.
The PHZ processing for DR1 is described in detail in Euclid Collaboration: Tucci et al. (2026), see the [References](#references) section below.
If you have questions about this notebook, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## The PHZ catalog tables

Like the MER catalog (see the [MER catalogs notebook](2_Euclid_intro_MER_catalog.md)), the DR1 PHZ catalogs are split by survey:

* `euclid_dr1_phz_photo_z_wide` -- the main PHZ catalog of the **Euclid Wide Survey**: photometric redshifts, unified fluxes, and classification.
* `euclid_dr1_phz_photo_z_deep_survey` -- the same for the **Euclid Deep Survey** (the deep tiles of the Euclid Deep Fields North, South, and Fornax).
* `euclid_dr1_phz_classification_wide` and `euclid_dr1_phz_classification_deep_survey` -- the detailed classification probabilities (star, galaxy, QSO, globular cluster) behind the `phz_classification` column of the main catalogs.

The PHZ tables have no coordinates; they are joined to the MER catalog of the same survey on `object_id`.
In this notebook we work in the Euclid Deep Field North, so we use the deep survey tables.

```{note}
Some column names differ from the Q1 PHZ catalog: for example the redshift intervals are `phz_90_int_1` and `phz_90_int_2` in DR1 (`phz_90_int1` and `phz_90_int2` in Q1).
```

+++

## Imports

```{important}
We rely on ``astroquery`` features that have been recently added, so please make sure you have version v0.4.10 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy matplotlib 'astropy>=5.3' 'astroquery>=0.4.10' fsspec s3fs firefly_client
```

```{code-cell} ipython3
import io
import json
import warnings

import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import QTable, join
from astropy import units as u
from astropy.utils.data import conf
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch, LogStretch, quantity_support
from astropy.wcs import WCS

from firefly_client import FireflyClient
from astroquery.ipac.irsa import Irsa

# Suppress warnings about deprecated units in the spectrum files
warnings.filterwarnings('ignore', message="The unit 'Angstrom' has been deprecated", category=u.UnitsWarning)
warnings.filterwarnings('ignore', message="The unit 'erg' has been deprecated", category=u.UnitsWarning)

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
Irsa.sia_url = 'https://irsadev.ipac.caltech.edu/SIA'
Irsa.ssa_url = 'https://irsadev.ipac.caltech.edu/SSA'
```

<!-- TODO: remove this cell once the SSA service serves the DR1 spectra again (IRSA temporarily removed the euclid_dr1_DpdSirCombinedSpectra collection for fixes on 2026-11-06). -->

```{code-cell} ipython3
# TEMPORARY: the SSA service does not return DR1 spectra yet, so we emulate
# Irsa.query_ssa with a TAP query on IRSA's CAOM tables, which already have them.
from astropy.table import Table


def query_ssa_from_caom(*, pos, radius, collection):
    query = f"""
        SELECT a.uri
        FROM euclid.plane_euclid_dr1 AS p
        JOIN euclid.observation_euclid_dr1 AS o ON o.obsid = p.obsid
        JOIN euclid.artifact_euclid_dr1 AS a ON a.planeid = p.planeid
        WHERE o.collection = '{collection}'
          AND CONTAINS(p.pt, CIRCLE('ICRS', {pos.ra.deg}, {pos.dec.deg}, {radius.to(u.deg).value})) = 1
    """
    uris = Irsa.query_tap(query, async_job=True).to_table()['uri']
    return Table({'s_ra': [pos.ra.deg] * len(uris), 's_dec': [pos.dec.deg] * len(uris),
                  'access_url': [f'https://irsadev.ipac.caltech.edu/{uri}' for uri in uris]})


Irsa.query_ssa = query_ssa_from_caom
```

## 1. Find the PHZ catalog tables at IRSA

The Euclid catalogs at IRSA can be listed with the `list_catalogs` method; here we filter on the DR1 tables.

```{note}
In the example below, we use ``include_metadata_tables=True``, as IRSA lists the PHZ and SPE tables among the metadata tables rather than the catalogs.
This parameter was introduced in astroquery 0.4.11. For earlier versions, please omit this parameter to receive the combined list of available catalogs and metadata tables.
```

```{code-cell} ipython3
Irsa.list_catalogs(filter='euclid_dr1', include_metadata_tables=True)
```

```{code-cell} ipython3
table_mer = 'euclid_dr1_mer_catalogue_deep_survey'
table_phz = 'euclid_dr1_phz_photo_z_deep_survey'
```

### Learn some information about the photo-z catalog

Let's look at the number of columns and their descriptions:

```{code-cell} ipython3
columns_info = Irsa.list_columns(catalog=table_phz)
print(f'{len(columns_info)} columns')
```

```{tip}
A few highlights of the PHZ catalog columns:

- `object_id`, the MER object ID used to join the PHZ catalog with the MER catalog;
- `phz_median`, the median of the photometric redshift probability distribution (PDF);
- `phz_90_int_1` and `phz_90_int_2`, the lower and upper values of the PDF interval containing 90% of the probability;
- `phz_mode_1` and `phz_mode_2`, the two highest modes of the PDF;
- `phz_classification`: 1 for a star, 2 for a galaxy, 4 for a QSO, 8 for a globular cluster (a sum, for example 3, when the source is accepted as several classes);
- `phz_flags`, the PHZ processing flags, 0 when all went well;
- `flux_vis_unif`, `flux_y_unif`, `flux_j_unif`, `flux_h_unif` and the ground-based bands: the "unified" fluxes, recomputed from the MER fluxes after correction for Galactic extinction and filter shifts, which are the fluxes used for the photometric redshifts.
```

```{code-cell} ipython3
# Full list of columns and their description
columns_info
```

## 2. Select galaxies between z = 1.4 and 1.6 with good photometric redshifts

We join the PHZ catalog with the MER catalog, which has the coordinates, and select galaxies in a 5 arcmin box around a position in the Euclid Deep Field North, with the following conditions:
- the unified fluxes in the four Euclid bands are greater than zero, so the galaxies are detected in all four Euclid MER mosaics;
- `phz_classification = 2`, so the source is accepted as a galaxy only;
- `phz_flags = 0`, so the PHZ processing had no problems;
- the 90% interval of the photometric redshift PDF is narrower than 20% of (1 + z), so the photometric redshift is well constrained;
- the median photometric redshift is between 1.4 and 1.6.

```{code-cell} ipython3
######################## User defined section ############################
# Set the center and size of the region to select galaxies from (and of the image cutout below)
ra_cutout = 267.8
dec_cutout = 66
im_cutout = 5 * u.arcmin

# Redshift range
z_min, z_max = 1.4, 1.6
##########################################################################

coords_cutout = SkyCoord(ra_cutout, dec_cutout, unit='deg', frame='icrs')
size_cutout = im_cutout.to(u.deg).value
```

```{code-cell} ipython3
adql = f"""
    SELECT mer.object_id, mer.ra, mer.dec,
           phz.flux_vis_unif, phz.flux_y_unif, phz.flux_j_unif, phz.flux_h_unif,
           phz.phz_classification, phz.phz_median, phz.phz_90_int_1, phz.phz_90_int_2
    FROM {table_mer} AS mer
    JOIN {table_phz} AS phz ON mer.object_id = phz.object_id
    WHERE CONTAINS(POINT('ICRS', mer.ra, mer.dec),
                   BOX('ICRS', {ra_cutout}, {dec_cutout}, {size_cutout}, {size_cutout})) = 1
      AND phz.flux_vis_unif > 0
      AND phz.flux_y_unif > 0
      AND phz.flux_j_unif > 0
      AND phz.flux_h_unif > 0
      AND phz.phz_classification = 2
      AND phz.phz_flags = 0
      AND (phz.phz_90_int_2 - phz.phz_90_int_1) / (1 + phz.phz_median) < 0.20
      AND phz.phz_median BETWEEN {z_min} AND {z_max}
"""

result_galaxies = Irsa.query_tap(adql, async_job=True).to_table()
print(f'Found {len(result_galaxies)} galaxies')
result_galaxies[:5]
```

## 3. Compare the photometric and spectroscopic redshifts

The Euclid Deep Fields are also covered by NISP slitless spectroscopy, so many of these galaxies have a spectroscopic redshift in the DR1 SPE catalog (see the [SPE catalogs notebook](5_Euclid_intro_SPE_catalog.md)).
We query the best redshift estimate of each galaxy (`spe_rank = 0`) with its probability `spe_z_prob` and reliability `spe_z_rel`, and join the result to our table on `object_id`, keeping the galaxies without a spectroscopic redshift (`join_type='left'`):

```{code-cell} ipython3
adql_spe = f"""
    SELECT object_id, spe_z, spe_z_prob, spe_z_rel
    FROM euclid_dr1_spectro_z_spe_galaxy_candidates_deep
    WHERE spe_rank = 0
      AND object_id IN ({', '.join(str(obj_id) for obj_id in result_galaxies['object_id'])})
"""
result_spe = Irsa.query_tap(adql_spe, async_job=True).to_table()

galaxies = join(result_galaxies, result_spe, keys='object_id', join_type='left')
has_spec_z = ~galaxies['spe_z'].mask

print(f'{has_spec_z.sum()} of the {len(galaxies)} galaxies have a spectroscopic redshift')
```

Not all of these spectroscopic redshifts can be trusted: most of the galaxies are faint, and when a single emission line is detected, it can be attributed to the wrong line, which gives a completely wrong redshift.
As in the [SPE catalogs notebook](5_Euclid_intro_SPE_catalog.md), we consider a spectroscopic redshift reliable when both its probability and its reliability are high (`spe_z_prob > 0.99` and `spe_z_rel > 0.05`):

```{code-cell} ipython3
reliable = has_spec_z & (galaxies['spe_z_prob'].filled(0) > 0.99) & (galaxies['spe_z_rel'].filled(0) > 0.05)

print(f'{reliable.sum()} of the {has_spec_z.sum()} spectroscopic redshifts are reliable')
```

Let's compare the two redshifts:

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(galaxies['spe_z'][has_spec_z], galaxies['phz_median'][has_spec_z],
           facecolors='none', edgecolors='grey', label='All spectroscopic redshifts')
ax.scatter(galaxies['spe_z'][reliable], galaxies['phz_median'][reliable],
           color='tab:blue', label='Reliable spectroscopic redshifts')
ax.plot([0, 5], [0, 5], color='grey', linestyle='--', linewidth=0.8)
ax.axhspan(z_min, z_max, color='tab:red', alpha=0.1)
ax.set_xlabel('Spectroscopic redshift (spe_z)')
ax.set_ylabel('Photometric redshift (phz_median)')
_ = ax.legend(loc='upper left')
```

The galaxies with a reliable spectroscopic redshift lie close to the one-to-one line, within the width of their photometric redshift intervals, while the unreliable ones are scattered over the whole redshift range.

## 4. Read in a cutout of the VIS mosaic and overplot the galaxies

We first search for the DR1 MER mosaics that cover our region.
The mosaics come in two collections, wide and deep tiles (see the [MER mosaics notebook](1_Euclid_intro_MER_images.md)); the Deep Survey catalogs are measured on the deep tiles, so we search the deep collection.

```{code-cell} ipython3
image_table = Irsa.query_sia(pos=(coords_cutout, 10 * u.arcsec), collection='euclid_dr1_deep_DpdMerBksMosaic', facility='Euclid')
```

```{note}
The result lists all of the MER mosaic images available at this position: the Euclid VIS, Y, J, H mosaics, and the mosaics of the ground-based surveys that have been put on the same pixel scale. We use the ``facility`` argument above to query for Euclid images only.
```

Various image types are returned, so we keep the science images only:

```{code-cell} ipython3
science_images = image_table[image_table['dataproduct_subtype'] == 'science']
science_images['obs_id', 'energy_bandpassname', 'access_url']
```

We keep the VIS mosaic, which we use for the cutouts and in Firefly below; the first part of its `obs_id` is its tile ID:

```{code-cell} ipython3
vis_mosaic = science_images[science_images['energy_bandpassname'] == 'VIS'][0]
filename = vis_mosaic['access_url']
tileID = vis_mosaic['obs_id'][:9]

print(f'The MER tile ID of this region is {tileID}')
```

The MER mosaics are large, so we read only a 5 x 5 arcmin cutout of the VIS mosaic directly from the cloud, without downloading the full file (see the [MER mosaics notebook](1_Euclid_intro_MER_images.md) for more details on the cloud access).
The `cloud_access` column of the SIA result gives the location of the mosaic in the cloud:

```{code-cell} ipython3
cloud_info = json.loads(vis_mosaic['cloud_access'])
s3_path = f"s3://{cloud_info['aws']['bucket_name']}/{cloud_info['aws']['key']}"
s3_path
```

```{code-cell} ipython3
# Open the file without downloading it, and read only the cutout
hdu = fits.open(s3_path, fsspec_kwargs={'anon': True})
header = hdu[0].header

cutout_image = Cutout2D(hdu[0].section, position=coords_cutout, size=im_cutout, wcs=WCS(header))
cutout_image.data.shape
```

```{tip}
We rely on astropy's WCSAxes framework for making plots of astronomical data in Matplotlib. Please note the usage of the `projection` and `transform` arguments in the code example below.

For more info, please visit the [WCSAxes documentation](https://docs.astropy.org/en/stable/visualization/wcsaxes/index.html).
```

```{code-cell} ipython3
ax = plt.subplot(projection=cutout_image.wcs)

ax.imshow(cutout_image.data, cmap='gray', origin='lower',
          norm=ImageNormalize(cutout_image.data, interval=PercentileInterval(99.9), stretch=LogStretch()))
ax.scatter(result_galaxies['ra'], result_galaxies['dec'], s=36, facecolors='none', edgecolors='red',
           transform=ax.get_transform('world'))

_ = ax.set_title(f'Galaxies with photometric redshifts between {z_min} and {z_max}')
```

## 5. Look at the spectrum and image of one of the brightest galaxies

Let's sort the galaxies with a reliable spectroscopic redshift by their H-band flux, and look at the five brightest ones:

```{code-cell} ipython3
brightest = galaxies[reliable]
brightest.sort(keys='flux_h_unif', reverse=True)
brightest['object_id', 'flux_h_unif', 'phz_median', 'phz_90_int_1', 'phz_90_int_2', 'spe_z', 'spe_z_prob', 'spe_z_rel'][:5]
```

We pick the brightest one:

```{code-cell} ipython3
galaxy = brightest[0]

obj_id = galaxy['object_id']
photo_z = galaxy['phz_median']
spec_z = galaxy['spe_z']
coords_galaxy = SkyCoord(galaxy['ra'], galaxy['dec'], unit='deg')

print(f'Object {obj_id}: photometric redshift {photo_z:.2f}, spectroscopic redshift {spec_z:.3f}')
```

### Retrieve the spectra

We search IRSA's Simple Spectral Access (SSA) service for the DR1 SIR 1D spectra at the position of the galaxy.
A galaxy in the Euclid Deep Fields can have several spectra, from the red and blue grisms and from the deep and wide tiles; as in the [SIR 1D spectra notebook](3_Euclid_intro_1D_spectra.md), we keep the spectrum with the most usable flux bins for each kind of grism.

```{code-cell} ipython3
ssa_result = Irsa.query_ssa(pos=coords_galaxy, radius=1 * u.arcsec, collection='euclid_dr1_DpdSirCombinedSpectra')
spectra = [QTable.read(url) for url in ssa_result['access_url']]

print(f'Found {len(spectra)} spectra')
```

```{code-cell} ipython3
def bad_flux_bins(spectrum):
    # Flux bins with an odd MASK value or a MASK value of 64 or more should not be used
    mask = spectrum['MASK'].value
    return (mask % 2 == 1) | (mask >= 64)


def grism_name(spectrum):
    # Blue grism spectra start below 1.1 micron
    return 'Blue grism' if spectrum['WAVELENGTH'].min() < 1.1 * u.micron else 'Red grism'


best = {}
for spectrum in sorted(spectra, key=lambda spectrum: np.sum(~bad_flux_bins(spectrum)), reverse=True):
    best.setdefault(grism_name(spectrum), spectrum)
list(best)
```

### Plot the spectra

We plot the usable flux bins of the best spectrum of each grism, and mark the strongest emission lines at the spectroscopic redshift, to see how well the photometric redshift agrees with it:

```{code-cell} ipython3
emission_lines = {'Hβ': 4862.7 * u.AA, '[OIII]': 5008.2 * u.AA, 'Hα': 6564.6 * u.AA}
grism_colors = {'Red grism': 'tab:red', 'Blue grism': 'tab:blue'}

fig, ax = plt.subplots(figsize=(12, 4.5))
for name, spectrum in best.items():
    bad = bad_flux_bins(spectrum)
    ax.plot(spectrum['WAVELENGTH'].to(u.micron), np.ma.masked_where(bad, spectrum['SIGNAL']),
            color=grism_colors[name], alpha=0.8, label=name)

for line, rest_wavelength in emission_lines.items():
    observed = (rest_wavelength * (1 + spec_z)).to(u.micron)
    ax.axvline(observed, color='grey', linestyle='--', linewidth=0.8)
    ax.text(observed.value, 0.98, f' {line}', transform=ax.get_xaxis_transform(), rotation=90, va='top', ha='right')

ax.set_title(f'Object {obj_id}: phz_median = {photo_z:.2f}, spe_z = {spec_z:.3f}')
_ = ax.legend(loc='upper left')
```

The Hα and [OIII] lines fall right at the spectroscopic redshift, which agrees well with the photometric redshift.
Not every feature is a line of this galaxy, though: slitless spectra are prone to contamination by neighboring sources, which is the likely origin of the narrow spike near 1.3 micron in the blue grism spectrum.

### Look at the galaxy on the VIS mosaic

We haven't closed the image file above, so we use `Cutout2D` again to cut out a small region around the galaxy:

```{code-cell} ipython3
size_galaxy_cutout = 3 * u.arcsec

cutout_galaxy = Cutout2D(hdu[0].section, position=coords_galaxy, size=size_galaxy_cutout, wcs=WCS(header))

ax = plt.subplot(projection=cutout_galaxy.wcs)
ax.imshow(cutout_galaxy.data, cmap='gray', origin='lower',
          norm=ImageNormalize(cutout_galaxy.data, interval=PercentileInterval(99.9), stretch=AsinhStretch()))
_ = ax.set_title(f'Object {obj_id}')
```

## 6. Load the image in Firefly for interactive exploration

First we initialize the Firefly client as explained [here](https://caltech-ipac.github.io/firefly_client/usage/initializing-vanilla.html) so that we can send our data to Firefly for visualization.

```{code-cell} ipython3
fc = FireflyClient.make_client('https://irsa.ipac.caltech.edu/irsaviewer')
```

### Visualize the FITS image with Firefly

After initializing the client, we can directly use the FITS image file URL to visualize it in Firefly.

Note this can take a while to upload the full MER image.

```{code-cell} ipython3
fc.show_fits_image(filename)
```

### Display the table as an overlay on the FITS image

Since our galaxies table is in memory, we write it to an IO stream (file-like object) for visualizing it in Firefly.
Besides displaying the table, Firefly will also overlay the positions of the galaxies on the FITS image we visualized earlier.
You can click on a table row to inspect the particular galaxy in the image, or click on an overlay marker to inspect the particular galaxy in the table.

```{code-cell} ipython3
tbl_stream = io.BytesIO()
galaxies.write(tbl_stream, format='votable')

fc.show_table(tbl_stream)
```

## References

<!-- TODO: confirm the author list and update the reference (journal, volume, arXiv ID) once the DR1 PHZ paper is public. -->

- Euclid Collaboration: Tucci, M., et al. 2026, *Euclid Data Release 1 (DR1): The PHZ processing*, A&A, in prep.

+++

## About this Notebook

**Updated**: 2026-11-12

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
