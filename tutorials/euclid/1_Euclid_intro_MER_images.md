---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.17.2
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

# Euclid DR1: MER Mosaics

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- understand the basic characteristics of Euclid DR1 MER mosaics;
- understand the difference between wide and deep MER tiles;
- know how to identify the official wide survey mosaic when a tile has several mosaics per band;
- know how to access full MER mosaics from the cloud;
- know how to make smaller cutouts of MER mosaics;
- know how to use matplotlib to plot a grid of cutouts;
- know how to make a color composite from the Euclid VIS and NISP cutouts;
- know how to identify sources in the cutouts and make basic measurements.

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

Among the data products included in DR1 are the Level 2 MER mosaics.
MER (the Merge processing function of the Euclid ground segment) combines the VIS, NISP, and ground-based images into mosaics on a common grid, and detects and measures the sources on them to build the MER catalog.
The mosaics are created from images taken with the Euclid instruments (VIS and NISP), as well as a variety of ground-based surveys (DES, CFIS, Pan-STARRS, and HSC-based surveys) that provide the u, g, r, i, and z bands.
All of the mosaics have been created according to a uniform tiling on the sky, and mapped to a common pixel scale of 0.1 arcsec.
The MER processing for DR1 is described in detail in Euclid Collaboration: Kümmel et al. (2026), see the [References](#references) section below.
This notebook provides a quick introduction to accessing MER mosaics from IRSA.
If you have questions about it, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## Wide and deep MER tiles

DR1 MER mosaics come in two kinds of tiles, which IRSA serves as two separate collections:

* **Wide tiles** (`euclid_dr1_wide_DpdMerBksMosaic`): tile IDs start with `102`, about 32 x 32 arcmin, the same tiling as in Q1.
* **Deep tiles** (`euclid_dr1_deep_DpdMerBksMosaic`): tile IDs start with `101`, about 17 x 17 arcmin, used for fields too deep to process as wide tiles: the Euclid Deep Fields and auxiliary fields such as COSMOS, CDFS, and AEGIS.
  In the Euclid Deep Fields, they stack all of the observations (about 8 times the VIS exposure time of a wide tile), so they reach much fainter sources.

The deep fields are also covered by wide tiles, and because they have been observed many times, these wide tiles have **several mosaics in every band**.
At most one of them is the official wide survey mosaic; the others are each built from a single pass over the deep field.
At our position, only one of the two overlapping wide tiles has its official mosaic in the collection, which the code below selects.
The `PPLANID` keyword in the FITS header tells them apart: the official wide survey mosaic has `-DATA_RELEASE-` in it, the others `DATA_RELEASE_DEEP`.

In the MER catalog (see the [MER catalogs notebook](2_Euclid_intro_MER_catalog.md)), the sources of the `wide_survey` and `wide_mode` tables are measured on the wide mosaics, and those of the `deep_survey` and `deep_mode` tables on the deep mosaics.

+++

## Data volume

Each wide MER mosaic is approximately 1.47 GB, and each deep MER mosaic is approximately 0.42 GB.
Downloading can take some time.
Running the full notebook takes about 20 minutes, most of it spent reading data from the cloud.

+++

## Imports

```{important}
We rely on ``astroquery`` and ``sep`` features that have been recently added, so please make sure you have the respective version v0.4.10 and v1.4 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy 'astropy>=5.3' matplotlib 'astroquery>=0.4.10' 'sep>=1.4' fsspec s3fs
```

```{code-cell} ipython3
import json

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch,  ZScaleInterval, SquaredStretch
from astropy.wcs import WCS
from astropy import units as u

from astroquery.ipac.irsa import Irsa
import s3fs
import sep
```

<!-- TODO: remove this cell when DR1 is public on the production IRSA services. -->

```{code-cell} ipython3
# DR1 is not public yet, so we query IRSA's development server.
Irsa.tap_url = 'https://irsadev.ipac.caltech.edu/TAP'
Irsa.sia_url = 'https://irsadev.ipac.caltech.edu/SIA'
```

## 1. Search for multiwavelength Euclid DR1 MER mosaics in the Euclid Deep Field North
Below we set a search radius of 10 arcsec around a bright edge-on spiral galaxy in the Euclid Deep Field North, an area covered by both deep and wide tiles.

```{code-cell} ipython3
search_radius = 10 * u.arcsec
coord = SkyCoord(273.1670, 64.3687, unit='deg', frame='icrs')
```

Use IRSA's Simple Image Access (SIA) API to search for all Euclid MER mosaics that overlap with the search region you have specified.

The IRSA SIA collections can be listed using the `list_collections` method; we filter on the ones containing "euclid" in the collection name:

```{code-cell} ipython3
Irsa.list_collections(filter='euclid')
```

The two DR1 MER collections are `euclid_dr1_wide_DpdMerBksMosaic` and `euclid_dr1_deep_DpdMerBksMosaic`; the others hold the Q1 and Early Release Observations data and other Euclid data products, such as calibrated frames and spectra.
We search both DR1 MER collections, one for the deep and one for the wide tiles.
For more information, please read the [IRSA User Guide](https://irsa.ipac.caltech.edu/data/Euclid/docs/euclid_archive_at_irsa_user_guide.pdf).

```{code-cell} ipython3
deep_table = Irsa.query_sia(pos=(coord, search_radius), collection='euclid_dr1_deep_DpdMerBksMosaic')
wide_table = Irsa.query_sia(pos=(coord, search_radius), collection='euclid_dr1_wide_DpdMerBksMosaic')
```

These tables list all MER mosaic images available in this search position. These mosaics include the Euclid VIS, Y, J, H images, as well as ground-based telescopes which have been put on the same pixel scale. For more information, see the [Euclid documentation at IPAC](https://euclid.caltech.edu/page/euclid-faq-tech/).

Note that there are various image types are returned as well, we filter out the `science` images from these:

```{code-cell} ipython3
deep_science_images = deep_table[deep_table['dataproduct_subtype'] == 'science']
wide_science_images = wide_table[wide_table['dataproduct_subtype'] == 'science']
```

The tile ID is the first part of the `obs_id`: deep tile IDs start with `101` and wide tile IDs with `102`.
Adjacent tiles overlap, and our position happens to lie in the overlap of two deep tiles and of two wide tiles:

```{code-cell} ipython3
for name, images in [('Deep', deep_science_images), ('Wide', wide_science_images)]:
    tile_ids, n_vis = np.unique([obs_id[:9] for obs_id in images[images['energy_bandpassname'] == 'VIS']['obs_id']], return_counts=True)
    print(f'{name} tiles and their number of VIS mosaics:', dict(zip(tile_ids.tolist(), n_vis.tolist())))
```

We see that each deep tile has a single VIS mosaic, while the wide tiles have several.
Among the two deep tiles, we keep the one whose center is closest to our position:

```{code-cell} ipython3
tile_centers = SkyCoord(deep_science_images['s_ra'], deep_science_images['s_dec'], unit='deg')
deep_tileID = deep_science_images['obs_id'][np.argmin(tile_centers.separation(coord))][:9]
deep_science_images = deep_science_images[[obs_id.startswith(deep_tileID) for obs_id in deep_science_images['obs_id']]]

print(f'Deep MER tile ID: {deep_tileID}')
deep_science_images
```

## 2. Locate the mosaics in the cloud

IRSA hosts the DR1 MER mosaics in an AWS S3 bucket.
The `cloud_access` column of the search results gives the bucket name and file location of each mosaic, so let's extract the S3 file path from it:

```{code-cell} ipython3
def get_s3_fpath(cloud_access):
    cloud_info = json.loads(cloud_access)  # converts str to dict
    bucket_name = cloud_info['aws']['bucket_name']
    key = cloud_info['aws']['key']

    return f'{bucket_name}/{key}'


deep_science_images['s3_fpath'] = [get_s3_fpath(row['cloud_access']) for row in deep_science_images]
wide_science_images['s3_fpath'] = [get_s3_fpath(row['cloud_access']) for row in wide_science_images]
```

A few files listed in the search results may not be available in the bucket yet, so we check that each file exists and skip the ones that do not:

```{code-cell} ipython3
s3 = s3fs.S3FileSystem(anon=True)


def drop_missing(images):
    exists = np.array([s3.exists(path) for path in images['s3_fpath']], dtype=bool)
    for obs_id in images['obs_id'][~exists]:
        print(f'Not in the bucket, skipping: {obs_id}')
    return images[exists]


deep_science_images = drop_missing(deep_science_images)
wide_science_images = drop_missing(wide_science_images)
```

## 3. Identify the official wide survey VIS mosaic

Because our position is in the Euclid Deep Field North, which has been observed many times, each wide tile has several VIS mosaics: the official wide survey mosaic and mosaics built from single passes over the deep field (see [Wide and deep MER tiles](#wide-and-deep-mer-tiles) above).
To find the official wide survey VIS mosaic, we read the `PPLANID` keyword from the header of each wide VIS mosaic returned by our search.
Reading only the header with a small block size avoids downloading the large image data.

```{code-cell} ipython3
def get_processing_plan(s3_fpath):
    with s3.open(s3_fpath, block_size=28800) as f:
        return fits.Header.fromfile(f).get('PPLANID', '')


wide_vis_images = wide_science_images[wide_science_images['energy_bandpassname'] == 'VIS']
wide_vis_images['plan'] = [get_processing_plan(path) for path in wide_vis_images['s3_fpath']]

wide_vis_images['official'] = ['-DATA_RELEASE-' in plan for plan in wide_vis_images['plan']]
official_vis = wide_vis_images[wide_vis_images['official']]
if len(official_vis) == 0:
    raise ValueError('No official wide survey VIS mosaic was found at this position; '
                     "try another position, or look at wide_vis_images['obs_id', 'plan'].")
wide_tileID = official_vis['obs_id'][0][:9]

wide_vis_images['obs_id', 'plan', 'official']
```

## 4. Retrieve a Euclid DR1 MER mosaic image in the VIS bandpass

+++

### Let's first look at one example full image, the deep VIS image

We open the full deep VIS mosaic directly from the S3 bucket. Note this file is about 0.42 GB.

```{code-cell} ipython3
filename = deep_science_images[deep_science_images['energy_bandpassname'] == 'VIS']['s3_fpath'][0]
print(filename)

hdu_mer_irsa = fits.open(f's3://{filename}', fsspec_kwargs={'anon': True})
hdu_mer_irsa.info()

header_mer_irsa = hdu_mer_irsa[0].header
```

If you would like to save the MER mosaic to disk, uncomment the following cell.
Please also define a suitable download directory; by default it will be `data` at the same location as your notebook.

```{code-cell} ipython3
# download_path = 'data'
# hdu_mer_irsa.writeto(os.path.join(download_path, 'MER_image_VIS.fits'), overwrite=True)
```

Have a look at the header information for this image.

```{code-cell} ipython3
header_mer_irsa
```

Lets extract just the primary image. This reads the full image from the cloud, so it may take a minute.

```{code-cell} ipython3
im_mer_irsa = hdu_mer_irsa[0].data
hdu_mer_irsa.close()

print(im_mer_irsa.shape)
```

Due to the large field of view of the MER mosaic, let's cut out a smaller section (2'x2') around our search position to inspect the image

```{code-cell} ipython3
im_2arcmin = Cutout2D(im_mer_irsa, position=coord, size=2 * u.arcmin, wcs=WCS(header_mer_irsa)).data

plt.imshow(im_2arcmin, cmap='gray', origin='lower',
           norm=ImageNormalize(im_2arcmin, interval=PercentileInterval(99.9), stretch=AsinhStretch()))
colorbar = plt.colorbar()
```

Let's look at the whole tile, plotting only every 20th pixel to keep it fast.

```{code-cell} ipython3
# Plot every 20th pixel, which is enough to see the whole tile and fast to display
im_binned = im_mer_irsa[::20, ::20]
plt.imshow(im_binned, cmap='gray', origin='lower',
           norm=ImageNormalize(im_binned, interval=PercentileInterval(99.9), stretch=AsinhStretch()))
plt.colorbar()
plt.title(f'Deep VIS mosaic of tile {deep_tileID} (every 20th pixel)')
plt.show()
```

## 5. Create multiwavelength Euclid DR1 MER cutouts of a region of interest

+++

```{note}
We'd like to take a look at the multiwavelength images of our object, but the full MER mosaics are very large, so we will inspect the multiwavelength cutouts.
```

Create short labels with the instrument and filter name so we can add them to the plots, and sort the bands by wavelength.

```{code-cell} ipython3
# Short labels for the plots, e.g. 'HSC G' or 'NISP H'
short_names = {'Hyper Suprime-Cam': 'HSC'}
deep_science_images['filters'] = [band if band == 'VIS' else f"{short_names.get(inst, inst)} {band}"
                             for inst, band in deep_science_images['instrument_name', 'energy_bandpassname']]

# Sort the bands by wavelength (u, g, r, VIS, i, z from the ground and Euclid, then the NISP bands)
band_order = ['U', 'G', 'R', 'R2', 'VIS', 'I', 'I2', 'Z', 'Y', 'J', 'H']
deep_science_images['band_order'] = [band_order.index(band) if band in band_order else len(band_order)
                                for band in deep_science_images['energy_bandpassname']]
deep_science_images.sort('band_order')
deep_science_images['filters']
```

## 6. Define cutout parameters for a smaller region of interest

Each cutout is read directly from the cloud, without downloading the full mosaic, but reading the cutouts in all of the bands can still take a few minutes.
Larger cutouts take longer.
A cutout cannot be larger than the tile itself (about 17 arcmin for deep tiles and 32 arcmin for wide tiles); near the edge of a tile it is trimmed to the part that overlaps the tile.

```{code-cell} ipython3
######################## User defined section ############################
# Set the image cutout size
im_cutout = 1.0 * u.arcmin

# Set the cutout center coordinates, here our search position
coords_cutout = coord

##########################################################################


def make_cutout(s3_fpath):
    # Read in only the cutout of the image that you want, without downloading the full file
    with fits.open(f's3://{s3_fpath}', fsspec_kwargs={'anon': True}) as hdu:
        header = hdu[0].header
        cutout_data = Cutout2D(hdu[0].section, position=coords_cutout, size=im_cutout, wcs=WCS(header))

    # Define a new fits file based on this smaller cutout, with accurate WCS based on the cutout size
    new_hdu = fits.PrimaryHDU(data=cutout_data.data, header=header)
    new_hdu.header.update(cutout_data.wcs.to_header())
    return new_hdu


# Iterate through each filter
cutout_list = [make_cutout(s3_fpath) for s3_fpath in deep_science_images['s3_fpath']]

# Combine all cutouts into a single HDUList and display information
final_hdulist = fits.HDUList(cutout_list)
final_hdulist.info()
```

## 7. Visualize multiwavelength Euclid DR1 MER cutouts

We need to determine the number of images for the grid layout, and then plot each cutout.

```{code-cell} ipython3
num_images = len(final_hdulist)
columns = 4
rows = -(-num_images // columns)

fig, axes = plt.subplots(rows, columns, figsize=(4 * columns, 4 * rows), subplot_kw={'projection': WCS(final_hdulist[0].header)})
axes = axes.flatten()

for idx, (ax, filt) in enumerate(zip(axes, deep_science_images['filters'])):
    image_data = final_hdulist[idx].data
    norm = ImageNormalize(image_data, interval=PercentileInterval(99.9), stretch=AsinhStretch())
    ax.imshow(image_data, cmap='gray', origin='lower', norm=norm)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.text(0.05, 0.05, filt, color='white', fontsize=14, transform=ax.transAxes, va='bottom', ha='left')

# Remove empty subplots if any
for ax in axes[num_images:]:
    fig.delaxes(ax)

plt.tight_layout()
plt.show()
```

### Make a color composite

Because all MER mosaics share the same pixel grid, we can combine the deep cutouts directly into a false-color RGB image, with NISP H as red, NISP J as green, and VIS as blue.

```{code-cell} ipython3
def scale_band(filt):
    image_data = final_hdulist[np.where(deep_science_images['filters'] == filt)[0][0]].data
    return ImageNormalize(image_data, interval=PercentileInterval(99.5), stretch=AsinhStretch(), clip=True)(image_data)


rgb = np.dstack([scale_band('NISP H'), scale_band('NISP J'), scale_band('VIS')])

plt.figure(figsize=(6, 6))
plt.imshow(rgb, origin='lower')
plt.title('Deep tile: H (red), J (green), VIS (blue)')
plt.axis('off')
plt.show()
```

### Compare the deep and wide VIS cutouts

Let's plot the same region in all of the VIS mosaics of the wide tile: first the official wide survey mosaic, then the mosaics built from single passes over the deep field, and finally the deep mosaic.
All of them are in the same units, so we display them with the same linear stretch from -2 to +5 times the sky noise of the official wide mosaic.
This makes the noise visible: the wide mosaics look grainy, while the deep mosaic is smoother and reveals more faint sources.
The sky noise of each mosaic, measured with sigma clipping, is given in the panel titles.
The official wide survey mosaic is itself built from a single wide survey visit, so its noise is similar to that of the single-pass mosaics; small differences between them come from the observing conditions of each visit.
Artifacts such as satellite trails can appear in any single visit.
This cell reads a cutout from each of the wide VIS mosaics of the tile, so it can take several minutes to run.

```{code-cell} ipython3
filt_index = np.where(deep_science_images['filters'] == 'VIS')[0][0]
deep_vis = final_hdulist[filt_index].data

# Wide VIS mosaics of the tile with the official mosaic, official first
wide_tile_vis = wide_vis_images[[obs_id.startswith(wide_tileID) for obs_id in wide_vis_images['obs_id']]]
wide_tile_vis.sort('official', reverse=True)

images = [make_cutout(s3_fpath).data for s3_fpath in wide_tile_vis['s3_fpath']] + [deep_vis]
titles = (['Wide (official)'] + [f'Wide (single pass {i})' for i in range(1, len(wide_tile_vis))]
          + [f'Deep tile {deep_tileID}'])

noise = [sigma_clipped_stats(image_data)[2] for image_data in images]
norm = ImageNormalize(vmin=-2 * noise[0], vmax=5 * noise[0])

columns = 4
rows = -(-len(images) // columns)
fig, axes = plt.subplots(rows, columns, figsize=(4 * columns, 4 * rows))
axes = axes.flatten()
for ax, image_data, title, sigma in zip(axes, images, titles, noise):
    ax.imshow(image_data, cmap='gray', origin='lower', norm=norm)
    ax.set_title(f'{title}\nnoise = {sigma:.1e}')
    ax.set_xticks([])
    ax.set_yticks([])
for ax in axes[len(images):]:
    fig.delaxes(ax)
plt.tight_layout()
plt.show()
```

## 8. Identify and measure sources in Euclid DR1 MER cutouts with sep

We extract sources on the deep VIS cutout.

```{code-cell} ipython3
img1 = deep_vis
```

### Extract some sources from the cutout using sep (python package based on source extractor)

Following the sep tutorial, first create a background for the cutout
https://sep.readthedocs.io/en/stable/tutorial.html

Need to do some initial steps (swap byte order) with the cutout to prevent sep from crashing. Then create a background model with sep, inspect it and its rms, and subtract it from the cutout.

```{code-cell} ipython3
img2 = img1.byteswap().view(img1.dtype.newbyteorder())
c_contiguous_data = np.array(img2, dtype=np.float32)

# Use a background mesh (bw, bh) larger than our galaxy (about 250 pixels long), so its light is not absorbed into the background
bkg = sep.Background(c_contiguous_data, bw=300, bh=300)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, image_data, title in zip(axes, [bkg.back(), bkg.rms()], ['Background', 'Background rms']):
    im = ax.imshow(image_data, interpolation='nearest', cmap='gray', origin='lower')
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
plt.show()

data_sub = img2 - bkg
```

### Source extraction via sep

```{code-cell} ipython3
######################## User defined section ############################

# Sigma threshold to consider this a detection above the global RMS
threshold= 3

# Minimum number of pixels required for an object. Default is 5.
minarea_0=2

# Minimum contrast ratio used for object deblending. Default is 0.005. To entirely disable deblending, set to 1.0.
deblend_cont_0= 0.005

flux_threshold= 0.01
##########################################################################


sources = sep.extract(data_sub, threshold, err=bkg.globalrms, minarea=minarea_0, deblend_cont=deblend_cont_0)
sources_thr = sources[sources['flux'] > flux_threshold]
print("Found", len(sources_thr), "objects above flux threshold")
```

## 9. Review detected sources on the VIS cutout


We plot the VIS cutout with the sources detected overplotted with a red ellipse

```{code-cell} ipython3
fig, ax = plt.subplots()
m, s = np.mean(data_sub), np.std(data_sub)
im = ax.imshow(data_sub, cmap='gray', origin='lower', norm=ImageNormalize(img2, interval=ZScaleInterval(), stretch=SquaredStretch()))

# Plot an ellipse for each object detected with sep

for i in range(len(sources_thr)):
    e = Ellipse(xy=(sources_thr['x'][i], sources_thr['y'][i]),
                width=6*sources_thr['a'][i],
                height=6*sources_thr['b'][i],
                angle=sources_thr['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
```

## References

<!-- TODO: update the Kümmel et al. reference (journal, volume, arXiv ID) and verify the DR1 DOI once the paper and DR1 are public. -->

- Euclid Collaboration: Kümmel, M., et al. 2026, *Euclid Data Release (DR1): The MERge processing*, A&A, in prep.
- Euclid Collaboration: Romelli, E., Kümmel, M., Dole, H., et al. 2026, A&A, 711, A4 (the Q1 MER processing)
- Euclid Data Release 1, 2026, [doi:10.57780/esa-d079c51](https://doi.org/10.57780/esa-d079c51)

## About this Notebook

**Updated**: 2026-11-12

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
