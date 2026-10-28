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
  - name: Lauren Aldoroty
---

# Following a simulated supernova through the Roman Time Domain Survey

+++

## Learning Goals:

By the end of this tutorial, you will:

1. learn more about the "observations" that make up the simulated Roman Time Domain Survey (TDS).
2. learn how to find the locations of simulated supernovae in the transient catalog.
3. learn how to ask IRSA which simulated Roman images cover a given position and time.
4. learn how to create aligned cutouts of simulated Roman images.
5. learn how to make an animated gif from these cutouts.

+++

## Introduction

The Roman Time Domain Survey revisits the same patch of sky over and over, which is what makes it possible to watch a transient appear and fade. This notebook picks one simulated Type Ia supernova out of the OpenUniverse2024 transient catalog, collects every Roman image that covers it while it is bright, and stacks those images into a short movie.

The survey stores its images by pointing and detector rather than by sky position, so which files contain a particular supernova is not something you can work out from the file paths. IRSA's Simple Image Access (SIA) service answers exactly that question, and we use it here to assemble the list of images to stack.

If you are new to OpenUniverse2024, the [Quickstart](openuniverse2024_quickstart) tutorial introduces the directory layout, the parquet catalogs, and the image search used below.

+++

## Install and Import required modules

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install astropy matplotlib numpy pandas pyarrow s3fs scipy astroquery hpgeom ipython
```

```{code-cell} ipython3
# Import modules
import warnings
import json

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import hpgeom
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.wcs import WCS, FITSFixedWarning
from astroquery.ipac.irsa import Irsa
from matplotlib import animation
from scipy.ndimage import rotate
from IPython.display import HTML

# Needed to access data in the cloud
import s3fs
s3 = s3fs.S3FileSystem(anon=True)  # create an S3 client

# Filter out the FITSFixedWarning, which is consequenceless and gets thrown every time you deal with a WCS
# in a Roman openuniverse simulated image using astropy.
warnings.simplefilter('ignore', category=FITSFixedWarning)
```

```{code-cell} ipython3
# Point the astroquery IRSA client at the simulated-data services, which are
# separate from the ones serving IRSA's observed data.
Irsa.sia_url = "https://irsa.ipac.caltech.edu/simulated/SIA"
Irsa.tap_url = "https://irsa.ipac.caltech.edu/simulated/TAP"

OU_ROMAN_SIA_COLLECTION = 'simulated_roman_openuniverse2024'
```

## Read in the Observation Sequence File to learn more about the "observations" that make up the simulated Roman Time Domain Survey.

```{code-cell} ipython3
# Read in the (simulated) Observation Sequence File.

BUCKET_NAME = 'nasa-irsa-simulations'
ROMAN_PREFIX = 'openuniverse2024/roman/full'

ROMAN_TDS_PATH = f'{ROMAN_PREFIX}/RomanTDS'
FILENAME = 'Roman_TDS_obseq_11_6_23.fits'
OBSEQ_PATH = f's3://{BUCKET_NAME}/{ROMAN_TDS_PATH}/{FILENAME}'

obseq_hdu = fits.open(OBSEQ_PATH, fsspec_kwargs={"anon": True})
obseq = pd.DataFrame(obseq_hdu[1].data)

print(obseq)
```

## What is the spatial and temporal coverage of the openuniverse2024 Roman TDS?

```{code-cell} ipython3
# Find the ranges of RA, Dec, and date listed in the observation sequence file.

ra_min, dec_min = obseq[['ra','dec']].min()
ra_max, dec_max = obseq[['ra','dec']].max()
mjd_min = obseq['date'].min()
mjd_max = obseq['date'].max()

print("ra_min, ra_max:", ra_min, ra_max)
print("mjd_min, mjd_max:", mjd_min, mjd_max)
```

## Read in the Supernova Analysis (SNANA) file.

The transient catalogs are split by HEALPix sky region (nside=32, RING ordering), with the region index in the filename. We convert the center of the Roman TDS into that index rather than guessing at the file name.

```{code-cell} ipython3
# The Roman Time-Domain Survey is centered near the LSST ELAIS-S1 Deep Drilling Field.
region = hpgeom.angle_to_pixel(32, 9.45, -44.02, lonlat=True, nest=False)

parquet_file = f's3://{BUCKET_NAME}/{ROMAN_PREFIX}/roman_rubin_cats_v1.1.2_faint/snana_{region}.parquet'
transients = pd.read_parquet(parquet_file, filesystem=s3)
```

## Let's find a relatively nearby SN Ia that the survey actually watched go off.

```{code-cell} ipython3
#List the unique models in the SNANA file.
unique_models = pd.Series(transients['model_name']).drop_duplicates().tolist()
unique_models
```

```{code-cell} ipython3
# Most of the models are non SNIa (NON1ASED).
# Choose only the SNIa
sn1a = transients[transients['model_name'] == 'SALT3.NIR_WAVEEXT'] # SNe Ia only.
print('Number of SN1a in SNANA file: ', len(sn1a))
```

```{code-cell} ipython3
# Choose the SNIa that overlap with the spatial and temporal extent of the survey.
ra_mask = np.logical_and(sn1a['ra'] > ra_min, sn1a['ra'] < ra_max)
dec_mask = np.logical_and(sn1a['dec'] > dec_min, sn1a['dec'] < dec_max)
mjd_mask = np.logical_and(sn1a['start_mjd'] > mjd_min, sn1a['end_mjd'] < mjd_max)
all_mask = np.logical_and.reduce((ra_mask,dec_mask,mjd_mask))
covered_sn1a = sn1a[all_mask]
print('Number of SNIa within the survey:', len(covered_sn1a))
```

```{code-cell} ipython3
# Choose the SNIa that are nearby, at redshifts less than 0.5, so they are bright enough
# to stand out clearly against their host galaxy.
nearby_sn1a = covered_sn1a[covered_sn1a['z_CMB'] < 0.5]
print('Number of nearby SNIa:', len(nearby_sn1a))
```

A supernova is only worth animating if Roman happened to be looking at that patch of sky while it was bright. The catalog records the date each one peaks in `peak_mjd`, and the survey visits any given field in bursts rather than continuously, so we pick an object whose peak falls inside a well-visited stretch of the survey.

The cuts above cannot check this for us.
They confirm that a supernova went off somewhere inside the survey's footprint and date range, not that the telescope was pointed at it while it was bright.
That combination is common, and when it happens the epoch window below comes back empty and there is nothing to animate.
So if you change `oid` to explore a different object, check that the dates it was bright (`start_mjd` to `end_mjd`) overlap the dates Roman visited its position (the `t_min` column of the image search below).

```{code-cell} ipython3
# Let's choose SN 20131477, which peaks while its field is being visited regularly.
oid = 20131477
chosen_object = nearby_sn1a[nearby_sn1a['id'] == oid].iloc[0]

ra, dec = chosen_object['ra'], chosen_object['dec']
peak_mjd = chosen_object['peak_mjd']
coord = SkyCoord(ra*u.deg, dec*u.deg)

print(f"SN {oid}: RA={ra:.6f}, Dec={dec:.6f}, z={chosen_object['z_CMB']:.3f}, peaks at MJD {peak_mjd:.1f}")
```

## Ask IRSA which simulated Roman images cover the chosen SNIa.

We hand the position to IRSA's image search, which returns one row per image along with the time it was taken and where the file lives in the cloud.

```{code-cell} ipython3
sia_results = Irsa.query_sia(pos=(coord, 1 * u.arcsec),
                             collection=OU_ROMAN_SIA_COLLECTION)

print(f"Images covering this position: {len(sia_results)}")
```

That covers every band and both Roman surveys, so we narrow it down to the Time Domain Survey images in a single band.

```{code-cell} ipython3
band = 'R062'

is_tds = np.array(['TDS_simple_model' in str(obs_id) for obs_id in sia_results['obs_id']])
is_band = np.char.strip(np.array(sia_results['energy_bandpassname'], dtype=str)) == band

instances = sia_results[is_tds & is_band]
instances.sort('t_min')

print(f"{band} images covering SN {oid}: {len(instances)}")
```

Finally we keep only the epochs around the peak. Frames from years before the explosion would add nothing to the movie beyond download time, so we take a window that starts shortly before the supernova appears and runs until it has faded.

```{code-cell} ipython3
epoch_mjd = np.asarray(instances['t_min'], dtype=float)
in_window = (epoch_mjd >= peak_mjd - 40) & (epoch_mjd <= peak_mjd + 80)

instances = instances[in_window]
epoch_mjd = epoch_mjd[in_window]

print(f"Epochs to animate: {len(instances)}, "
      f"spanning MJD {epoch_mjd.min():.1f} to {epoch_mjd.max():.1f}")
```

The cloud location of each image arrives as a JSON string, which we unpack into an S3 path.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def get_s3_fpath(cloud_access):
    """Extract the S3 URI from the cloud_access JSON string in an image search result."""
    cloud_info = json.loads(cloud_access)['aws']
    return f"s3://{cloud_info['bucket_name']}/{cloud_info['key']}"
```

```{code-cell} ipython3
image_paths = [get_s3_fpath(row['cloud_access']) for row in instances]
image_paths[:3]
```

## Create cutouts of the chosen SNIa.

Roman's focal plane sits at a different angle on the sky at each visit, so the same patch of sky arrives rotated differently in every image. Before the frames can be stacked into a movie we rotate each one so that north points the same way throughout. The rotation angle is recorded in the header, and detectors in every third slot are mounted flipped, which we correct for as well.

```{code-cell} ipython3
#Make the cutouts; this will take a few minutes.
stamps = []
mjd = []
for imgpath, epoch in zip(image_paths, epoch_mjd):
    print(imgpath)
    with fits.open(imgpath, fsspec_kwargs={"anon": True}) as hdu:
        img = hdu[1].data
        header = hdu[0].header
        wcs = WCS(header)
        x, y = wcs.world_to_pixel(coord)

        # Manually rotate the images so they are all aligned.
        CDmat = np.array([header['CD1_1'], header['CD1_2'],
                          header['CD2_1'], header['CD2_2']]).reshape(2,2)

        orientation = header['ORIENTAT']

        # These chips are "flipped".
        if header['SCA_NUM'] % 3 == 0:
            orientation += 180

        # Build rotation matrix.
        CD1_1_rot = np.cos(-orientation*np.pi/180)
        CD1_2_rot = -np.sin(-orientation*np.pi/180)
        CD2_1_rot = np.sin(-orientation*np.pi/180)
        CD2_2_rot = np.cos(-orientation*np.pi/180)

        RotMat = np.array([CD1_1_rot, CD1_2_rot,
                          CD2_1_rot, CD2_2_rot]).reshape(2,2)

        RotMat_inv = np.array([CD1_1_rot, -CD1_2_rot,
                              -CD2_1_rot, CD2_2_rot]).reshape(2,2)

        # Apply rotation to the CDi_j header keywords.
        CDmat_rot = np.dot(CDmat,RotMat_inv)

        # Update header.
        header['CD1_1'], header['CD1_2'] = CDmat_rot[0]
        header['CD2_1'], header['CD2_2'] = CDmat_rot[1]
        header['ORIENTAT'] -= orientation

        # Rotate the image.
        rot_img = rotate(img,angle=orientation,reshape=False,cval=np.nan)

        rot_wcs = WCS(header)

        try:
            # Make cutout around SN Ia location.
            cutout = Cutout2D(rot_img,coord,100,wcs=rot_wcs,mode='partial')
            stamps.append(cutout.data)
            mjd.append(epoch)
        except NoOverlapError:
            pass

print(f"Collected {len(stamps)} cutouts")
```

## Define a module to create an animated gif from a collection of cutouts.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def animate_stamps(stamps, savepath, no_whitespace=True,
                   labels=[],labelxy=(0.05,0.95)):
    """Make an animation of a sequence of image stamps.

    Parameters
    ----------
    stamps : list of `~numpy.ndarray`
        Image stamps, in chronological order.
    savepath : str
        Path to save the gif to.
    no_whitespace : bool, optional
        Drop any stamp that is entirely NaN, along with its label.
    labels : list of str, optional
        Per-frame text, drawn in the corner of each frame.
    labelxy : tuple of float, optional
        Position of the label, in axes fractions.

    Returns
    -------
    `~matplotlib.animation.FuncAnimation`
        The animation, for displaying with playback controls.
    """

    if no_whitespace:
        with_whitespace = np.invert(np.any((np.isnan(np.array(stamps))), axis=(1,2))) # NOTE: Your first axis (first indexing value) should return one stamp. e.g. stamps[0] is the first stamp.
        idx_whitespace = np.where(with_whitespace)[0]
        stamps = np.array(stamps)[idx_whitespace]
        if len(labels) != 0:
            labels = np.array(labels)[idx_whitespace]

    fig, ax = plt.subplots(figsize=(5,5))
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    plt.xticks([])
    plt.yticks([])

    im = ax.imshow(stamps[0], animated=True)

    if len(labels) != 0:
        txt = ax.text(labelxy[0],labelxy[1],labels[0],animated=True,color='white',transform=ax.transAxes,va='top',ha='left')

    def animate(i):
        im.set_array(stamps[i])
        if len(labels) != 0:
            txt.set_text(labels[i])

            return [im] + [txt]
        else:
            return [im]

    writer = animation.PillowWriter()
    anim = animation.FuncAnimation(fig, animate, interval=600, frames=len(stamps))
    anim.save(savepath, writer=writer)

    # Close the figure so the notebook does not also print a static copy of the first frame,
    # and hand the animation back so it can be displayed with playback controls.
    plt.close(fig)
    return anim
```

## Make an animated gif out of the cutouts.

Watch the supernova brighten near the middle of the sequence and fade away again.

```{code-cell} ipython3
savepath = f'SN{oid}.gif'
savepath
anim = animate_stamps(stamps, savepath, labels=[f'MJD {m:.1f}' for m in mjd])
```

The saved gif loops without stopping, which makes a brief brightening hard to follow. Displaying the animation instead gives playback controls: pause it, step one epoch at a time with the arrows, drag the slider to any frame, and read the date in the corner as you go.

```{code-cell} ipython3
HTML(anim.to_jshtml(default_mode='once'))
```

***

## About this notebook

**Updated:** 2026-08-05

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
