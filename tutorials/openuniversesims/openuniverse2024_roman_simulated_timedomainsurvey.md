---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: science_demo
  language: python
  name: conda-env-science_demo-py
---

# Explore OpenUniverse 2024 Data Preview: Simulated Roman Time Domain Survey Data

+++

## Learning Goals:

By the end of this tutorial, you will:

1. learn more about the "observations" that make up the simulated Roman TDS preview.
2. learn how to find the locations of simulated supernovae in the preview data.
3. learn how to create aligned cutouts of simulated Roman images.
4. learn how to make an animated gif from these cutouts.

+++

# Import required modules

```{code-cell} ipython3
# Install libraries if necessary
!pip install astropy matplotlib numpy pandas pyarrow s3fs scipy
```

```{code-cell} ipython3
# Import modules
import warnings

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning
from matplotlib import animation
from scipy.ndimage import rotate

# Needed to access data in the cloud
import s3fs
s3 = s3fs.S3FileSystem(anon=True) # create an S3 client

# Filter out the FITSFixedWarning, which is consequenceless and gets thrown every time you deal with a WCS
# in a Roman openuniverse simulated image using astropy.
warnings.simplefilter('ignore',category=FITSFixedWarning)
```

## Define a module to get the date (mjd) of a particular pointing.

```{code-cell} ipython3
def get_mjd(pointing,
            obseq_path=f's3://nasa-irsa-simulations/openuniverse2024/roman/preview/RomanTDS/Roman_TDS_obseq_11_6_23.fits'):

    """
    Retrieve MJD of a given pointing.

    :param pointing: Pointing ID.
    :type pointing: int
    :param obseq_path: Path to obseq file Roman_TDS_obseq_11_6_23.fits.
    :type obseq_path: str, optional
    :return: MJD of specified pointing.
    :rtype: float
    """

    with fits.open(obseq_path, fsspec_kwargs={"anon": True}) as obs:
        obseq = Table(obs[1].data)
    mjd = float(obseq['date'][int(pointing)])

    return mjd
```

## Define a module to create an animated gif from a collection of cutouts.

```{code-cell} ipython3
def animate_stamps(stamps,savepath,no_whitespace=True,
                   labels=[],labelxy=(0.05,0.95),
                   **kwargs):
    """
    Make an animation of a sequence of image stamps.

    :param stamps: Must be in chronological order.
    :type stamps: List of stamps from get_stamps or get_object_instances.
    :param savepath: Path to save gif.
    :type savepath: str
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
        txt = ax.text(labelxy[0],labelxy[1],labels[0],animated=True,color='white',transform=ax.transAxes,va='top',ha='left',**kwargs)

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
```

## Read in the Observation Sequence File to learn more about the "observations" that make up the simulated Roman Time Domain Survey.

```{code-cell} ipython3
# Read in the (simulated) Observation Sequence File.

BUCKET_NAME = 'nasa-irsa-simulations'
ROMAN_PREFIX = 'openuniverse2024/roman/preview'

ROMAN_TDS_PATH = f'{ROMAN_PREFIX}/RomanTDS'
FILENAME = 'Roman_TDS_obseq_11_6_23.fits'
OBSEQ_PATH = f's3://{BUCKET_NAME}/{ROMAN_TDS_PATH}/{FILENAME}'

obseq_hdu = fits.open(OBSEQ_PATH, fsspec_kwargs={"anon": True})
obseq = pd.DataFrame(obseq_hdu[1].data)

print(obseq)
```

## What is the spatial and temporal coverage of the openuniverse2024 Roman TDS data preview?

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

```{code-cell} ipython3
parquet_file = f's3://{BUCKET_NAME}/{ROMAN_PREFIX}/roman_rubin_cats_v1.1.2_faint/snana_10307.parquet'
transients = pd.read_parquet(parquet_file, filesystem=s3)
```

## Let's find a relatively nearby SN Ia that lies within the region of the data preview. 

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
# Choose the SNIa that overlap with the spatial extent of the OpenUniverse2024 Roman TDS data preview.
ra_mask = np.logical_and(sn1a['ra'] > ra_min, sn1a['ra'] < ra_max)
dec_mask = np.logical_and(sn1a['dec'] > dec_min, sn1a['dec'] < dec_max)
mjd_mask = np.logical_and(sn1a['start_mjd'] > mjd_min, sn1a['end_mjd'] < mjd_max)
all_mask = np.logical_and.reduce((ra_mask,dec_mask,mjd_mask))
preview_sn1a = sn1a[all_mask]
print('Number of SNIa in OpenUniverse2024 data preview:', len(preview_sn1a))
```

```{code-cell} ipython3
# Choose the SNIa in the data preview that are nearby, at redshifts less than 0.7.
nearby_preview_sn1a = preview_sn1a[preview_sn1a['z_CMB'] < 0.7]
print('Number of nearby SNIa in OpenUniverse2024 data preview:', len(nearby_preview_sn1a))
```

```{code-cell} ipython3
# Let's choose SN 20000808.
oid = 20000808
chosen_object = nearby_preview_sn1a[nearby_preview_sn1a['id'] == oid]
ra = chosen_object.get('ra')
dec = chosen_object.get('dec')
ra, dec = 9.619282, -44.313894
coord = SkyCoord(ra*u.deg, dec*u.deg)
```

## Read in the auxiliary file that lists the simulated Roman TDS images covering the chosen SNIa.

```{code-cell} ipython3
# The auxiliary file contains all the images this thing is in.
# If you need to download this file, see https://irsa.ipac.caltech.edu/docs/notebooks/.
csvfile = './openuniverse2024_roman_demodata_20000808_instances.csv'
instances = pd.read_csv(csvfile, usecols=['filter','pointing','sca'])
instances
```

## Create cutouts of the chosen SNIa in the band of your choice.

```{code-cell} ipython3
band = 'R062'
instances = instances[instances['filter'] == band]
```

```{code-cell} ipython3
#Make the cutouts; this will take a couple of minutes.
stamps = []
mjd = []
for i, row in enumerate(instances.itertuples()):
    band, pointing, sca = row[1], row[2], row[3]
    imgpath = f's3://{BUCKET_NAME}/{ROMAN_TDS_PATH}/images/simple_model/{band}/{pointing}/Roman_TDS_simple_model_{band}_{pointing}_{sca}.fits.gz'
    print(imgpath)
    with fits.open(imgpath, fsspec_kwargs={"anon": True}) as hdu:
        img = hdu[1].data
        header = hdu[0].header
        wcs = WCS(header)
        x, y = wcs.world_to_pixel(coord)

        # Manually rotate the images so they are all aligned.
        CDmat = np.array([header['CD1_1'], header['CD1_2'],
                          header['CD2_1'], header['CD2_2']]).reshape(2,2)

        orientation = hdu[0].header['ORIENTAT']

        # These chips are "flipped".
        if sca % 3 == 0:
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
        hdu[1].data = rot_img

        rot_wcs = WCS(header)

        try:
            # Make cutout around SN Ia location.
            cutout = Cutout2D(rot_img,coord,100,wcs=rot_wcs,mode='partial')
            stamps.append(cutout.data)
            mjd.append(get_mjd(pointing))
        except NoOverlapError:
            pass
```

## Make an animated gif out of the cutouts.

```{code-cell} ipython3
savepath = f'SN{oid}.gif'
savepath
animate_stamps(stamps,savepath,labels=mjd)
```

![animated gif](SN20000808.gif)

+++

## About this notebook

- Author: Lauren Aldoroty (laurenaldoroty@gmail.com) with minor subsequent modifications to match repository style
- Contact: https://irsa.ipac.caltech.edu/docs/help_desk.html
- Updated: 2024-06-10

```{code-cell} ipython3

```
