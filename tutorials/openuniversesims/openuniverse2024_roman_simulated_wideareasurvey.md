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
  name: python3
---

# Explore OpenUniverse 2024 Data Preview: Simulated Roman Coadds

+++

## Learning Goals

By the end of this tutorial, you will learn how to do the following:

- Identify the row and column that contains a particular ra, dec coordinate.
- Browse the bucket containing the simulated Roman coadds.
- Identify a simulated Roman coadd by filter, row, and column.
- Take a closer look at a simulated coadd file.

```{code-cell} ipython3
# Install libraries if necessary
!pip install astropy numpy s3fs
```

```{code-cell} ipython3
#Import modules
from astropy.io import fits
import numpy as np
import s3fs  # browse buckets
```

## Identify the row and column that contains a particular ra, dec coordinate.

The full simulation covers a 1 deg x 1 deg area centered around RA = 9.5 deg and
Dec = -44.1 deg. This region is divided into 1296 blocks (36 rows and 36 columns),
each 100 arcsec across.

The data preview presented here covers the central 20x20 arcmin, corresponding to 144 blocks (12 rows and 12 columns).

```{code-cell} ipython3
#Choose an RA, DEC of interest.
ra = 9.5981595
dec = -44.2026950

#Centers of data preview blocks. Do not alter.
ra_block_centers = np.array([9.76330352298415, 9.724522605135252, 9.68574158906671,
                    9.646960496603766, 9.608179349571955, 9.56939816979703,
                    9.530616979104877, 9.491835799321422, 9.453054652272561,
                    9.414273559784032, 9.375492543681393, 9.336711625789874])
dec_block_centers = np.array([-44.252584927082495, -44.22480733304182, -44.197029724175756,
                              -44.16925210374898, -44.14147447502621, -44.11369684127218,
                              44.08591920575162, -44.05814157172923, -44.03036394246976,
                              -44.0025863212379, -43.974808711298394, -43.94703111591591])

ra_difference_array = np.absolute(ra_block_centers-ra)
ra_block_centers_index = ra_difference_array.argmin()
closest_ra_center = ra_block_centers[ra_block_centers_index]
ra_dist = 3600. * ra_difference_array[ra_block_centers_index]
if ra_dist > 50:
    print("Chosen ra not covered by OpenUniverse 2024 data preview simulated Roman coadds")
else:
    COLUMN = ra_block_centers_index + 12
    print("COLUMN:", COLUMN)
    print("")

dec_difference_array = np.absolute(dec_block_centers-dec)
dec_block_centers_index = dec_difference_array.argmin()
closest_dec_center = dec_block_centers[dec_block_centers_index]
dec_dist = 3600. * dec_difference_array[dec_block_centers_index]

if dec_dist > 50:
    print("Chosen dec not covered by OpenUniverse 2024 data preview simulated Roman coadds")
else:
    ROW = dec_block_centers_index + 12
    print("ROW:", ROW)
```

## Browse the bucket containing the simulated Roman coadds.

```{code-cell} ipython3
s3 = s3fs.S3FileSystem(anon=True) # create an S3 client

BUCKET_NAME = "nasa-irsa-simulations"
ROMAN_PREFIX = "openuniverse2024/roman/preview"
COADD_PATH = f"{ROMAN_PREFIX}/RomanWAS/images/coadds"

s3.ls(f"{BUCKET_NAME}/{COADD_PATH}")
```

## Identify a Roman simulated coadd by filter, row, and column.

A simulated coadd can be uniquely identified by filter, row, and column.

```{code-cell} ipython3
#Choose a filter, row, and column
FILTER = 'H158' #Filters F184, H158, J129, K213, and Y106 are available in the data preview.
#ROW = 12 #Rows 12-23 are available in the data preview.
#COLUMN = 12 #Columns 12-23 are available in the data preview.

#Construct the coadd filename from the chosen filter, row, and column.
filename_root = f"prod_{FILTER[0]}_{COLUMN}_{ROW}_map.fits"

#Construct the full coadd path from the chosen filter, row, and column.
s3_uri = f"s3://{BUCKET_NAME}/{COADD_PATH}/{FILTER}/Row{ROW}/{filename_root}"

#List this filename to make sure it is found.
s3.ls(s3_uri)
```

## Take a closer look at the simulated coadd file you identified.

```{code-cell} ipython3
#Show a summary of extensions for this file.

with fits.open(s3_uri, fsspec_kwargs={"anon": True}) as hdul:
    hdul.info()
```

The Primary HDU for the coadded image is a cube with 15 layers, i.e., its shape is 1x15x2688x2688. The layers are as follows:

0 = simulated "Science" image (Roman+Rubin simulation, units of e/(0.11 arcsec)^2/exposure)

1 = lab noise: based on dark frames from the April 2023 test, masked at 3 e/p/s. Units: e/(0.11 arcsec)^2/s

2 = GalSim stars, on HEALPix resolution 14 grid, normalized to total flux of 1

3 = noisy stars, on HEALPix resolution 14 grid, normalized to total flux of 2.4e5 e with self-Poisson noise, including 86 e^2/input pixel background variance

4 = stars, on HEALPix resolution 14 grid, total flux 1, but on in only one of the passes (to test transient response)

5 = stars, on HEALPix resolution 14 grid, with total flux that varies by 5% from center to edge of the focal plane (to test what happens when the filter bandpass varies; 5% is highly exaggerated)

6 = GalSim extended objects, on HEALPix resolution 14 grid, right now exponential profiles. The scale radius is log-distributed between 0.125 and 0.500 arcsec, and the ellipticity (g1,g2) is uniformly distributed in the disc of radius 0.5, i.e., g1^2+g2^2<0.5^2.

7,8,9 = same objects as layer 6, but with applied shear of 0.02. The shear orientations are spaced by 60° in tangent vector space, so that in the (g1,g2)-space they are spaced by 120° and can be used for finite differences. Specifically, the directions are: layer 7 -> in East-West direction (shear PA = 270°). (g1,g2) = (0.02,0) layer 8 -> in NNW-SSE direction (shear PA = 330°). (g1,g2) = (-0.02/2,0.02√3/2) layer 9 -> in NNE-SSW direction (shear PA = 30°). (g1,g2) = (-0.02/2,-0.02√3/2)

10 = coadded 1/f noise map, normalized to variance per ln f of 1

11,12,13,14 = coadded white noise maps, different seeds, normalized to variance of 1 in each input pixel

The following HDUs contain additional information:

CONFIG = the configuration file

INDATA = the input images used, as a binary table. The columns are: obsid (int32) -> observation ID sca (int16) -> SCA (1 through 18, inclusive) ra (float64) -> right ascension of pointing center in degrees dec (float64) -> declination of pointing center in degrees pa (float64) -> position angle of pointing in degrees valid (logical) -> input science data file is valid (should be True)

INWEIGHT = the mean input weights for how much each 1.25x1.25 arcsec postage stamp depends on each input exposure. The shape is 1 x Nin x 84 x 84, where Nin is the number of input images listed in INDATA. Note that each postage stamp is 32 output pixels, so 84x32=2688.

If summed on axis 1, this will normally be something close to 1. Deviations of ~10% are common, due to plate scale variations and the normalization issues introduced by diffraction spikes.

INWEIGHTFLAT = a reshape of INWEIGHT suitable for display as a single image in a viewer such as DS9. The contributions from the Nin input exposures are rearranged into a 1 x 84 x (N_in*84) array.

FIDELITY, SIGMA, INWTSUM, EFFCOVER = maps of U_alpha/C, S_alpha, sum_i T_{alpha i}, and the effective coverage as rescaled int16 maps. See Rowe et al. (2011) for details on the definitions of these quantities. The comment in the 'UNIT' keyword indicates how to rescale these.

+++

## About this notebook

- Author: Vandana Desai (Science Lead, IRSA) in conjunction with the IPAC Science Platform team
- Contact: https://irsa.ipac.caltech.edu/docs/help_desk.html
- Updated: 2024-06-10

```{code-cell} ipython3

```
