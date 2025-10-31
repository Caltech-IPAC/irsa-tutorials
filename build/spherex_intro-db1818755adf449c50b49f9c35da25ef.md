---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.17.2
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(spherex-intro)=
# Introduction to SPHEREx Spectral Images

+++

## 1. Learning Goals

- Query for SPHEREx Spectral Image Multi-Extension FITS files (MEFs) that overlap a given coordinate on the sky.
- Read in a SPHEREx Spectral Image MEF using Astropy and understand its structure.
- Interactively visualize a SPHEREx Spectral Image MEF using Firefly.
- Explore each extension. For the Spectral Image extension, this includes understanding the astrometric and spectral WCS systems.

+++

## 2. SPHEREx Overview

SPHEREx is a NASA Astrophysics Medium Explorer mission that launched in March 2025. During its planned two-year mission, SPHEREx will obtain 0.75-5 micron spectroscopy over the entire sky, with deeper data in the SPHEREx Deep Fields. SPHEREx data will be used to:

* **constrain the physics of inflation** by measuring its imprints on the three-dimensional large-scale distribution of matter,
* **trace the history of galactic light production** through a deep multi-band measurement of large-scale clustering,
* **investigate the abundance and composition of water and biogenic ices** in the early phases of star and planetary disk formation.

The community will also mine SPHEREx data and combine it with synergistic data sets to address a variety of additional topics in astrophysics.

More information is available in the [SPHEREx Explanatory Supplement](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp_QR_v1.0.pdf).

+++

## 3. Requirements
The following packages must be installed to run this notebook.

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# %pip install numpy matplotlib astropy astroquery firefly-client
```

## 4. Imports

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

from astroquery.ipac.irsa import Irsa

from firefly_client import FireflyClient

# The time it takes to read SPHEREx files can exceed
# astropy's default timeout limit. Increase it.
from astropy.utils.data import conf
conf.remote_timeout = 120
```

## 5. Search for SPHEREx Spectral Image MEFs that overlap coordinates you are interested in

+++

Define some coordinates of interest.

```{code-cell} ipython3
ra_deg = 304.693508808
dec_deg = 42.4436872991

coord = SkyCoord(ra_deg, dec_deg, unit='deg')
search_radius = 1 * u.arcsec
```

Query IRSA for a list of Spectral Image MEFs that overlap this position. We use the [IRSA module in astroquery](https://astroquery.readthedocs.io/en/latest/ipac/irsa/irsa.html) and the Simple Image Access (SIA) API.

+++

```{tip}
The IRSA SIA collections can be listed using using the ``list_collections`` method, we can filter on the ones containing "spherex" in the collection name:

    Irsa.list_collections(filter='spherex')
```

+++

The collections are documented at [SPHEREx Data Access: Application Program Interfaces (APIs)](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access#application-program-interfaces-apis)
There are currently three collections available for the second Quick Release:

* `'spherex_qr2'` -- Quick Release 2 Spectral Image MEFs that are part of the SPHEREx **Wide Survey**
* `'spherex_qr2_cal'` -- Quick Release 2 **Calibration files**
* `'spherex_qr2_deep'` -- Quick Release 2 Spectral Image MEFs that are part of the SPHEREx **Deep Survey**

```{code-cell} ipython3
results = Irsa.query_sia(pos=(coord, search_radius), collection='spherex_qr2')
```

:::{note}
SPHEREx data are ingested on a weekly basis.
Due to the nature of the ingestion process, availability via SIA will lag on the order of a day.
To avoid this delay, users can access data through the browsable directories or the SPHEREx Data Explorer GUI (see [SPHEREx Data Access](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access)), or do a TAP query as shown in {ref}`spherex-cutouts`.
:::

Each row of the results of your query represents a different spectral image.
Because SPHEREx data will be released on a weekly basis, the number of rows returned will change
depending on when you submit the query.
Let's see how many images are returned today.

```{code-cell} ipython3
len(results)
```

The query results provide a lot of metadata about each spectral image. These columns have standard names as defined by the IVOA. Let's list them:

```{code-cell} ipython3
results.colnames
```

The `'access_url'` column is particularly important because it tells you how to access the data. Let's look at the `'access_url'` value for the first row:

```{code-cell} ipython3
spectral_image_url = results['access_url'][0]
print(spectral_image_url)
```

You can put this URL into a browser to download the file. Or you can work with it in Python, as shown below.

+++

## 6. Examine the header of one of the SPHEREx Spectral Image MEFs

+++

Use Astropy to examine the header of the URL from the previous step.

```{code-cell} ipython3
hdulist = fits.open(spectral_image_url)
hdulist.info()
```

You can see that the Level 2 Spectral Image files are multi-extension FITS files (MEFs) with the following extensions:

* **IMAGE:** Calibrated fluxes for a detector array in scientific units of MJy/sr.
* **FLAGS:** A bitmap with per-pixel status and processing flags.
* **VARIANCE:** A per-pixel estimate of the variance.
* **ZODI:** An model of the zodiacal dust background signal. Note that this is not subtracted from the IMAGE extension.
* **PSF:** An image cube (3D array) where each plane represents a Point Spread Function (PSF) for a cutout in over-sampled pixel space.
* **WCS-WAVE:** Spectral WCS lookup table that maps pixel coordinates to the central wavelength and bandwidth of each pixel.

We will tour each extension in the sections below.

+++

## 7. Visualize a SPHEREx Spectral Image MEF using the Firefly Python Client.

+++

We will use the open-source astronomy data visualization software Firefly visualize SPHEREx data. Firefly has a Python client and it understands the alternative WCS coordinates of SPHEREx, allowing users to see how the wavelength and bandwidth vary with spectral image pixels.

Open a Firefly viewer in a separate browser tab and initialize it.

```{code-cell} ipython3
fc = FireflyClient.make_client(url="https://irsa.ipac.caltech.edu/irsaviewer")

fc.reinit_viewer()
```

Visualize a spectral image MEF by sending its URL to the viewer.

```{code-cell} ipython3
fc.show_fits_image(file_input=spectral_image_url,
             plot_id="spectral_image",
             Title="Spectral Image"
             )
```

Try use the interactive tools in the viewer to explore the data.

+++

## 8. Explore the first extension: IMAGE

The first extension of the MEF is the the calibrated surface brightness flux density in units of MJy/sr, stored as a 2040 x 2040 image. No zodiacal light subtraction is applied.

The SPHEREx focal plane is split with a dichroic to three short-wavelength and three long-wavelength detector arrays. Two focal plane assemblies (FPAs) simultaneously image the sky through a dichroic beam splitter. Each FPA contains three 2K x 2K detector arrays placed behind a set of linear variable filters (LVFs), providing narrow-band response with a band center that varies along one axis of the array. SPHEREx obtains spectra through multiple exposures, placing a given source at multiple positions in the field of view, where it is measured at multiple wavelengths by repointing the spacecraft.

* Band 1: λ= 0.75 - 1.09 µm; R=39
* Band 2: λ= 1.10 - 1.62 µm; R=41
* Band 3: λ= 1.63 - 2.41 µm; R=41
* Band 4: λ= 2.42 - 3.82 µm; R=35
* Band 5: λ= 3.83 - 4.41 µm; R=112
* Band 6: λ= 4.42 - 5.00 µm; R=128

+++

Examine the header of the first extension, printing out select keywords.

```{code-cell} ipython3
spectral_image_header = hdulist[1].header

keywords_to_print = ['EXTNAME', 'NAXIS1', 'NAXIS2', 'BUNIT', 'DETECTOR', 'OBSID', 'DATE', 'PSF_FWHM']

for keyword in keywords_to_print:
    value = spectral_image_header.get(keyword, 'Keyword not found')
    print(f"{keyword}: {value}")
```

We can see that this image was taken with Detector 2, so we can expect the wavelength to vary from 1.10 - 1.62 µm. Try examining the Firefly visualization to confirm this.

+++

Notice that there is more than one WCS!

```{code-cell} ipython3
spectral_image_header['WCSNAME*']
```

The main WCS describes the astrometric registration of the image, including optical distortion parameters, based on the FITS pixel convention starting with 1.

There are two alternative WCS systems:
- WCSNAMEA describes zero-based pixel coordinates.
- WCSNAMEW describes spectral coordinates 'Wavelength' and 'Bandwidth'. This WCS contains a reference to the lookup table in the 'WCS-WAVE' extension.

+++

### 8a. Spatial WCS

+++

Load the Spatial WCS

```{code-cell} ipython3
wcs = WCS(spectral_image_header)
wcs
```

Use standard Astropy methods to resolve the central ra and dec (given by crval1 and crval2 in the header) into image pixel coordinates.

```{code-cell} ipython3
ra, dec = wcs.wcs.crval
x, y = wcs.world_to_pixel(SkyCoord(ra=ra, dec=dec, unit="deg"))
print(ra, dec, x, y)
```

### 8b. Spectral WCS

The use of Linear Variable Filters (LVFs) is a key component of SPHEREx imaging. Each pixel of the detector corresponds to a slightly different wavelength and bandwidth due to the LVF’s gradual variation in spectral transmission across its surface. Contours of constant wavelength are curved due to the method of filter fabrication, so the wavelength-vs-pixel function is inherently two-dimensional.

A compact, approximate representation of the wavelength and bandwidth per pixel is included in each image using the WAVE-TAB lookup-table mechanism defined in the FITS standard.

Below we illustrate how to use Spectral WCS to find an approximate wavelength at each pixel.

```{code-cell} ipython3
# Load the Spectral WCS from the header.
# Note that we need to provide a reference to HDU List, which contains a lookup table.
spectral_wcs = WCS(header=hdulist["IMAGE"].header, fobj=hdulist, key="W")
```

Note: The previous line triggers an Astropy INFO printout,
which implies that the SIP distortion coefficients from the main WCS are preserved in the alternative WCS.
This is because the SIP convention, not formally part of the FITS standard,
is ambiguous as to whether it is meant to apply to 'alternative' (lettered) WCSes in addition to the primary WCS.
See [astropy/astropy#13105](https://github.com/astropy/astropy/issues/13105).

The wavelength per pixel is a property of the detector-filter combination and is independent of optical distortion in the telescope,
and is modeled accordingly in WCS 'W', so we turn the SIP distortion off for this WCS.

```{code-cell} ipython3
# SIP distortions must be turned off for the spectral WCS.
spectral_wcs.sip = None
```

```{code-cell} ipython3
# The standard Astropy methods for converting pixel coordinates to world coordinates can also be used to obtain spectral coordinates.
# Take the pixel coordinates that we determined for the image center and resolve them to the wavelength and bandwidth for that pixel
wl, bw = spectral_wcs.pixel_to_world(x, y)
wl, bw
```

### 8c. How does wavelength vary across the detector?

```{code-cell} ipython3
# Read in the spectral image data.
spectral_image = hdulist["IMAGE"]
spectral_image_data = spectral_image.data

# Get arrays of pixel coordinates for the image.
(y, x) = np.indices(spectral_image.shape)

# Use the spectral WCS to convert these pixel coordinates to spectral coordinates.
spectral_coords = spectral_wcs.pixel_to_world(x, y)

# Break out the two spectral coordinates (wavelength and bandwidth) from the spectral coordinates, and print the results.
wavelength, bandwidth = spectral_coords
print("Wavelength: \n", wavelength)
print("Bandwidth: \n", bandwidth)
```

```{code-cell} ipython3
# Plot the wavelength as a greyscale color map across the image pixels.
plt.imshow(wavelength.value, origin='lower', cmap='gray')
plt.colorbar(label=f"Wavelength [{wavelength.unit}]")

detector = hdulist["IMAGE"].header["DETECTOR"]
plt.title(f"Wavelength - Detector {detector}")
plt.show()
```

As expected, the wavelengths range from approximately 1.1=16 micron for Detector 2. You can see that the longest wavelengths are at the bottom and the shortest wavelengths are at the top. You can verify this is the case in the Firefly visualization, as well.

+++

### 8d. Number of flagged pixels in this image

Now let's take a look at some header keywords that provide information about how many pixels have been flagged during processing.

```{code-cell} ipython3
spectral_image_header['L2 N_*']
```

There are 14 flags in total. Typically, most of the pixels are identified as SOURCE pixels, which are pixels mapped to a known source. The remaining flags are described in Table 8 of the [SPHEREx Explanatory Supplement](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spherex_explanatory_supplement.pdf).

+++

## 9. Explore the second extension: FLAGS

The second extension (FLAGS) provides a bitmap of per-pixel status and processing flags, stored as a 2040 x 2040 image.

Let's take a look at the header of the FLAGS extension, and print out some header keywords of interest.

```{code-cell} ipython3
flags_header = hdulist[2].header

keywords_to_print = ['EXTNAME', 'NAXIS1', 'NAXIS2', 'BUNIT']

for keyword in keywords_to_print:
    value = flags_header.get(keyword, 'Keyword not found')
    print(f"{keyword}: {value}")
```

Conveniently, the definitions of the flags are also provided in the FLAGS header.

```{code-cell} ipython3
flags_header['MP*']
```

Tip: In the Firefly visualization, if you are looking at the first extension (IMAGE), you can open up the layers icon (top right) to enable overlaying a visualization of the flags.

+++

## 10. Explore the third extension: VARIANCE

The third extension of the MEF is the variance of the calibrated surface brightness flux in units of (MJy/sr)^2, stored as a 2,040 x 2,040 image. Let's look at some specific header keywords:

```{code-cell} ipython3
variance_header = hdulist[3].header

keywords_to_print = ['EXTNAME', 'NAXIS1', 'NAXIS2', 'BUNIT']

for keyword in keywords_to_print:
    value = variance_header.get(keyword, 'Keyword not found')
    print(f"{keyword}: {value}")
```

## 11. Explore the fourth extension: ZODI

The fourth extension of the MEF is the modeled zodiacal light background flux in units of MJy/sr, stored as a 2040 x 2040 image. This has not been subtracted from the IMAGE extension. Let's examine some header keywords:

```{code-cell} ipython3
zodi_header = hdulist[4].header

keywords_to_print = ['EXTNAME', 'NAXIS1', 'NAXIS2', 'BUNIT']

for keyword in keywords_to_print:
    value = zodi_header.get(keyword, 'Keyword not found')
    print(f"{keyword}: {value}")
```

## 12. Explore the fifth extension: PSF

The fifth extension of the MEF contains 121 Point-spread functions (PSFs); each PSF is represented as a 101 x 101 image and all 121 are assembled together into a cube. Each of the 121 layers represents a "super-resolution" PSF estimate in a different region (defined by an 11x11 grid) of the detector. Each PSF is a two-dimensional array with size of 101 × 101 pixels. The PSFs are oversampled such that 10 PSF pixels cover the same spatial extent as one spectral image pixel (0.615 arcsec). Let's look at some specific header keywords for the PSF:

```{code-cell} ipython3
psf_header = hdulist[5].header

keywords_to_print = ['EXTNAME', 'NAXIS1', 'NAXIS2', 'NAXIS3', 'BUNIT']

for keyword in keywords_to_print:
    value = psf_header.get(keyword, 'Keyword not found')
    print(f"{keyword}: {value}")
```

## 13. Explore the sixth extension: WCS-WAVE

The sixth extension is a FITS-compliant spectral World Coordinate System (WCS) lookup table that maps spectral image pixel coordinates to central wavelengths and bandwidths. The lookup table consists of 1 row with 3 columns (X, Y, VALUES). X and Y are each arrays defining a grid of control points in spectral image pixel space. For each (X, Y) control point, VALUES defines a two-element array containing the central wavelength and the corresponding bandwidth. Originally adopted to support the unique nature of the SPHEREx LVF filters, this rarely-used part of the FITS standard has yet to be implemented by all readers. The Firefly client we use in this notebook does correctly interpret this lookup table.

Let's look at the header of the WCS-WAVE extension:

```{code-cell} ipython3
wcs_wave_header = hdulist[6].header
wcs_wave_header
```

Unlike the other extensions, this is a binary table. Let's read it iinto an Astropy table.

```{code-cell} ipython3
wcs_wave_table = Table(hdulist[6].data)

print("Number of rows:", len(wcs_wave_table))
print("Column names:", wcs_wave_table.colnames)
```

This table consists of just 1 row with three columns. Let's inspect the columns:

```{code-cell} ipython3
x_array      = wcs_wave_table["X"][0]
y_array      = wcs_wave_table["Y"][0]
wavelengths  = wcs_wave_table["VALUES"][0]

print("X:",  x_array)
print("Y:",y_array)
print("Dimensions of VALUES array:",  wavelengths.shape)
```

## Acknowledgements

- [IPAC-IRSA](https://irsa.ipac.caltech.edu/)

+++

## About this notebook

**Authors:** IPAC Science Platform Team, including Jessica Krick, Troy Raen, Brigitta Sipőcz, Jaladh Singhal,
Andreas Faisst, Shoubaneh Hemmati, Vandana Desai

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions
or problems.

**Updated:** 24 October 2025

**Runtime:** approximately 30 seconds
