---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.17.3
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
authors:
  - name: Vandana Desai
  - name: Jessica Krick
  - name: Andreas Faisst
  - name: Troy Raen
  - name: Brigitta Sipőcz
---

(spherex-cutouts)=
# Download a collection of SPHEREx Spectral Image cutouts as a multi-extension FITS file

## 1. Learning Goals

- Perform a query for the list of SPHEREx Spectral Image Multi-Extension FITS files (MEFs) that overlap a given coordinate.
- Retrieve cutouts for every entry in this list and package the cutouts as a new MEF.
- Learn how to use parallel or serial processing to retrieve the cutouts.

## 2. SPHEREx Overview

SPHEREx is a NASA Astrophysics Medium Explorer mission that launched in March 2025.
During its planned two-year mission, SPHEREx will obtain 0.75-5 micron spectroscopy over the entire sky, with deeper data in the SPHEREx Deep Fields.
SPHEREx data will be used to:

* **constrain the physics of inflation** by measuring its imprints on the three-dimensional large-scale distribution of matter,
* **trace the history of galactic light production** through a deep multi-band measurement of large-scale clustering,
* **investigate the abundance and composition of water and biogenic ices** in the early phases of star and planetary disk formation.

The community will also mine SPHEREx data and combine it with synergistic data sets to address a variety of additional topics in astrophysics.

More information is available in the [SPHEREx Explanatory Supplement](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp_QR.pdf).

## 3. Imports

The following packages must be installed to run this notebook.

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy astropy astroquery matplotlib
```

```{code-cell} ipython3
import concurrent.futures
import http.client
import time
import urllib.error

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import quantity_support
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa

# Suppress logging temporarily to prevent astropy
# from repeatedly printing out warning notices related to alternate WCSs
import logging
logging.getLogger('astropy').setLevel(logging.ERROR)

# The time it takes to read SPHEREx files can exceed
# astropy's default timeout limit. Increase it.
from astropy.utils.data import conf
conf.remote_timeout = 120

# Let matplotlib recognize astropy Quantities and auto-label axes with their units.
quantity_support()
```

## 4. Specify inputs and outputs

Specify a right ascension, declination, cutout size, and a SPHEREx bandpass (e.g., 'SPHEREx-D2').

In this example, we are creating cutouts of the Pinwheel galaxy (M101) for the SPHEREx detector D2.

```{code-cell} ipython3
# Choose a position.
coord = SkyCoord(ra=210.80227 * u.degree, dec=54.34895 * u.degree, frame="icrs")

# Choose a cutout size.
size = 0.1 * u.degree

# Choose the bandpass of interest.
bandpass = 'SPHEREx-D2'

# Choose an output filename root for the output MEF file.
output_filename = 'spherex_cutouts_mef.fits'
```

## 5. Query IRSA for a list of cutouts that satisfy the criteria specified above.

Here we show how to use `astroquery`'s IRSA TAP interface to retrieve all images that overlap with the position defined above.
This query will retrieve a table of URLs that link to the MEF cutouts.
Each row in the table corresponds to a single cutout and includes the data access URL and an observation timestamp.
The results are sorted from oldest to newest.

```{code-cell} ipython3
# Define a query that will search the appropriate SPHEREx metadata tables
# for spectral images that cover the chosen coordinates and match the
# specified bandpass. Return the cutout data access URL and the time of observation.
# Sort by observation time.
query = f"""
SELECT
    'https://irsa.ipac.caltech.edu/' || a.uri || '?center={coord.ra.deg},{coord.dec.deg}d&size={size.value}' AS uri,
    p.time_bounds_lower
FROM spherex.artifact a
JOIN spherex.plane p ON a.planeid = p.planeid
WHERE 1 = CONTAINS(POINT('ICRS', {coord.ra.deg}, {coord.dec.deg}), p.poly)
        AND p.energy_bandpassname = '{bandpass}'
ORDER BY p.time_bounds_lower
"""

# Execute the query and return as an astropy Table.
t1 = time.time()
results = Irsa.query_tap(query)
print("Time to do TAP query: {:2.2f} seconds.".format(time.time() - t1))
print("Number of images found: {}".format(len(results)))
```

:::{note}
SPHEREx data are also available via SIA which can provide a simpler interface for many queries, as demonstrated in {ref}`spherex-intro`.
An advantage of the method shown above is that it provides access to data immediately after ingestion (which occurs weekly) and is not subject to the same ~1 day delay as SIA.
:::

## 6. Define a function that processes a list of SPHEREx Spectral Image Cutouts

This function takes in a row of the catalog that we created above and does the following:

- It downloads the cutout
- It computes the wavelength of the center pixel of the cutout (in micro-meters)
- It keeps a copy of the cutout's original primary header and `VERSION` keyword, so that metadata from IBE is not lost when many cutouts are later combined into one MEF
- It combines the image HDUs into a new HDU and adds it to the table row.

Note that the values of the rows are being added in place.

```{code-cell} ipython3
def process_cutout(row, coord, cache):
    '''
    Downloads the cutouts given in a row of the table including all SPHEREx images overlapping with a position.

    Parameters:
    ===========

    row : astropy.table row
        Row of a table that will be changed in place by this function. The table
        is created by the SQL TAP query.
    coord : astropy.coordinates.SkyCoord
        Coordinates (same as used for the TAP query).
    cache : bool
        If set to `True`, the output of cached and the cutout processing will run faster next time.
        Turn this feature off by setting ``cache=False``.
    '''

    with fits.open(row["uri"], cache=cache) as hdulist:
        # There are seven HDUs:
        # 0 contains minimal metadata in the header and no data.
        # 1 through 6 are: IMAGE, FLAGS, VARIANCE, ZODI, PSF, WCS-WAVE
        header = hdulist[1].header

        # Compute pixel coordinates corresponding to cutout position.
        spatial_wcs = WCS(header)
        x, y = spatial_wcs.world_to_pixel(coord)

        # Compute wavelength at cutout position.
        spectral_wcs = WCS(header, fobj=hdulist, key="W")
        spectral_wcs.sip = None
        wavelength, bandwidth = spectral_wcs.pixel_to_world(x, y)
        row["central_wavelength"] = wavelength.to(u.micrometer).value

        # IBE's primary header carries per-file metadata (e.g. VERSION) that isn't repeated
        # in the extension headers. Keep a copy of it and the version string so they aren't
        # silently dropped when we build the combined MEF below.
        row["primary_header"] = hdulist[0].header.copy()
        row["version"] = hdulist[0].header.get("VERSION", "")

        # Collect the HDUs for this cutout and append the row's cutout_index to the EXTNAME.
        hdus = []
        for hdu in hdulist[1:]:  # skip the primary header
            hdu.header["EXTNAME"] = f"{hdu.header['EXTNAME']}{row['cutout_index']}"
            hdus.append(hdu.copy())  # Copy so the data is available after the file is closed
        row["hdus"] = hdus
```

We provide a small convenience wrapper around `process_cutout` that is used in the rest of this notebook.
It catches transient read errors and simply skips those cutouts, which is sufficient for this demonstration.
For science use cases, users may instead want to implement their own retry logic or error handling strategy.

```{code-cell} ipython3
def process_cutout_with_error_handling(row, coord, cache):
    '''
    Call `process_cutout` while catching transient read errors.

    Parameters:
    ===========

    row : astropy.table row
        Row of a table that will be changed in place by this function. The table
        is created by the SQL TAP query.
    coord : astropy.coordinates.SkyCoord
        Coordinates (same as used for the TAP query).
    cache : bool
        If set to `True`, the output of cached and the cutout processing will run faster next time.
        Turn this feature off by setting ``cache=False``.
    '''
    try:
        process_cutout(row, coord, cache=cache)
    except (TimeoutError, urllib.error.HTTPError, http.client.IncompleteRead):
        # Transient read errors. Skip this cutout.
        row["hdus"] = None
```

## 7. Download the Cutouts

This process can take a while.
If run in series, it can take about 5 minutes for 700 images on a typical laptop machine.
Here, we therefore exploit two different methods.
First we show the serial approach and next we show how to parallelize the methods.
The latter can be run on many CPUs and is therefore significantly faster.

### 7.1 Serial Approach

First, we implement the serial approach -- a simple `for` loop.
Before that, we turn the results into an astropy table and add some place holders that will be filled in by the `process_cutout()` function.

```{warning}
Running the cell below may take a while for a large number of cutouts.
Approximately 5-7 minutes for 700 images of cutout size 0.01 degree on a typical machine.
```

```{tip}
The astropy `fits.open()` supports a caching argument.
This can be passed through in the `process_cutout()` function.
If `cache=True` is set, the images are cached and the cutout creation is sped up next time the code is run (even if the Jupyter kernel is restarted!).
The downside is that the images are saved on the machine where this notebook is run (usually in `~/.astropy/cache/`).
If many cutouts are created, this can sum up to a large cached data volume, in which case `cache=False` is preferred.

To learn more about the cache please read the [astropy cache management documentation](https://docs.astropy.org/en/stable/utils/data.html#cache-management).
```

```{code-cell} ipython3
results_table_serial = results.to_table()
results_table_serial["cutout_index"] = range(1, len(results_table_serial) + 1)
results_table_serial["central_wavelength"] = np.full(len(results_table_serial), np.nan)
results_table_serial["primary_header"] = np.full(len(results_table_serial), None)
results_table_serial["version"] = np.full(len(results_table_serial), "", dtype="U32")
results_table_serial["hdus"] = np.full(len(results_table_serial), None)

t1 = time.time()
for row in results_table_serial:
    process_cutout_with_error_handling(row, coord, cache=False)
print("Time to create cutouts in serial mode: {:2.2f} minutes.".format((time.time() - t1) / 60))

# Drop rows that failed to download.
failed_cutouts_serial = results_table_serial[results_table_serial["hdus"] == None]
print("Failed to get the following cutouts:\n", failed_cutouts_serial["uri"])
results_table_serial = results_table_serial[results_table_serial["hdus"] != None]
```

### 7.2 Parallel Approach

Next, we implement parallel processing, which will make the cutout creation faster.
The maximal number of workers can be limited by setting the `max_workers` argument.
The choice of this value depends on the number of cores but also on the number of parallel calls that can be digested by the IRSA server.

```{tip}
A good value for the maximum number of workers is between 7 and 12 for a machine with 8 cores.
See the caching tip in section 7.1 above -- the same `cache` argument applies here.
```

Again, before running the cutout processing we define some place holders.

```{code-cell} ipython3
results_table_parallel = results.to_table()
results_table_parallel["cutout_index"] = range(1, len(results_table_parallel) + 1)
results_table_parallel["central_wavelength"] = np.full(len(results_table_parallel), np.nan)
results_table_parallel["primary_header"] = np.full(len(results_table_parallel), None)
results_table_parallel["version"] = np.full(len(results_table_parallel), "", dtype="U32")
results_table_parallel["hdus"] = np.full(len(results_table_parallel), None)

t1 = time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
    futures = [executor.submit(process_cutout_with_error_handling, row, coord, False)
               for row in results_table_parallel]
    concurrent.futures.wait(futures)
print("Time to create cutouts in parallel mode: {:2.2f} minutes.".format((time.time() - t1) / 60))

# Drop rows that failed to download.
failed_cutouts_parallel = results_table_parallel[results_table_parallel["hdus"] == None]
print("Failed to get the following cutouts:\n", failed_cutouts_parallel["uri"])
results_table_parallel = results_table_parallel[results_table_parallel["hdus"] != None]
```

## 8. Create a summary table HDU with renamed columns

In the following, we continue to use the output of the parallel mode.
The following cell does the following:

- Create a summary FITS table.
- Create the final FITS HDU including the summary table.

```{code-cell} ipython3
# Create a summary table HDU with renamed columns
cols = fits.ColDefs([
    fits.Column(name="cutout_index", format="J", array=results_table_parallel["cutout_index"], unit=""),
    fits.Column(name="observation_date", format="D", array=results_table_parallel["time_bounds_lower"], unit="d"),
    fits.Column(name="central_wavelength", format="D", array=results_table_parallel["central_wavelength"], unit="um"),
    fits.Column(name="version", format="A32", array=results_table_parallel["version"], unit=""),
    fits.Column(name="access_url", format="A200", array=results_table_parallel["uri"], unit=""),
])
table_hdu = fits.BinTableHDU.from_columns(cols)
table_hdu.header["EXTNAME"] = "CUTOUT_INFO"
```

## 9. Create the final MEF

Now, we create a primary HDU and combine it with the summary table HDU and the cutout image and wavelength HDUs to a final HDU list.
A blank `fits.PrimaryHDU()` would silently drop IBE's per-file metadata (e.g. `VERSION`), so instead we base the primary header on the first cutout's -- the `version` column in `CUTOUT_INFO` above records each individual cutout's processing version in case they differ.

```{code-cell} ipython3
primary_hdu = fits.PrimaryHDU(header=results_table_parallel["primary_header"][0])
primary_hdu.header["HISTORY"] = (
    "Primary header copied from the first cutout in this MEF. "
    "See the VERSION column in the CUTOUT_INFO table for each cutout's processing version."
)
hdulist_list = [primary_hdu, table_hdu]
hdulist_list.extend(hdu for fits_hdulist in results_table_parallel["hdus"] for hdu in fits_hdulist)
combined_hdulist = fits.HDUList(hdulist_list)
```

## 10. Write the final MEF

Finally, we save the full HDU list to a multi-extension FITS file to disk.
We use `overwrite=False` so that re-running this cell does not silently clobber an existing output file; delete the file first (or change `output_filename` above) if you want to regenerate it.

```{code-cell} ipython3
# Write the final MEF
combined_hdulist.writeto(output_filename, overwrite=False)
```

## 11. Test and visualize the final result

We can now open the new MEF FITS file that we have created above.

Loading the summary table (including wavelength information of each cutout) is straight forward:

```{code-cell} ipython3
summary_table = Table.read(output_filename , hdu=1)
```

Let's also extract the first 10 images (or all of the extracted cutouts if less than 10).

```{code-cell} ipython3
nbr_images = np.nanmin([10, len(summary_table)])
imgs = []
with fits.open(output_filename) as hdul:
    for ii in range(nbr_images):
        extname = "IMAGE{}".format(summary_table["cutout_index"][ii])
        imgs.append(hdul[extname].data)
```

Plot the images with the wavelength corresponding to their central pixel.
We read `central_wavelength` as an astropy `Quantity` (via the column's `.quantity` property, which picks up the `um` unit written into the FITS header) rather than a bare float, so the unit travels with the value instead of being assumed.

```{code-cell} ipython3
central_wavelength = summary_table["central_wavelength"].quantity

fig = plt.figure(figsize=(15, 5))
axs = [fig.add_subplot(2, 5, ii + 1) for ii in range(10)]

for ii, img in enumerate(imgs):
    axs[ii].imshow(imgs[ii], norm="log", origin="lower")
    wl = central_wavelength[ii].to(u.micrometer)
    axs[ii].text(0.05, 0.05, rf"$\lambda_{{\rm center}} = {wl.value:2.4g} \,{{\rm \mu m}}$",
                 va="bottom", ha="left", color="white", transform=axs[ii].transAxes)

plt.show()
```

## Acknowledgements

- [Caltech/IPAC-IRSA](https://irsa.ipac.caltech.edu/)

## About this notebook

**Updated:** 5 March 2026

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.

**Runtime:** As of the date above, this notebook takes about 3 minutes to run to completion on a machine with 8GB RAM and 4 CPU.
