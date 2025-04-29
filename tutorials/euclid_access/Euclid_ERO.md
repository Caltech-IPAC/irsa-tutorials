---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Exploring Star Clusters in the Euclid ERO Data


## Learning Goals
By the end of this tutorial, you will be able to:

 &bull; Access the Euclid ERO images using `astroquery`

 &bull; Create cutouts on the large Euclid ERO images directly

 &bull; Extract sources on the Euclid image and run photometry tools

 &bull; Compare photometry to the Gaia catalog

 &bull; Visualize the Euclid ERO image and the Gaia catalog in `Firefly`


## Introduction
Euclid is a European Space Agency (ESA) space mission with NASA participation, to study the geometry and nature of the dark Universe. As part of its first observations, Euclid publicly released so-called *Early Release Observations* (ERO) targeting some press-worthy targets on the sky such as star clusters or local galaxies.

In this notebook, we will focus on the ERO data of **NGC 6397**, a globular cluster 78,000 light years away. (Note that there is another globular cluster, **NGC 6254** that can also be used for this analysis - in fact the user can choose which one to use) The goal of this analysis is to extract the Euclid photometry of the stars belonging to the cluster and compare them to the photometry of Gaia. Due to the different pixel scales of the visible (VIS, $0.1^{\prime\prime}$) and near-IR (NISP, $0.3^{\prime\prime}$) wavelengths, we will first detect and extract the stars in the VIS filter and use their position to subsequently extract the photometry in the NISP (Y, J, H) filters (a method also referred to as *forced photometry*).

One challenge with Euclid data is their size due to the large sky coverage and small pixel size.
This notebook demonstrates how to download only a cutout of the large Euclid ERO observation image (namely focussing only on the position of the globular cluster. Furthermore, this notebook demonstrates how to extract sources on a large astronomical image and how to measure their photometry across different pixel scales using a very simple implementation of the forced photometry method with positional priors.
Finally, we also demonstrate how to visualize the Euclid images and catalog in `Firefly`, an open-source web-based UI library for astronomical data archive access and visualization developed at Caltech (https://github.com/Caltech-IPAC/firefly). `Firefly` is a convenient tool to visualize images similar to DS9 on your local computer, but it can run on a cloud-based science platform.


This notebook is written to be used in Fornax which is a cloud based computing platform using AWS. The advantage of this is that the user does not need to download actuall data, hence the analysis of large datasets is not limited by local computing resources. It also allows to access data across all archives fast and easy.


## Data Volume
The total data volume required for running this notebook is less than 20 MB.

+++

## Imports

```{important}
We rely on astroquery, firefly_client, photutils, and sep features that have been recently added, so please make sure you have the respective versions v0.4.10, v3.2, v2.0, and v1.4 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install tqdm numpy matplotlib astropy 'photutils>=2.0' 'astroquery>=0.4.10' 'sep>=1.4' 'firefly_client>=3.2'
```

First, we import all necessary packages.


```{code-cell} ipython3
# General imports
import os
import numpy as np
from tqdm import tqdm

# Astroquery
from astroquery.ipac.irsa import Irsa
from astroquery.gaia import Gaia

# Astropy
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

# Photometry tools
import sep
from photutils.detection import DAOStarFinder
from photutils.psf import PSFPhotometry, IterativePSFPhotometry, CircularGaussianSigmaPRF, make_psf_model_image
from photutils.background import LocalBackground, MMMBackground

# Firefly
from firefly_client import FireflyClient

# Matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
```

Next, we define some parameters for `Matplotlib` plotting.

```{code-cell} ipython3
## Plotting stuff
mpl.rcParams['font.size'] = 14
mpl.rcParams['axes.labelpad'] = 7
mpl.rcParams['xtick.major.pad'] = 7
mpl.rcParams['ytick.major.pad'] = 7
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['xtick.minor.top'] = True
mpl.rcParams['xtick.minor.bottom'] = True
mpl.rcParams['ytick.minor.left'] = True
mpl.rcParams['ytick.minor.right'] = True
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['hatch.linewidth'] = 1

def_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
```

## Setting up the Environment

Next, we set up the environment. This includes
* setting up an output data directory (will be created if it does not exist)
* define the search radius around the target of interest to pull data from the Gaia catalog
* define the cutout size that will be used to download a certain part of the Euclid ERO images

Finally, we also define the target of interest here. We can choose between two globular clusters, **NGC 6254** and **NGC6397**, both covered by the Euclid ERO data.

Note that `astropy` units can be attached to the `search_radius` and `cutout_size`.

```{code-cell} ipython3
# create output directory
if os.path.exists("./data/"):
    print("Output directory already created.")
else:
    print("Creating data directory.")
    os.mkdir("./data/")

search_radius = 1.5 * u.arcmin # search radius
cutout_size = 1.5 * u.arcmin # cutout size

coord = SkyCoord.from_name('NGC 6397')
```

## Search Euclid ERO Images

Now, we search for the Euclid ERO images using the `astroquery` package.
Note that the Euclid ERO images are no in the cloud currently, but we access them directly from IRSA using IRSA's *Simple Image Access* (SIA) methods.

```{note}
The following only works for combined images (either extended or point source stacks). This would not work if there are multiple, let's say, H-band images of Euclid at a given position. Therefore, no time domain studies here (which is anyway not one of the main goals of Euclid).
```

The IRSA SIA products can be listed via
```
Irsa.list_collections(servicetype='SIA')
```

Here we use the collection *euclid_ero*, containing the Euclid ERO images. We first create a `SkyCoord` object and then query the SIA.

```{code-cell} ipython3
image_tab = Irsa.query_sia(pos=(coord, search_radius), collection='euclid_ero')
print("Number of images available: {}".format(len(image_tab)))
```

Sort the queried table by wavelength (column `em_min`). This allows us later when we plot the images to keep them sorted by wavelength (VIS, Y, J, H).

```{code-cell} ipython3
image_tab.sort('em_min')
```

Let's inspec the table that was downloaded.

```{code-cell} ipython3
image_tab[0:3]
```

Next, we add additional restrictions to narrow down the images.

The Euclid ERO images come in two different flavors:
* *Flattened* images: idealized for compact sources (for example stars)
* *LSB* images: idealized for low surface brightness objects (for example galaxies)

Maybe counter-intuitively, we select the *LSB* images here by checking if the file name given in the URL (`access_url` column) contains that key word. We found that the *Flattened* images contain many masked pixels.

```{code-cell} ipython3
#sel_basic = np.where( ["Flattened" in tt["access_url"] for tt in image_tab] )[0]
sel_basic = np.where( ["LSB" in tt["access_url"] for tt in image_tab] )[0]
image_tab = image_tab[sel_basic]
```

Next, we want to check what filteres are available. This can be done by `np.unique()`, however, in that case it would sort the filters alphabetically. We want to keep the sorting based on the wavelength, therefore we opt for a more complicated way.

```{code-cell} ipython3
idxs = np.unique(image_tab["energy_bandpassname"], return_index=True)[1]
filters = [image_tab["energy_bandpassname"][idx] for idx in sorted(idxs)]
print("Filters: {}".format(filters))
```

We can now loop throught the filters and collect the images to create a handy summary table with all the data we have. This way we will have an easier time to later access the data.

```{code-cell} ipython3
# Create a dictionary with all the necessary information (science, weight, noise, mask)
summary_table = Table(names=["filter","products","facility_name","instrument_name"] , dtype=[str,str,str,str])
for filt in filters:
    sel = np.where(image_tab['energy_bandpassname'] == filt)[0]
    products = list( np.unique(image_tab["dataproduct_subtype"][sel].value) )
    if "science" in products: # sort so that science is the first in the list. This is the order we create the hdu extensions
        products.remove("science")
        products.insert(0,"science")
    else:
        print("WARNING: no science image found!")
    print(products)

    summary_table.add_row( [filt , ";".join(products), str(np.unique(image_tab["facility_name"][sel].value)[0]), str(np.unique(image_tab["instrument_name"][sel].value)[0])] )
```

Let's check out the summary table that we have created. We see that we have all the 4 Euclid bands and what data products are available for each of them.

```{code-cell} ipython3
summary_table
```

## Create Cutout Images

Now that we have a list of data products, we can create the cutouts. This is important as the full Euclid ERO images would be too large to run extraction and photometry software on them (they would simply fail due to memory issues).

For each image, we create a cutout around the target of interest, using the `cutout_size` defined above. The cutouts are then collected into HDUs. That way we can easily access the different data products. Note that we only use the *science* products as the *ancillary* products are not needed.

We save the HDU to disk as it will be later used when we visualize the Euclid ERO FITS images in `Firefly`.

```{note}
You will notice that `Cutout2D` can be applied to an URL. That way, it we do not need to download the full image to create a cutout. This is a useful trick to keep in mind when analyzing large images. This makes creating cutout images very fast.
```

```{code-cell} ipython3
%%time
for ii,filt in tqdm(enumerate(filters)):
    products = summary_table[summary_table["filter"] == filt]["products"][0].split(";")

    for product in products:
        sel = np.where( (image_tab["energy_bandpassname"] == filt) & (image_tab["dataproduct_subtype"] == product) )[0]

        with fits.open(image_tab['access_url'][sel[0]], use_fsspec=True) as hdul:
            tmp = Cutout2D(hdul[0].section, position=coord, size=cutout_size, wcs=WCS(hdul[0].header)) # create cutout


            if (product == "science") & (ii == 0): # if science image, then create a new hdu.
                hdu0 = fits.PrimaryHDU(data = tmp.data, header=hdul[0].header)
                hdu0.header.update(tmp.wcs.to_header()) # update header with WCS info
                hdu0.header["EXTNAME"] = "{}_{}".format(filt,product.upper())
                hdu0.header["PRODUCT"] = product.upper()
                hdu0.header["FILTER"] = filt.upper()
                hdulcutout = fits.HDUList([hdu0])
            elif (product == "science") & (ii > 0):
                hdu = fits.ImageHDU(data = tmp.data, header=hdul[0].header)
                hdu.header.update(tmp.wcs.to_header()) # update header with WCS info
                hdu.header["EXTNAME"] = "{}_{}".format(filt,product.upper())
                hdu.header["PRODUCT"] = product.upper()
                hdu.header["FILTER"] = filt.upper()
                hdulcutout.append(hdu)

## Save the HDUL cube:
hdulcutout.writeto("./data/euclid_images_test.fits", overwrite=True)
```

Let's look at the HDU that we created to check what information we have. You see that all filters are collected in different extensions. Also note the different dimensions of the FITS layers as the pixel scale of VIS and NISP are different.

```{code-cell} ipython3
hdulcutout.info()
```

We can now plot the image cutouts that we generated. The globular cluster is clearly visible.

```{code-cell} ipython3
nimages = len(filters) # number of images
ncols = int(4) # number of columns
nrows = int( nimages // ncols ) # number of rows

fig = plt.figure(figsize = (5*ncols,5*nrows) )
axs = [fig.add_subplot(nrows,ncols,ii+1) for ii in range(nimages)]

for ii,filt in enumerate(filters):
    img = hdulcutout["{}_SCIENCE".format(filt)].data
    axs[ii].imshow(img , origin="lower")
    axs[ii].text(0.05,0.05 , "{} ({}/{})".format(summary_table["facility_name"][ii],summary_table["instrument_name"][ii],filt) , fontsize=14 , color="white",
                 va="bottom", ha="left" , transform=axs[ii].transAxes)

plt.show()
```

## Extract Sources and Measure their Photometry on the VIS Image

Now that we have the images in memory (and on disk - but we do not need them, yet), we can measure the fluxes of the individual stars.
Our simple photometry pipeline has different parts:

* First, we are using the Python package `sep` (similar to SExtractor) to extract the position of the sources. We do that on the VIS image, which provides the highest resolution.
* Second, we use `sep` to perform aperture measurements of the photometry. We will use that to compare the obtained fluxes to our forced photometry method
* Third, we apply a PSF fitting technique (using the `photutils` Python package) to improve the photometry measurement

+++

We start by extracting the sources using `sep`. We first isolate the data that we want to look at (the VIS image only).

```{code-cell} ipython3
## Get Data (this will be replaced later)
img = hdulcutout["VIS_SCIENCE"].data
hdr = hdulcutout["VIS_SCIENCE"].header
img[img == 0] = np.nan
```

There are some NaN value on the image that we need to mask out. For this we create a mask image that we later feed to `sep`.

```{code-cell} ipython3
mask = np.isnan(img)
```

Next, we compute the background statistics what will be used by `sep` to extract the sources above a certain threshold.

```{code-cell} ipython3
mean, median, std = sigma_clipped_stats(img, sigma=3.0)
print(np.array((mean, median, std)))
```

Finally, we perform object detection with `sep`. There are also modules in `photutils` to do that, however, we found that `sep` works best here. We output the number of objects found on the image.

```{code-cell} ipython3
objects = sep.extract(img-median, thresh=1.2, err=std, minarea=3, mask=mask, deblend_cont=0.0002, deblend_nthresh=64 )
print("Number of sources extracted: ", len(objects))
```

Next, we perform simple aperture photometry on the detected sources. Again, we use the `sep` package for this. We will use these aperture photometry later to compare to the PSF photometry.

```{code-cell} ipython3
flux, fluxerr, flag = sep.sum_circle(img-median, objects['x'], objects['y'],r=3.0, err=std, gain=1.0)
```

Now, we use the `photutils` Python package to perform PSF fitting. Here we assume a simple Gaussian with a FWHM given by `psf_fwhm` as PSF.

```{note}
We use a Gaussian PSF here for simplicity. The photometry can be improved by using a pixelated PSF measured directly from the Euclid images (for example by stacking stars).
```

```{code-cell} ipython3
psf_fwhm = 0.16 # PSF FWHM in arcsec
pixscale = 0.1 # VIS pixel scale in arcsec/px

init_params = Table([objects["x"],objects["y"]] , names=["x","y"]) # initial positions
psf_model = CircularGaussianSigmaPRF(flux=1, sigma=psf_fwhm/pixscale / 2.35)
psf_model.x_0.fixed = True
psf_model.y_0.fixed = True
psf_model.sigma.fixed = False
fit_shape = (5, 5)
psfphot = PSFPhotometry(psf_model,
                        fit_shape,
                        finder = DAOStarFinder(fwhm=0.1, threshold=3.*std, exclude_border=True), # not needed because we are using fixed initial positions.
                        aperture_radius = 4,
                        progress_bar = True)
```

After initiating the `PSFPhotometry` object, we can run the PSF photometry measurement. This can take a while (typically between 1 and 2 minutes).

```{code-cell} ipython3
phot = psfphot(img-median, error=None, mask=mask, init_params=init_params)
```

Once this is done, we can create a residual image.

```{code-cell} ipython3
resimage = psfphot.make_residual_image(data = img-median, psf_shape = (9, 9))
```

We now want to add the best-fit coordinates (R.A. and Decl.) to the VIS photometry catalog. For this, we have to convert the image coordinates into sky coordinates using the WCS information. We will need these coordinates because we want to use them as positional priors for the photometry measurement on the NISP images.

```{code-cell} ipython3
## Add coordinates to catalog
wcs1 = WCS(hdr) # VIS
radec = wcs1.all_pix2world(phot["x_fit"],phot["y_fit"],0)
phot["ra_fit"] = radec[0]
phot["dec_fit"] = radec[1]
```

Finally, we plot the VIS image and the residual with the extracted sources overlaid.

```{code-cell} ipython3
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax1.imshow(np.log10(img), cmap="Greys", origin="lower")
ax1.plot(phot["x_fit"], phot["y_fit"] , "o", markersize=8 , markeredgecolor="red", fillstyle="none")

ax2.imshow(resimage,vmin=-5*std, vmax=5*std, cmap="RdBu", origin="lower")
ax2.plot(phot["x_fit"], phot["y_fit"] , "o", markersize=8 , markeredgecolor="red", fillstyle="none")

plt.show()
```

As an additional check, we can compare the aperture photometry and the PSF photometry. We should find a good agreement between those two measurement methods. However, note that the PSF photometry should do a better job in deblending the fluxes of sources that are close by.

```{code-cell} ipython3
x = objects["flux"]
y = phot["flux_fit"]

fig = plt.figure(figsize=(6,5))
ax1 = fig.add_subplot(1,1,1)

ax1.plot(x , y , "o", markersize=2, alpha=0.5)
minlim = np.nanmin(np.concatenate((x,y)))
maxlim = np.nanmax(np.concatenate((x,y)))

ax1.fill_between(np.asarray([minlim,maxlim]),np.asarray([minlim,maxlim])/1.5,np.asarray([minlim,maxlim])*1.5, color="gray", alpha=0.2, linewidth=0)
ax1.fill_between(np.asarray([minlim,maxlim]),np.asarray([minlim,maxlim])/1.2,np.asarray([minlim,maxlim])*1.2, color="gray", alpha=0.4, linewidth=0)
ax1.plot(np.asarray([minlim,maxlim]),np.asarray([minlim,maxlim]), ":", color="gray")

ax1.set_xlabel("Aperture Photometry [flux]")
ax1.set_ylabel("PSF forced-photometry [flux]")
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.show()
```

## Measure the Photometry on the NISP Images

We now have the photometry and the position of sources on the VIS image. We can now proceed with similar steps on the NISP images. Because the NISP PSF and pixel scale are larger that those of the VIS images, we utilize the advantage of position prior-based forced photometry.
For this, we use the positions of the VIS measurements and perform PSF fitting on the NISP image using these priors.

The steps below are almost identical to the ones applied to the VIS images.

+++

We first isolate the data, which is in this case the NISP *H*-band filter. Note that this is an arbitrary choice and you should be encouraged to try other filters as well!

```{code-cell} ipython3
img2 = hdulcutout["H_SCIENCE"].data
hdr2 = hdulcutout["H_SCIENCE"].header
img2[img2 == 0] = np.nan
```

We again need to create a mask that will be fed to the PSF fitting module.

```{code-cell} ipython3
mask2 = np.isnan(img2)
```

... and we also get some statistics on the sky background.

```{code-cell} ipython3
mean2, median2, std2 = sigma_clipped_stats(img2, sigma=3.0)
print(np.array((mean2, median2, std2)))
```

Now, we need to obtain the (x,y) image coordinates on the NISP image that correspond to the extracted sources on the VIS image. We use the WCS information from the NISP image for this case and apply it to the sky coordinates obtained on the VIS image.

```{code-cell} ipython3
wcs = WCS(hdr) # VIS
wcs2 = WCS(hdr2) # NISP
radec = wcs.all_pix2world(objects["x"],objects["y"], 0)
xy = wcs2.all_world2pix(radec[0],radec[1],0)
```

Having all this set up, we can now again perform the PSF photometry using `PSFPhotometry()` from the `photutils` package. This again can take a while, typically 1 minute.

```{code-cell} ipython3
psf_fwhm = 0.3 # arcsec
pixscale = 0.3 # arcsec/px

init_params = Table([xy[0],xy[1]] , names=["x","y"]) # initial positions
psf_model = CircularGaussianSigmaPRF(flux=1, sigma=psf_fwhm/pixscale / 2.35)
psf_model.x_0.fixed = True
psf_model.y_0.fixed = True
psf_model.sigma.fixed = False
fit_shape = (3, 3)
psfphot2 = PSFPhotometry(psf_model,
                        fit_shape,
                        finder = DAOStarFinder(fwhm=0.1, threshold=3.*std2, exclude_border=True), # not needed because we are using fixed initial positions.
                        aperture_radius = 4,
                        progress_bar = True)
phot2 = psfphot2(img2-median2, error=None, mask=mask2, init_params=init_params)
resimage2 = psfphot2.make_residual_image(data = img2-median2, psf_shape = (3, 3))
```

Again, we convert the pixel coordinates to sky coordinates and add them to the catalog.

```{code-cell} ipython3
wcs2 = WCS(hdr2) # NISP
radec = wcs2.all_pix2world(phot2["x_fit"],phot2["y_fit"],0)
phot2["ra_fit"] = radec[0]
phot2["dec_fit"] = radec[1]
```

Finally, we create the same figure as above, showing the NISP image and the residual with the (VIS-extracted) sources overlaid.

```{code-cell} ipython3
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax1.imshow(np.log10(img2), cmap="Greys", origin="lower")
ax1.plot(phot2["x_fit"], phot2["y_fit"] , "o", markersize=8 , markeredgecolor="red", fillstyle="none")

ax2.imshow(resimage2,vmin=-20*std2, vmax=20*std2, cmap="RdBu", origin="lower")
ax2.plot(phot2["x_fit"], phot2["y_fit"] , "o", markersize=8 , markeredgecolor="red", fillstyle="none")

plt.show()
```

## Load Gaia Catalog

We now load the Gaia sources at the location of the globular clusters. The goal is to compare the photometry of Gaia to the one derived above for the Euclid VIS and NISP images. This is scientifically useful, for example we can compute the colors of the stars in the Gaia optical bands and the Euclid near-IR bands.
To search for Gaia sources, we use `astroquery` again.

We first have to elimiate the row limit for the Gaia query by setting

```{code-cell} ipython3
Gaia.ROW_LIMIT = -1
```

Next, we request the Gaia catalog around the position of the globular cluster. We use the same size as the cutout size.

```{code-cell} ipython3
gaia_objects = Gaia.query_object_async(coordinate=coord, radius = cutout_size/2)
print("Number of Gaia stars found: {}".format(len(gaia_objects)))
```

We then convert the sky coordinates of the Gaia stars to (x,y) image coordinates for VIS and NISP images using the corresponding WCS. This makes it more easy to plot the Gaia sources later on the images.

```{code-cell} ipython3
wcs = WCS(hdr) # VIS
wcs2 = WCS(hdr2) # NISP
xy = wcs.all_world2pix(gaia_objects["ra"],gaia_objects["dec"],0)
xy2 = wcs2.all_world2pix(gaia_objects["ra"],gaia_objects["dec"],0)

gaia_objects["x_vis"] = xy[0]
gaia_objects["y_vis"] = xy[1]
gaia_objects["x_nisp"] = xy2[0]
gaia_objects["y_nisp"] = xy2[1]
```

We save the Gaia table to disk as we will later use it for the visualization in `Firefly`.

```{code-cell} ipython3
gaia_objects.write("./data/gaiatable.csv", format="csv", overwrite=True)
```

Now we can overlay the Gaia sources on the VIS and NISP images (here the x/y coordinates become handy).

```{code-cell} ipython3
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax1.imshow(np.log10(img), cmap="Greys", origin="lower")
ax1.plot(gaia_objects["x_vis"], gaia_objects["y_vis"] , "o", markersize=8 , markeredgecolor="red", fillstyle="none")
ax1.set_title("VIS")

ax2.imshow(np.log10(img2), cmap="Greys", origin="lower")
ax2.plot(gaia_objects["x_nisp"], gaia_objects["y_nisp"] , "o", markersize=8 , markeredgecolor="red", fillstyle="none")
ax2.set_title("NISP")

plt.show()
```

## Match the Gaia Catalog to the VIS and NISP Catalogs

Now, we match the Gaia source positions to the extracted sources in the VIS and NISP images.

We first define which Gaia columns to copy to the matched catalog as well as the matching distance.

```{code-cell} ipython3
gaia_keys = ["source_id", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag","ra","dec","pmra","pmdec"]
matching_distance = 0.6*u.arcsecond
```

First match to the VIS image. We use the `astropy` *SkyCoord()* function for matching in sky coordinates.

```{code-cell} ipython3
c = SkyCoord(ra=phot["ra_fit"]*u.degree, dec=phot["dec_fit"]*u.degree )
catalog = SkyCoord(ra=gaia_objects["ra"].data*u.degree, dec=gaia_objects["dec"].data*u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)

sel_matched = np.where(d2d.to(u.arcsecond) < (matching_distance))[0]
print("Gaia Sources matched to VIS: {}".format( len(sel_matched) ) )
phot["gaia_distance"] = d2d.to(u.arcsecond)

for gaia_key in gaia_keys:
    phot["gaia_{}".format(gaia_key)] = 0.0
    phot["gaia_{}".format(gaia_key)][sel_matched] = gaia_objects[gaia_key][idx[sel_matched]]
```

And then we add the NISP sources. Note that we do not have to perform matching here because by design the VIS and NISP sources are matched (spatial prior forced photometry).

```{code-cell} ipython3
phot2["gaia_distance"] = d2d.to(u.arcsecond)

for gaia_key in gaia_keys:
    phot2["gaia_{}".format(gaia_key)] = 0.0
    phot2["gaia_{}".format(gaia_key)][sel_matched] = gaia_objects[gaia_key][idx[sel_matched]]
```

Once matched, we can now compare the Gaia and Euclid/NISP magnitudes of the stars.

```{code-cell} ipython3
# Data
x = phot["gaia_phot_rp_mean_mag"]
y = -2.5*np.log10(phot["flux_fit"]) + hdr["ZP_STACK"]

# selection
sel_good = np.where(phot["gaia_source_id"] > 0)[0]
x = x[sel_good]
y = y[sel_good]

fig = plt.figure(figsize=(6,5))
ax1 = fig.add_subplot(1,1,1)

ax1.plot(x , y , "o", markersize=2)
minlim = np.nanmin(np.concatenate((x,y)))
maxlim = np.nanmax(np.concatenate((x,y)))

ax1.fill_between(np.asarray([minlim,maxlim]),np.asarray([minlim,maxlim])/1.5,np.asarray([minlim,maxlim])*1.5, color="gray", alpha=0.2, linewidth=0)
ax1.fill_between(np.asarray([minlim,maxlim]),np.asarray([minlim,maxlim])/1.2,np.asarray([minlim,maxlim])*1.2, color="gray", alpha=0.4, linewidth=0)
ax1.plot(np.asarray([minlim,maxlim]),np.asarray([minlim,maxlim]), ":", color="gray")

ax1.set_xlabel("Gaia Rp [mag]")
ax1.set_ylabel("I$_E$ [mag]")
plt.show()
```

## Visualization with Firefly

At the end of this Notebook, we demonstrate how we can visualize the images and catalogs created above in `Firefly`.

We start by initializing the Firefly client.
The following line will open a new `Firefly` GUI in a separate tab **inside** the Jupyter Notebook environment. The user can drag the tab onto the currently open tab to create a "split tab". This the user to see the code and images side-by-side.

```{code-cell} ipython3
# Uncomment when using within Jupyter Lab with jupyter_firefly_extensions installed
# fc = FireflyClient.make_lab_client()

# Uncomment for contexts other than the above
fc = FireflyClient.make_client(url="https://irsa.ipac.caltech.edu/irsaviewer")
```

In order to display in image or catalog in `Firefly`, it needs to be uploaded to the `Firefly` server. We do this here using the `upload_file()` function.
We first upload the FITS image that we created above.

```{code-cell} ipython3
fval = fc.upload_file('./data/euclid_images_test.fits')
```

Once the image is uploaded we can use the `show_fits()` function to display it.
Note that our FITS image has multiple extensions (VIS, and NISP bands). We can open them separately in new `Firefly` tabs by looping over the HDUs and specifying the plot ID by the extension's name.

```{code-cell} ipython3
for hh,hdu in enumerate(hdulcutout):
    fc.show_fits(fval, MultiImageIdx=hh, plot_id=hdu.header["EXTNAME"] )
```

We can lock the WCS between the images (allowing the user to pan and zoom the images simultaneously) by running:

```{code-cell} ipython3
fc.align_images(lock_match=True)
```

In the same way, we can upload a table, in this case our Gaia table. We again use `upload_file()` but in this case we use `show_table()` to show it in `Firefly`.

```{code-cell} ipython3
tval = fc.upload_file('./data/gaiatable.csv')
fc.show_table(tval, tbl_id = "gaiatable")
```

Now, check out the `Firefly` GUI. You can zoom the images, click on sources, filter the table, display different selection, and much more!

+++

***

## About this Notebook

**Author**: Andreas Faisst (IPAC Scientist)

**Updated**: 2025-03-17

**Contact**: the [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
