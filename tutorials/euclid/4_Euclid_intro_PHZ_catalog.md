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
- know how to find and plot the spectra of one of the galaxies from its object ID, and compare its photometric and spectroscopic redshifts;
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

Among the data products included in DR1 are the catalogs created by the PHZ processing function, which estimates the photometric redshift of every source in the MER catalog from its photometry, and classifies it as a star or a galaxy.
In DR1, the photometric redshifts are measured with a nearest-neighbour method (NNPZ), which compares the photometry of each source with that of a reference sample of galaxies with precise redshifts in the COSMOS field.
This notebook provides an introduction to the main PHZ catalog, which contains the photometric redshift probability distribution, the unified fluxes, and the classification for each source.
We follow a typical path through the catalog: we select galaxies with good photometric redshifts in a redshift range, look at them on a MER mosaic, and compare the photometric redshift of one of them with its spectrum and spectroscopic redshift.
The PHZ processing for DR1 is described in detail in Euclid Collaboration: Hartley et al. (2026), and the galaxy physical properties it also provides in Euclid Collaboration: Enia et al. (2026); see the [References](#references) section below.
If you have questions about this notebook, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## The PHZ catalog tables

Like the MER catalog (see the [MER catalogs notebook](2_Euclid_intro_MER_catalog.md)), the DR1 PHZ catalogs are split by survey:

* `euclid_dr1_phz_photo_z_wide` -- the main PHZ catalog of the **Euclid Wide Survey**: photometric redshifts, unified fluxes, and classification.
* `euclid_dr1_phz_photo_z_deep_survey` -- the same for the **Euclid Deep Survey** (the deep tiles of the Euclid Deep Fields North, South, and Fornax).
* `euclid_dr1_phz_classification_wide` and `euclid_dr1_phz_classification_deep_survey` -- the star and galaxy probabilities, and the thresholds above which a source is accepted as a star or a galaxy, behind the `phz_classification` column of the main catalogs.

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
from astropy.table import QTable
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
table_spe = 'euclid_dr1_spectro_z_spe_galaxy_candidates_deep'
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
- `phz_mode_1`, the mean redshift of the most probable peak of the PDF, and `phz_mode_2`, the mean redshift of the PDF of the closest neighbour in the reference sample (despite their names, neither is exactly the mode of the PDF);
- `phz_classification`: 1 for a star, 2 for a galaxy, and 3 when the source is accepted as both. The data model also defines 4 for a QSO and 8 for a globular cluster, but DR1 only classifies sources as stars or galaxies: QSO candidates are the point-like sources that are not accepted as stars;
- `phz_flags`, the PHZ processing flags, 0 when all went well (for example, 1 when no close neighbour was found in the reference sample);
- `flux_vis_unif`, `flux_y_unif`, `flux_j_unif`, `flux_h_unif` and the ground-based bands: the "unified" fluxes, the MER fluxes corrected for Galactic extinction and put on a common filter system (the DECam filters for the ground-based bands), so that they can be compared across the whole survey. They are computed for all sources, but should not be used for stars.
```

```{code-cell} ipython3
# Full list of columns and their description
columns_info
```

## 2. Select galaxies between z = 0.85 and 1.85 with good photometric redshifts

We join the PHZ catalog with the MER catalog, which has the coordinates, and select galaxies in a 5 arcmin box around a position in the Euclid Deep Field North, with the following conditions:
- the unified fluxes in the four Euclid bands are greater than zero;
- `phz_classification = 2`, so the source is accepted as a galaxy only;
- `phz_flags = 0`, so the PHZ processing had no problems;
- the 90% interval of the photometric redshift PDF is narrower than 20% of (1 + z), so the photometric redshift is well constrained;
- the median photometric redshift is between 0.85 and 1.85, the range in which the Hα line falls in the NISP red grism spectra (1.2 to 1.9 micron), so that the galaxies can also have spectroscopic redshifts (see below).

In the same query, we also get the best spectroscopic redshift of each galaxy (`spe_rank = 0`), with its probability `spe_z_prob` and reliability `spe_z_rel`, from the DR1 SPE catalog (see the [SPE catalogs notebook](5_Euclid_intro_SPE_catalog.md)).
We use a left outer join, which keeps the galaxies without a spectroscopic redshift; the condition `spe_rank = 0` goes in the `ON` clause of the join rather than in the `WHERE` clause, which would drop them.

```{code-cell} ipython3
######################## User defined section ############################
# Set the center and size of the region to select galaxies from (and of the image cutout below)
ra_cutout = 269.46
dec_cutout = 66.0
im_cutout = 5 * u.arcmin

# Redshift range
z_min, z_max = 0.85, 1.85
##########################################################################

coords_cutout = SkyCoord(ra_cutout, dec_cutout, unit='deg', frame='icrs')
size_cutout = im_cutout.to(u.deg).value
```

```{code-cell} ipython3
adql = f"""
    SELECT mer.object_id, mer.ra, mer.dec,
           phz.flux_vis_unif, phz.flux_y_unif, phz.flux_j_unif, phz.flux_h_unif,
           phz.phz_classification, phz.phz_median, phz.phz_90_int_1, phz.phz_90_int_2,
           spe.spe_z, spe.spe_z_prob, spe.spe_z_rel
    FROM {table_mer} AS mer
    JOIN {table_phz} AS phz ON mer.object_id = phz.object_id
    LEFT OUTER JOIN {table_spe} AS spe ON spe.object_id = mer.object_id AND spe.spe_rank = 0
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

galaxies = Irsa.query_tap(adql, async_job=True).to_table()
print(f'Found {len(galaxies)} galaxies')
galaxies[:5]
```

## 3. Compare the photometric and spectroscopic redshifts

The Euclid Deep Fields are also covered by NISP slitless spectroscopy, so many of these galaxies have a spectroscopic redshift.
The galaxies without one have a masked `spe_z`:

```{code-cell} ipython3
has_spec_z = ~np.ma.getmaskarray(galaxies['spe_z'])

print(f'{has_spec_z.sum()} of the {len(galaxies)} galaxies have a spectroscopic redshift')
```

Not all of these spectroscopic redshifts can be trusted: most of the galaxies are faint, and when a single emission line is detected, it can be attributed to the wrong line, which gives a completely wrong redshift.
As in the [SPE catalogs notebook](5_Euclid_intro_SPE_catalog.md), we consider a spectroscopic redshift reliable when both its probability and its reliability are high (`spe_z_prob > 0.99` and `spe_z_rel > 0.05`):

```{code-cell} ipython3
reliable = has_spec_z & (galaxies['spe_z_prob'].filled(0) > 0.99) & (galaxies['spe_z_rel'].filled(0) > 0.05)

print(f'{reliable.sum()} of the {has_spec_z.sum()} spectroscopic redshifts are reliable')
```

Let's compare the two redshifts of the galaxies with a reliable spectroscopic redshift:

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(galaxies['spe_z'][reliable], galaxies['phz_median'][reliable], color='tab:blue')
ax.plot([0.5, 2.5], [0.5, 2.5], color='grey', linestyle='--', linewidth=0.8)
ax.axhspan(z_min, z_max, color='tab:red', alpha=0.1, label='Selected photometric redshift range')
ax.set_xlim(0.5, 2.5)
ax.set_ylim(0.5, 2.5)
ax.set_xlabel('Spectroscopic redshift (spe_z)')
ax.set_ylabel('Photometric redshift (phz_median)')
_ = ax.legend(loc='upper left')
```

Most of the galaxies lie close to the one-to-one line.
Let's quantify this with the usual metrics: the normalized median absolute deviation (NMAD) of (z<sub>phot</sub> - z<sub>spec</sub>) / (1 + z<sub>spec</sub>), and the fraction of outliers, for which this difference is larger than 0.15.
We also count how many spectroscopic redshifts fall within the 90% interval of the photometric redshift:

```{code-cell} ipython3
spec_z_reliable = np.asarray(galaxies['spe_z'][reliable], dtype=float)
dz = (np.asarray(galaxies['phz_median'][reliable]) - spec_z_reliable) / (1 + spec_z_reliable)
nmad = 1.48 * np.median(np.abs(dz - np.median(dz)))
inside_90 = ((spec_z_reliable >= galaxies['phz_90_int_1'][reliable])
             & (spec_z_reliable <= galaxies['phz_90_int_2'][reliable]))

print(f'NMAD = {nmad:.3f}, outlier fraction = {np.mean(np.abs(dz) > 0.15):.1%}')
print(f'{inside_90.mean():.0%} of the spectroscopic redshifts are within the 90% interval of the photometric redshift')
```

The scatter and outlier fraction are similar to those found by Euclid Collaboration: Hartley et al. (2026) for the whole Euclid Deep Field North.
However, far fewer than 90% of the spectroscopic redshifts fall within the 90% interval: the photometric redshift PDFs of the deep fields are too narrow in DR1, because no systematic uncertainty was added to the fluxes in their processing (Euclid Collaboration: Hartley et al. 2026).
Keep this in mind when using the PDF intervals, as we did in section 2 to select galaxies with a well-constrained photometric redshift.
Try plotting the unreliable spectroscopic redshifts as well: they are scattered over the whole redshift range, because a single detected line attributed to the wrong transition gives a completely wrong redshift.

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

We keep the VIS mosaic, which we use for the cutouts and in Firefly below:

```{code-cell} ipython3
vis_mosaic = science_images[science_images['energy_bandpassname'] == 'VIS'][0]
filename = vis_mosaic['access_url']
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

We overplot the galaxies with a reliable spectroscopic redshift, colored by their spectroscopic redshift.
In section 5, we look at the spectrum of the brightest of them in the H band, which we mark with a star:

```{code-cell} ipython3
brightest = galaxies[reliable]
if len(brightest) == 0:
    raise ValueError('No galaxy has a reliable spectroscopic redshift; try a larger region or a different redshift range.')
brightest.sort(keys='flux_h_unif', reverse=True)
```

```{code-cell} ipython3
fig = plt.figure(figsize=(8, 7))
ax = plt.subplot(projection=cutout_image.wcs)

ax.imshow(cutout_image.data, cmap='gray', origin='lower',
          norm=ImageNormalize(cutout_image.data, interval=PercentileInterval(99.9), stretch=LogStretch()))
# Open circles, with their edge colored by the spectroscopic redshift
spec_z_values = np.asarray(brightest['spe_z'], dtype=float)
norm = plt.Normalize(spec_z_values.min(), spec_z_values.max())
cmap = plt.get_cmap('cool')
ax.scatter(brightest['ra'], brightest['dec'], s=60, facecolors='none', edgecolors=cmap(norm(spec_z_values)),
           linewidths=1.5, transform=ax.get_transform('world'))
ax.scatter(brightest['ra'][0], brightest['dec'][0], s=300, marker='*', facecolors='none', edgecolors='yellow',
           linewidths=1.5, transform=ax.get_transform('world'))

fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='Spectroscopic redshift (spe_z)')
ax.set_title('Galaxies with a reliable spectroscopic redshift')
plt.show()
```

## 5. Look at the spectrum and image of one of the brightest galaxies

Here are the five brightest galaxies with a reliable spectroscopic redshift, sorted by their H-band flux:

```{code-cell} ipython3
brightest['object_id', 'flux_h_unif', 'phz_median', 'phz_90_int_1', 'phz_90_int_2', 'spe_z', 'spe_z_prob', 'spe_z_rel'][:5]
```

### Retrieve the spectra

The DR1 SIR 1D spectra are stored in large FITS files with one extension (HDU) per spectrum.
The table `euclid.objectid_spectrafile_association_dr1` at IRSA lists, for each object ID, the file (`path`) and the HDU of each of its spectra.
To read a spectrum, we don't download its whole file: IRSA's spectrumDM service reads a single HDU of a spectrum file in the cloud and returns it as a VOTable following the IVOA Spectrum Data Model, which astropy can read directly.
Its URL is built from the file path and the HDU number (see the [SIR 1D spectra notebook](3_Euclid_intro_1D_spectra.md)):

```{code-cell} ipython3
SPECTRUMDM_URL = 'https://irsacloud.ipac.caltech.edu/api/spectrumdm/convert/'
SPECTRUM_DATASET_ID = 'euclid_combspec_dr1'


def spectrum_url(path, hdu):
    """Build the spectrumDM URL of one spectrum.

    Parameters
    ----------
    path : str
        Path of the spectrum file, from the ``path`` column of the association table.
    hdu : int
        Number of the HDU of the spectrum in the file.

    Returns
    -------
    str
        URL that returns the spectrum as a VOTable.
    """
    return f'{SPECTRUMDM_URL}{path}?dataset_id={SPECTRUM_DATASET_ID}&hdu={hdu}'


def find_spectra(object_id):
    """Find the spectra of an object in the spectrum association table.

    Parameters
    ----------
    object_id : int
        MER object ID of the source.

    Returns
    -------
    astropy.table.Table
        One row per spectrum, with the file path and HDU number of the spectrum, and its spectrumDM URL.
    """
    query = f"""
        SELECT objectid, tileid, path, hdu
        FROM euclid.objectid_spectrafile_association_dr1
        WHERE objectid = {object_id}
    """
    spectrum_files = Irsa.query_tap(query).to_table()
    spectrum_files['url'] = [spectrum_url(path, hdu) for path, hdu in zip(spectrum_files['path'], spectrum_files['hdu'])]
    return spectrum_files
```

We find the spectra of the brightest galaxy; if a galaxy has no spectrum listed, we move on to the next brightest one:

```{code-cell} ipython3
for bright_galaxy in brightest:
    result_spectra = find_spectra(bright_galaxy['object_id'])
    if len(result_spectra) > 0:
        break
    print(f"No spectra listed for object {bright_galaxy['object_id']}, trying the next galaxy")
else:
    raise ValueError('None of the galaxies with a reliable spectroscopic redshift has a spectrum listed; '
                     'try a larger region or a different redshift range.')

obj_id = bright_galaxy['object_id']
photo_z = bright_galaxy['phz_median']
spec_z = bright_galaxy['spe_z']
coords_galaxy = SkyCoord(bright_galaxy['ra'], bright_galaxy['dec'], unit='deg')

print(f'Object {obj_id}: photometric redshift {photo_z:.2f}, spectroscopic redshift {spec_z:.3f}')
result_spectra['objectid', 'tileid', 'hdu']
```

Each URL returns a VOTable containing the spectrum, which we read into an astropy table:

```{code-cell} ipython3
spectra = [QTable.read(url) for url in result_spectra['url']]
```

A galaxy in the Euclid Deep Fields can have spectra from both the red and the blue grism; as in the [SIR 1D spectra notebook](3_Euclid_intro_1D_spectra.md), we use the MASK column to find the flux bins not to use, and keep the spectrum with the most usable flux bins for each kind of grism:

```{code-cell} ipython3
def bad_flux_bins(spectrum):
    """Find the flux bins of a spectrum that should not be used.

    Parameters
    ----------
    spectrum : astropy.table.QTable
        Spectrum, with a MASK column.

    Returns
    -------
    numpy.ndarray
        Boolean array, True for the flux bins with an odd MASK value or a MASK value of 64 or more.
    """
    mask = spectrum['MASK'].value
    return (mask % 2 == 1) | (mask >= 64)


def grism_name(spectrum):
    """Tell which kind of grism a spectrum comes from.

    Parameters
    ----------
    spectrum : astropy.table.QTable
        Spectrum, with a WAVELENGTH column.

    Returns
    -------
    str
        ``'Blue grism'`` if the spectrum starts below 1.1 micron, ``'Red grism'`` otherwise.
    """
    return 'Blue grism' if spectrum['WAVELENGTH'].min() < 1.1 * u.micron else 'Red grism'


def best_spectra(spectra):
    """Pick the spectrum with the most usable flux bins for each kind of grism.

    Parameters
    ----------
    spectra : list of astropy.table.QTable
        Spectra of one object.

    Returns
    -------
    dict
        Best spectrum of each kind of grism, keyed by the name of the grism.
    """
    best = {}
    for spectrum in sorted(spectra, key=lambda spectrum: np.sum(~bad_flux_bins(spectrum)), reverse=True):
        best.setdefault(grism_name(spectrum), spectrum)
    return best


best = best_spectra(spectra)
list(best)
```

### Plot the spectra

We plot the usable flux bins of the best spectrum of each grism, and mark the strongest emission lines that fall in the wavelength range of the usable flux bins at the spectroscopic redshift, to see how well the photometric redshift agrees with it:

```{code-cell} ipython3
emission_lines = {'Hβ': 4862.7 * u.AA, '[OIII]': 5008.2 * u.AA, 'Hα': 6564.6 * u.AA}
grism_colors = {'Red grism': 'tab:red', 'Blue grism': 'tab:blue'}

fig, ax = plt.subplots(figsize=(12, 4.5))
for name, spectrum in best.items():
    bad = bad_flux_bins(spectrum)
    # Draw the red grism spectrum on top, as it is less noisy where the two grisms overlap
    ax.plot(spectrum['WAVELENGTH'].to(u.micron), np.ma.masked_where(bad, spectrum['SIGNAL']),
            color=grism_colors[name], alpha=0.8, label=name, zorder=3 if name == 'Red grism' else 2)

usable_wavelengths = [spectrum['WAVELENGTH'][~bad_flux_bins(spectrum)] for spectrum in best.values()]
wave_min = min(wavelengths.min() for wavelengths in usable_wavelengths)
wave_max = max(wavelengths.max() for wavelengths in usable_wavelengths)
for line, rest_wavelength in emission_lines.items():
    observed = (rest_wavelength * (1 + spec_z)).to(u.micron)
    if wave_min < observed < wave_max:
        ax.axvline(observed, color='grey', linestyle='--', linewidth=0.8)
        ax.text(observed.value, 0.98, f' {line}', transform=ax.get_xaxis_transform(), rotation=90, va='top', ha='right')

ax.set_title(f'Object {obj_id}: phz_median = {photo_z:.2f}, spe_z = {spec_z:.3f}')
ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize='small')
plt.show()
```

This galaxy has spectra from both grisms.
At its spectroscopic redshift, which agrees well with the photometric redshift, the Hα line falls at 1.30 micron, where the two grisms overlap, and it is seen in both spectra.
This bright galaxy has a strong continuum, and its emission lines are relatively weak: the Hβ and [OIII] lines, near 1 micron in the blue grism spectrum, are not clearly detected.

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

<!-- TODO: update the Hartley et al. and Enia et al. references (journal, volume, arXiv ID) once the papers are public. -->

- Euclid Collaboration: Hartley, W., Paltani, S., Dubath, F., et al. 2026, *Euclid Data Release 1 (DR1): Photometric-redshift pipeline overview and validation*, A&A, in prep.
- Euclid Collaboration: Enia, A., Hartley, W. G., Bolzonella, M., et al. 2026, *Euclid Data Release 1 (DR1): The PHZ processing function for the physical properties of galaxies*, A&A, in prep.
- Euclid Collaboration: Tucci, M., Paltani, S., Hartley, W. G., et al. 2026, A&A, 711, A5, which describes the PHZ processing of Q1.

+++

## About this Notebook

**Updated**: 2026-11-12

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
