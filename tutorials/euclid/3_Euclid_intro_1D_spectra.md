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

# Euclid DR1: SIR 1D Spectra

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- know how to find a galaxy and its spectroscopic redshift in the DR1 MER and SPE catalogs with IRSA's Table Access Protocol (TAP) service;
- know how to find and retrieve the spectra of a galaxy from its object ID;
- know how to use the MASK column to pick the flux bins to plot;
- understand the difference between the red and blue grism spectra, and how to pick the best spectrum of each grism;
- know how to plot a Euclid spectrum and mark its emission lines;
- know how to explore the spectra of randomly selected galaxies.

+++

## Introduction

+++

Euclid launched in July 2023 as a European Space Agency (ESA) mission with involvement by NASA.
The primary science goals of Euclid are to better understand the composition and evolution of the dark Universe.
The Euclid mission is providing space-based imaging and spectroscopy as well as supporting ground-based imaging to achieve these primary goals.
These data are archived by multiple global repositories, including IRSA, where they will support transformational work in many areas of astrophysics.

Euclid Data Release 1 (DR1) comprises about 1900 sq deg of imaging, spectroscopy, and catalogs from the Euclid Wide Survey (EWS), which will eventually cover about 14,000 sq deg, and from the Euclid Deep Survey (EDS), which covers about 53 sq deg in the three Euclid Deep Fields (North, South, and Fornax).
It follows Quick Release 1 (Q1), which covered about 63 sq deg with a single pass over each of the three Euclid Deep Fields and the dark cloud LDN 1641; all of the DR1 data, including the areas observed in Q1, have been reprocessed with updated pipelines and calibrations.
Further information about Euclid DR1 can be found at https://euclid.caltech.edu/page/euclid-dr1.

Among the data products included in DR1 are the one-dimensional (1D) spectra extracted from the NISP slitless spectroscopic observations by the SIR (Spectroscopic Imaging Reduction) processing function.
The SIR processing for DR1 is described in detail in Euclid Collaboration: Scodeggio et al. (2026), see the [References](#references) section below.
This notebook provides an introduction to the Euclid DR1 SIR 1D spectra.
If you have questions about it, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## Imports

```{important}
We rely on ``astroquery`` features that have been recently added, so please make sure you have version v0.4.10 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy matplotlib astropy 'astroquery>=0.4.10'
```

```{code-cell} ipython3
import warnings
from urllib.error import HTTPError

import numpy as np

import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.utils.data import conf
from astropy.visualization import quantity_support
from astropy import units as u

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
```

## 1. Find a galaxy and its redshift in the DR1 catalogs

We start with a galaxy in the Euclid Deep Field North.
To find its spectra, we need its Euclid object ID, and to identify the emission lines in the spectra, we need its redshift.
Both come from the DR1 catalogs: the MER catalog gives the object ID of the source at a position, and the SPE catalog gives the spectroscopic redshift measured from its spectra (see the [MER catalogs notebook](2_Euclid_intro_MER_catalog.md) and the [SPE catalogs notebook](5_Euclid_intro_SPE_catalog.md)).
The SPE catalog has no coordinates, so we join it with the MER catalog on `object_id`.
Our galaxy is in a deep tile, so we use the deep versions of both catalogs, and keep the best redshift estimate (`spe_rank = 0`):

```{code-cell} ipython3
coord = SkyCoord(272.3240, 68.2962, unit='deg', frame='icrs')
search_radius = 1 * u.arcsec
```

```{code-cell} ipython3
def find_galaxy(coord, search_radius=1 * u.arcsec):
    query = f"""
        SELECT mer.object_id, mer.ra, mer.dec, spe.spe_z, spe.spe_z_prob
        FROM euclid_dr1_mer_catalogue_deep_survey AS mer
        JOIN euclid_dr1_spectro_z_spe_galaxy_candidates_deep AS spe ON spe.object_id = mer.object_id
        WHERE spe.spe_rank = 0
          AND CONTAINS(POINT('ICRS', mer.ra, mer.dec),
                       CIRCLE('ICRS', {coord.ra.deg}, {coord.dec.deg}, {search_radius.to(u.deg).value})) = 1
    """
    return Irsa.query_tap(query, async_job=True).to_table()


galaxy = find_galaxy(coord)
galaxy
```

```{code-cell} ipython3
object_id = galaxy['object_id'][0]
redshift = galaxy['spe_z'][0]
```

## 2. Find the spectra of the galaxy

The DR1 SIR 1D spectra are stored in large FITS files, one per MER tile, with one extension (HDU) per spectrum.
The table `euclid.objectid_spectrafile_association_dr1` at IRSA lists, for each object ID, the file (`path`) and the HDU of each of its spectra.
IRSA's spectrumDM service converts a given HDU of a spectrum file into a VOTable; its URL is built from the file path and HDU, which we do directly in the query:

```{code-cell} ipython3
def find_spectra(object_id):
    query = f"""
        SELECT objectid, tileid, hdu,
               'https://irsacloud.ipac.caltech.edu/api/spectrumdm/convert/' || path
               || '?dataset_id=euclid_combspec_dr1&hdu=' || hdu AS url
        FROM euclid.objectid_spectrafile_association_dr1
        WHERE objectid = {object_id}
    """
    return Irsa.query_tap(query).to_table()


spectrum_files = find_spectra(object_id)
spectrum_files
```

```{note}
IRSA's Simple Spectral Access (SSA) service will also serve the DR1 SIR 1D spectra, which will allow searching for them by position; this notebook will be updated when it is available.
```

## 3. Read in the spectra and plot one of them

Each URL returns a VOTable containing the spectrum, which we read into an astropy table:

```{code-cell} ipython3
spectra = [QTable.read(url) for url in spectrum_files['url']]
```

```{note}
Some highlights of the columns of the 1D spectra:

- WAVELENGTH is in Angstroms by default.
- SIGNAL is the flux. The values are scaled and the scaling factor is included in the column's units. This value corresponds to the `'FSCALE'` entry in the HDU header of the original FITS file.
- VAR is the variance of the flux, and UNCERTAINTY is its square root, the flux uncertainty.
- MASK values can be used to determine which flux bins to discard: flux bins with an odd MASK value or a MASK value of 64 or more should not be used.
- NDITH is the number of dithered exposures that contribute to each flux bin.
```

Before plotting, we use the MASK column to find the flux bins not to use:

```{code-cell} ipython3
def bad_flux_bins(spectrum):
    mask = spectrum['MASK'].value
    return (mask % 2 == 1) | (mask >= 64)
```

Let's plot the spectrum with the most usable flux bins, with the flux as a solid line and its uncertainty as a thin line:

```{code-cell} ipython3
spectrum = max(spectra, key=lambda spectrum: np.sum(~bad_flux_bins(spectrum)))

bad = bad_flux_bins(spectrum)
wavelength = spectrum['WAVELENGTH'].to(u.micron)

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(wavelength, np.ma.masked_where(bad, spectrum['SIGNAL']), color='tab:red', label='Flux')
ax.plot(wavelength, np.ma.masked_where(bad, spectrum['UNCERTAINTY']), color='tab:red', linewidth=0.5, label='Uncertainty')
ax.set_title(f'Euclid SIR 1D spectrum of object {object_id}')
ax.legend()
plt.show()
```

The strong emission line near 1.76 micron is Hα, and the lines near 1.3 micron are Hβ and [OIII]; we will confirm this with the redshift of the galaxy in section 5.

## 4. Red and blue grism spectra

NISP takes slitless spectra with two kinds of grisms:

* **Red grisms** cover about 1.2 to 1.9 micron and are used in both the Euclid Wide Survey and the Euclid Deep Survey.
* **The blue grism** covers about 0.9 to 1.4 micron and is used only in the Euclid Deep Survey.

A galaxy in the Euclid Wide Survey has only a red grism spectrum, while a galaxy in the Euclid Deep Fields, like ours, has spectra from both kinds of grism, each of them a separate spectrum in the files.
We tell them apart by their wavelength range (blue grism spectra start below 1.1 micron):

```{code-cell} ipython3
def grism_name(spectrum):
    return 'Blue grism' if spectrum['WAVELENGTH'].min() < 1.1 * u.micron else 'Red grism'
```

Let's compare the spectra we found: their grism, wavelength range, the largest number of dithered exposures in a flux bin, and the number of usable flux bins:

```{code-cell} ipython3
summary = QTable({
    'grism': [grism_name(spectrum) for spectrum in spectra],
    'min_wavelength': [spectrum['WAVELENGTH'].min().to(u.micron) for spectrum in spectra],
    'max_wavelength': [spectrum['WAVELENGTH'].max().to(u.micron) for spectrum in spectra],
    'max_ndith': [spectrum['NDITH'].max() for spectrum in spectra],
    'usable_bins': [np.sum(~bad_flux_bins(spectrum)) for spectrum in spectra],
})
summary.sort(['grism', 'usable_bins'], reverse=True)
summary
```

The spectrum we plotted above is the red grism one.
A galaxy can also have several spectra of the same kind of grism, for example when it is covered by several tiles; these can differ in depth, so in general we keep, for each kind of grism, the spectrum with the most usable flux bins:

```{code-cell} ipython3
def best_spectra(spectra):
    best = {}
    for spectrum in sorted(spectra, key=lambda spectrum: np.sum(~bad_flux_bins(spectrum)), reverse=True):
        best.setdefault(grism_name(spectrum), spectrum)
    return best


best = best_spectra(spectra)
list(best)
```

## 5. Plot the best spectra and mark the emission lines

These are the rest-frame vacuum wavelengths of some of the strongest emission lines of galaxies; at the redshift of the galaxy, they are shifted by a factor of (1 + z).
Hα is blended with the [NII] lines, and the [SII] and [OII] doublets are not resolved at the resolution of the Euclid spectra, so we mark each of them with a single label.

```{code-cell} ipython3
emission_lines = {
    '[OII]': 3728.5 * u.AA,
    '[NeIII]': 3870.9 * u.AA,
    'Hγ': 4341.7 * u.AA,
    'Hβ': 4862.7 * u.AA,
    '[OIII]': 5008.2 * u.AA,
    'Hα': 6564.6 * u.AA,
    '[SII]': 6725.5 * u.AA,
}
```

We now wrap the plotting into a function that plots the best spectrum of each grism.
We plot the flux bins to use as solid lines, the flux bins not to use as dotted lines, and the flux uncertainty as thin lines.
We set the y-axis range from the flux bins to use, so that the bad flux bins do not stretch it, and mark the emission lines that fall in the wavelength range of the spectra.

```{code-cell} ipython3
grism_colors = {'Red grism': 'tab:red', 'Blue grism': 'tab:blue'}


def plot_spectra(spectra, title, redshift=None):
    fig, ax = plt.subplots(figsize=(12, 4.5))
    for name, spectrum in spectra.items():
        wavelength = spectrum['WAVELENGTH'].to(u.micron)
        bad = bad_flux_bins(spectrum)
        color = grism_colors[name]

        ax.plot(wavelength, np.ma.masked_where(bad, spectrum['SIGNAL']), color=color, alpha=0.8, label=name)
        ax.plot(wavelength, np.ma.masked_where(~bad, spectrum['SIGNAL']), color=color, linestyle=':', alpha=0.6)
        ax.plot(wavelength, spectrum['UNCERTAINTY'], color=color, linewidth=0.5, alpha=0.5)

    good_flux = np.concatenate([spectrum['SIGNAL'][~bad_flux_bins(spectrum)].filled(np.nan)
                                for spectrum in spectra.values()])
    low, high = np.nanpercentile(good_flux, [1, 100])
    ax.set_ylim(low - 0.1 * (high - low), high + 0.25 * (high - low))

    # Mark the emission lines at the redshift of the galaxy
    if redshift is not None:
        wave_min = min(spectrum['WAVELENGTH'].min() for spectrum in spectra.values())
        wave_max = max(spectrum['WAVELENGTH'].max() for spectrum in spectra.values())
        for line, rest_wavelength in emission_lines.items():
            observed = (rest_wavelength * (1 + redshift)).to(u.micron)
            if wave_min < observed < wave_max:
                ax.axvline(observed, color='grey', linestyle='--', linewidth=0.8)
                ax.text(observed.value, 0.98, f' {line}', transform=ax.get_xaxis_transform(),
                        rotation=90, va='top', ha='right', fontsize=9)
        title = f'{title}, z = {redshift:.3f}'

    # Legend entries for the line styles
    ax.plot([], [], color='grey', linestyle=':', label='Do not use')
    ax.plot([], [], color='grey', linewidth=0.5, label='Uncertainty')

    ax.set_title(title)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize='small')
    plt.show()


plot_spectra(best, f'Euclid SIR 1D spectra of object {object_id}', redshift)
```

## 6. Look closer at the spectrum table and the MASK column

Now that we have seen the spectra, let's look at the usable flux bins of the red grism spectrum:

```{code-cell} ipython3
red_spectrum = best['Red grism']
red_spectrum[~bad_flux_bins(red_spectrum)]
```

The MASK values are bit flags, so each value combines several conditions.
Let's plot them as a function of wavelength, to see where the flux bins not to use are:

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(12, 3))
for name, spectrum in best.items():
    ax.plot(spectrum['WAVELENGTH'].to(u.micron), spectrum['MASK'], color=grism_colors[name], alpha=0.6, label=name)
ax.set_ylabel('MASK value')
ax.set_title('Values of MASK by flux bin')
ax.legend()
plt.show()
```

Most of the flux bins not to use are at the edges of the wavelength range of each grism, where its throughput is low.

## 7. Explore the spectra of random galaxies

Now it's your turn to explore!
We select galaxies from the DR1 SPE and MER deep catalogs in a region of the Euclid Deep Field North, pick a few of them at random, and plot their spectra with their emission lines.
To avoid stars and problematic sources, we keep only sources with a reliable redshift (`spe_z_prob` and `spe_z_rel`, see the [SPE catalogs notebook](5_Euclid_intro_SPE_catalog.md)), a low probability of being point-like (`point_like_prob`), and no MER detection quality flags set (`det_quality_flag = 0`), for example for blending with a neighbor.
Change the parameters below to look at other regions, redshift ranges, or numbers of galaxies; change or remove the random seed to get a different selection.
Each time you rerun the plotting cell, it picks new galaxies from the same selection.
Most of these galaxies are much fainter than our first example, so their spectra are noisier, and their emission lines can be weak or missing.

```{code-cell} ipython3
######################## User defined section ############################
# Center and radius of the region to select galaxies from
region_center = SkyCoord(271.0, 66.0, unit='deg', frame='icrs')
region_radius = 6 * u.arcmin

# Redshift range, the minimum probability and reliability of the redshift estimate,
# and the maximum probability of the source being point-like
z_min, z_max = 0.9, 1.8
min_z_prob = 0.99
min_z_rel = 0.05
max_point_like_prob = 0.5

# Number of galaxies to plot, and the random seed
n_galaxies = 3
rng = np.random.default_rng(seed=42)
##########################################################################

query = f"""
    SELECT mer.object_id, mer.ra, mer.dec, spe.spe_z, spe.spe_z_prob, spe.spe_z_rel
    FROM euclid_dr1_mer_catalogue_deep_survey AS mer
    JOIN euclid_dr1_spectro_z_spe_galaxy_candidates_deep AS spe ON spe.object_id = mer.object_id
    WHERE spe.spe_rank = 0
      AND spe.spe_z BETWEEN {z_min} AND {z_max}
      AND spe.spe_z_prob > {min_z_prob}
      AND spe.spe_z_rel > {min_z_rel}
      AND mer.point_like_prob < {max_point_like_prob}
      AND mer.det_quality_flag = 0
      AND CONTAINS(POINT('ICRS', mer.ra, mer.dec),
                   CIRCLE('ICRS', {region_center.ra.deg}, {region_center.dec.deg}, {region_radius.to(u.deg).value})) = 1
"""
galaxies = Irsa.query_tap(query, async_job=True).to_table()
print(f'Found {len(galaxies)} galaxies')
```

We wrap sections 2 to 4 into a function, which returns the best spectrum of each grism of an object.
If a spectrum file cannot be read, we skip it and use the other spectra of the galaxy:

```{code-cell} ipython3
def get_best_spectra(object_id):
    spectra = []
    for url in find_spectra(object_id)['url']:
        try:
            spectra.append(QTable.read(url))
        except HTTPError as error:
            print(f'Skipping a spectrum that could not be read ({error}): {url}')
    return best_spectra(spectra)
```

```{code-cell} ipython3
for galaxy in galaxies[rng.choice(len(galaxies), size=min(n_galaxies, len(galaxies)), replace=False)]:
    galaxy_spectra = get_best_spectra(galaxy['object_id'])
    if not galaxy_spectra:
        print(f"No spectra found for object {galaxy['object_id']}")
        continue
    plot_spectra(galaxy_spectra, f"Object {galaxy['object_id']}", galaxy['spe_z'])
```

## References

<!-- TODO: update the Scodeggio et al. reference (journal, volume, arXiv ID) once the paper is public. -->

- Euclid Collaboration: Scodeggio, M., et al. 2026, *Euclid Data Release 1 (DR1): The SIR spectroscopic processing*, A&A, in prep.

+++

## About this Notebook

**Updated**: 2026-11-12

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
