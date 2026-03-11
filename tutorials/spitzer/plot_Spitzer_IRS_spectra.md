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
---

# Plot Spitzer IRS Spectra

## Learning Goals

By the end of this tutorial, you will be able to :

- Discover which spectroscopy catalogs are available through IRSA’s Simple Spectral Access (SSA) service
- Use Irsa.query_ssa() to search for spectra at specific sky positions
- Retrieve and visualize infrared spectra from the IRSA archive
- Adapt the workflow shown here to your own science use cases
  
## Introduction

In this tutorial we use a published sample of debris disk host stars to demonstrate how to search for and retrieve infrared spectra using IRSA’s SSA service. 
The targets are drawn from the catalog of debris disks compiled by Mittal et al. (2015, ApJ, 798, 87), which identifies stars exhibiting infrared excesses indicative of circumstellar dust. 
Debris disks trace the remnants of planet formation and are valuable probes of disk evolution, dust composition, and dynamical interactions with planets.

Infrared spectroscopy is particularly powerful for studying debris disks because many of the diagnostic features of dust grains—such as silicate emission bands—appear at mid‑infrared wavelengths. 
The Spitzer Infrared Spectrograph (IRS) provides sensitive, low‑ and moderate‑resolution spectra in this regime, enabling measurements of dust temperature, composition, and structure. 
Querying the IRSA archive for Spitzer IRS observations allows us to connect published debris‑disk samples with archival spectroscopy and explore these systems in greater physical detail.

### Instructions

Feel free to adapt this notebook to your own targets and science goals. 
The functions and overall structure shown here are intended to be reusable for searching IRSA for spectra of interest using SSA.

### Input

- A list of sky positions (RA, Dec)

### Output

- Tables describing available spectra
Optional plots of flux versus wavelength

## Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install astropy "astroquery>=0.4.10"
```

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table

from astroquery.vizier import Vizier
from astroquery.ipac.irsa import Irsa

import warnings
#suppress warnings about cache 
warnings.filterwarnings(
    "ignore",
    message="XDG_CACHE_HOME is set",
    category=UserWarning,
)
```

## 1. Exploring SSA catalogs available at IRSA
Before querying for spectra, it is often useful to see which spectral collections are available through SSA. IRSA hosts spectra from multiple missions and instruments, each grouped into collections.

```{code-cell} ipython3
Irsa.list_collections()
```

Each entry corresponds to a distinct spectral data collection (for example, Spitzer IRS enhanced products). You can use these collection names with query_ssa(collection=...) to control which archive holdings are searched. Users interested in other instruments or wavelength ranges are encouraged to explore this list and substitute a different collection name below.

+++

## 2. Define the targets

We begin by loading a published catalog of debris disk host stars from VizieR.

```{code-cell} ipython3
vizier = Vizier() # this instantiates Vizier with its default parameters
vizier.ROW_LIMIT = 150
# VizieR catalog identifier for Mittal et al. (2015)
mittal = "J/ApJ/798/87"

debris_disks = vizier.get_catalogs(mittal)
```

### 2.1 Explore the target catalog

Let’s inspect the contents of the returned catalog and identify the columns we need.

```{code-cell} ipython3
debris_disks
```

```{code-cell} ipython3
#we really want the first table, not the refs table
debris_disks = debris_disks[0]
debris_disks
```

```{code-cell} ipython3
debris_disks.colnames
```

## 3. Find IRSA spectroscopy

Lets see if any of these debris disks have spectra in the IRSA archive

```{code-cell} ipython3
# IRSA queries require sky coordinates, 
# so we convert the RA and Dec columns into a vectorized SkyCoord object.
coords = SkyCoord(ra=debris_disks["_RA"],
                           dec=debris_disks["_DE"],
                           unit=u.deg,
                           frame='icrs')

# Just to make this tutorial run faster, we will limit the number of debris disks
coords = coords[0:10]
```

The function below queries IRSA’s SSA service for Spitzer IRS enhanced spectra near each target position. It optionally plots the retrieved spectra and includes inline comments explaining each step.

We provide this as a function so it can be easily lifted from this tutorial and used in your own work.

```{code-cell} ipython3
def query_and_plot_spectra(positions, *, plot=True, verbose = True):
    """
    Query IRSA for Spitzer IRS spectra near each target position using SSA.

    
    Parameters
    ----------
    positions : astropy.coordinates.SkyCoord
        Vectorized SkyCoord object containing target positions.
    plot : bool, optional
        If True, generate plots of flux versus wavelength.
    verbose : bool, optional
        If True, print status messages about query results.
    """
    #for each set of coordinates
    for i, sc in enumerate(positions):
        # Retrieve the target name from the debris disk catalog
        target_name = debris_disks["Name"][i]
        
        # Query the IRSA SSA service around this position
        result = Irsa.query_ssa(pos=sc, 
                                radius=5*u.arcsec,
                                collection='spitzer_irsenh')

        # Handle cases with no available spectra
        if result is None or len(result) == 0:
            if verbose:
                print(f"No IRS spectra available for target {target_name}")
                continue

        # Let the user know we have a winner        
        if verbose:
            print(f"Found {len(result)} spectrum(s) for {target_name}")

        # Loop through each spectrum returned for the object
        for j, row in enumerate(result):
            
            # Each SSA row includes an access URL pointing to the spectrum
            spectrum_url = row['access_url']
            
            # Read the spectrum into an Astropy Table
            single_spec = Table.read(spectrum_url, format="ipac")
                                                    
            # If plotting is enabled, plot the spectrum.
            if plot:
                plt.figure()
                if len(result) > 1:
                    label = row['ObsID'] if 'ObsID' in row.colnames else f"Spectrum {j+1}"
                    plt.plot(single_spec['wavelength'], single_spec['flux_density'], label=label)
                else:
                    plt.plot(single_spec['wavelength'], single_spec['flux_density'])
             
                plt.xlabel("Wavelength (μm))")
                plt.ylabel("Flux ")
                plt.title(f"Spitzer IRS Spectrum for {target_name}")
            
                # If more than one spectrum was found, add a legend
                if len(result) > 1:
                    plt.legend()
                
                plt.tight_layout()
                plt.show()
```

```{code-cell} ipython3
query_and_plot_spectra(coords)
```

## Acknowledgements

- [IPAC-IRSA](https://irsa.ipac.caltech.edu/)

## About this notebook

**Authors:** IPAC Science Platform Team, including Troy Raen, Brigitta Sipőcz, Jessica Krick, Andreas Faisst, Vandana Desai

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions
or problems.

**Updated:** 2026-01-13

**Runtime:** As of the date above, this notebook takes about 3 minutes to run to completion on a machine with 8GB RAM and 2 CPU. This runtime is
heavily dependent on archive servers which means runtime will vary for users".)

**AI:** AI-based tools were used to assist in drafting and refining this tutorial, with all content reviewed and validated by the authors.
