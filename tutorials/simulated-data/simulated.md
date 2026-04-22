# Simulated Data Tutorial Notebooks

IRSA hosts a diverse collection of simulated astronomical datasets spanning multiple missions and science domains; designed to support survey planning, algorithm development, and scientific exploration.
Because this collection is heterogeneous in coverage, structure, and intended use, the simulated products are released with detailed documentation.
Access methods are tailored to the structure and scale of each product.
These tutorials are designed to help users get started with accessing, visualizing, and analyzing simulated datasets hosted at IRSA.

```{notebook-gallery} notebook_metadata.yml
tutorials/simulated-data/roman_hlss_number_density.md
tutorials/simulated-data/cosmoDC2_TAP_access.md
tutorials/simulated-data/OpenUniverse2024/quickstart.md
tutorials/simulated-data/OpenUniverse2024/OpenUniverse2024Preview_Firefly.md
tutorials/simulated-data/OpenUniverse2024/openuniverse2024_roman_simulated_wideareasurvey.md
tutorials/simulated-data/OpenUniverse2024/TDE_light_curve.md
tutorials/simulated-data/OpenUniverse2024/SED_fit.md
tutorials/simulated-data/OpenUniverse2024/openuniverse2024_roman_simulated_timedomainsurvey.md
```

## Roman High Latitude Spectroscopic Survey — Number Density as a Function of Redshift

Demonstrates how to query the simulated Roman HLSS catalog and derive galaxy number density as a function of redshift, illustrating large-scale structure and the distribution of matter over cosmic time.

## CosmoDC2 — Querying Mock Catalogs via TAP

Demonstrates how to access and query the CosmoDC2 Mock v1 catalogs using IRSA's TAP service and the PyVO library, with queries written in ADQL.

## OpenUniverse2024

The [OpenUniverse2024](https://arxiv.org/abs/2501.05632) dataset is a suite of large-scale cosmological simulations designed to support joint survey planning between the Nancy Grace Roman Space Telescope and the Vera C. Rubin Observatory (LSST). Covering roughly 70 square degrees of sky with matched optical and infrared imaging, it provides realistic synthetic catalogs and images incorporating detailed extragalactic modeling, transient populations, and instrument effects.

The following six tutorials cover the OpenUniverse2024 dataset — from data access and visualization to photometric analysis and spectral energy distribution (SED) fitting.

### Quickstart — Accessing OpenUniverse2024 Data

A focused introduction to the three main data access patterns: browsing the S3 directory structure for Roman and Rubin FITS images, reading the parquet catalogs (transient, galaxy, and galaxy-flux tables), and querying which images cover a given sky position using the IRSA Simple Image Access (SIA) service.

### Firefly Visualization — Exploring OpenUniverse2024 with Firefly

Shows how to access cloud-hosted Roman and Rubin simulated images and visualize them interactively using the Firefly JupyterLab extension, including overplotting DS9 regions, Parquet catalogs, and creating 3-color composites.

### Roman Coadds — Analyzing Simulated Roman Coadded Images

Demonstrates how to browse the S3 bucket containing simulated Roman wide-area coadded images, identify a coadd by filter, row, and column given a sky coordinate, and inspect the structure of the resulting FITS file.

### TDE Light Curve — Building a Multi-Epoch Light Curve

Demonstrates an end-to-end science workflow for transient astronomy: locating a simulated Tidal Disruption Event (TDE) in the OpenUniverse2024 transient catalog, identifying its host galaxy, retrieving Roman images via the IRSA SIA service, and performing aperture photometry to construct a multi-epoch light curve.

### SED Fitting — Fitting Galaxy SEDs with Prospector

Demonstrates how to build a full science workflow: starting from OpenUniverse2024 photometric catalogs, constructing spectral energy distributions (SEDs), and fitting them using the Prospector Bayesian SED fitting code. This example focuses on supernova host galaxies, comparing stellar populations between Type Ia and core-collapse supernovae.

### Roman Time Domain Survey — Supernovae Cutouts and Animations

Introduces the simulated Roman TDS preview data, demonstrates how to locate simulated supernovae, create aligned image cutouts, and produce animated GIFs showing transient evolution across epochs.
