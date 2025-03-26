# Caltech/IPAC--IRSA Python Notebook Tutorials


These Python Jupyter Notebook tutorials demonstrate access methods and techniques for working with data served by the [NASA/IPAC Infrared Science Archive (IRSA)](https://irsa.ipac.caltech.edu).
They cover topics like querying IRSA, working with catalogs in Parquet format, visualizing with Firefly, and general other techniques.


## Accessing IRSA's on-premises holdings using VO protocols

```{toctree}
---
maxdepth: 1
caption: VO on-prem data access
---
tutorials/irsa-sia-examples/sia_2mass_allsky
tutorials/irsa-sia-examples/sia_allwise_atlas
tutorials/irsa-sia-examples/sia_cosmos
tutorials/irsa-sia-examples/siav2_seip
tutorials/cosmodc2/cosmoDC2_TAP_access.md

```

## Accessing IRSA's cloud holdings

These notebooks demonstrate how to access the IRSA-curated datasets that available in Amazon Web Services (AWS) S3 cloud storage buckets.


```{toctree}
---
maxdepth: 1
caption: Cloud data access
---

tutorials/cloud_access/cloud-access-intro
tutorials/cloud_access/euclid-cloud-access
tutorials/parquet-catalog-demos/euclid-hats-parquet
tutorials/parquet-catalog-demos/wise-allwise-catalog-demo
tutorials/parquet-catalog-demos/neowise-source-table-strategies
tutorials/parquet-catalog-demos/neowise-source-table-lightcurves
tutorials/openuniversesims/openuniverse2024_roman_simulated_timedomainsurvey
tutorials/openuniversesims/openuniverse2024_roman_simulated_wideareasurvey

```

## Accessing Euclid data

### Euclid Early Release Observation (ERO)

```{toctree}
---
maxdepth: 1
caption: Euclid Early Release Observations
---

tutorials/euclid_access/Euclid_ERO

```

### Euclid Quick Release 1 (Q1)

```{toctree}
---
maxdepth: 1
caption: Euclid Quick Release 1
---

tutorials/euclid_access/1_Euclid_intro_MER_images
tutorials/euclid_access/2_Euclid_intro_MER_catalog
tutorials/euclid_access/3_Euclid_intro_1D_spectra
tutorials/euclid_access/4_Euclid_intro_PHZ_catalog
tutorials/euclid_access/5_Euclid_intro_SPE_catalog
tutorials/cloud_access/euclid-cloud-access
tutorials/parquet-catalog-demos/euclid-hats-parquet

```

## Interactive visualization in Python with Firefly

These notebooks demonstrate how to use the Firefly visualization tools from Python.
[Firefly](https://github.com/Caltech-IPAC/firefly) is an open-source toolkit based on IVOA standards and designed to enable astronomical data archive access, exploratory data analysis, and visualization.

It is used in archive user interfaces at [IRSA](https://irsa.ipac.caltech.edu), the [NASA Exoplanet Science Institute (NExScI)](https://nexsci.caltech.edu/), the [NASA/IPAC Extragalactic Database (NED)](https://ned.ipac.caltech.edu/), and the [Vera C. Rubin Observatory](https://www.lsst.org/).

```{toctree}
---
maxdepth: 1
caption: Visualizations with Firefly
---

tutorials/firefly/SEDs_in_Firefly
tutorials/firefly/NEOWISE_light_curve_demo
tutorials/firefly/OpenUniverse2024Preview_Firefly

```

## Simulated Data

```{toctree}
---
maxdepth: 1
caption: Simulated Data
---

tutorials/roman_simulations/roman_hlss_number_density.md

```

## Generally useful techniques

These notebooks  cover miscellaneous topics that users might find useful in their analysis of IRSA-curated data.

```{toctree}
---
maxdepth: 1
caption: Generic techniques
---

tutorials/parallelize/Parallelize_Convolution

```

***

## About these notebooks

**Authors:** IRSA Scientists and Developers wrote and maintain these notebooks.

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
