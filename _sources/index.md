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
tutorials/parquet-catalog-demos/wise-allwise-catalog-demo
tutorials/parquet-catalog-demos/neowise-source-table-strategies
tutorials/parquet-catalog-demos/neowise-source-table-lightcurves
tutorials/openuniversesims/openuniverse2024_roman_simulated_timedomainsurvey
tutorials/openuniversesims/openuniverse2024_roman_simulated_wideareasurvey

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