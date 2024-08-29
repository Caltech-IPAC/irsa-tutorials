# Caltech/IPAC -- IRSA notebooks


These Jupyter Notebook tutorials demonstrate access methods and techniques for working with data served by IRSA.
They cover topics like querying IRSA, working with catalogs in Parquet format, and general parallelization techniques.


## Accessing IRSA's on Premises Holding using VO protocols

These notebooks use the tools: PyVO, Astropy.

```{toctree}
---
maxdepth: 1
caption: IRSA's on premises
---
tutorials/irsa-sia-examples/sia_2mass_allsky
tutorials/irsa-sia-examples/sia_allwise_atlas
tutorials/irsa-sia-examples/sia_cosmos
tutorials/irsa-sia-examples/siav2_seip
tutorials/cosmodc2/cosmoDC2_TAP_access.md

```

## Accessing IRSA's Cloud Holdings

These notebooks demonstrate basic access to the IRSA-curated datasets available in AWS S3 cloud storage buckets.
They show examples for the Parquet version of the AllWISE Source Catalog, located in AWS S3 cloud storage and the exploration of simulated Roman observations store in AWS S3 cloud storage.

They use the tools: Pandas, PyArrow, Astropy, Astroquery, PyVO, S3FS, matplotlib


```{toctree}
---
maxdepth: 1
caption: IRSA in the cloud
---

tutorials/cloud_access/cloud-access-intro
tutorials/parquet-catalog-demos/wise-allwise-catalog-demo
tutorials/openuniversesims/openuniverse2024_roman_simulated_timedomainsurvey
tutorials/openuniversesims/openuniverse2024_roman_simulated_wideareasurvey

```


## Interactive Visualization in Python with Firefly

```{toctree}
---
maxdepth: 1
caption: Visualizations with Firefly
---

tutorials/firefly/SEDs_in_Firefly
tutorials/firefly/NEOWISE_light_curve_demo

```


## Generally Useful Techniques

```{toctree}
---
maxdepth: 1
caption: Generic Techniques
---

tutorials/parallelize/Parallelize_Convolution

```
