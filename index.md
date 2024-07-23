# Caltech/IPAC -- IRSA notebooks


These Jupyter Notebook tutorials demonstrate access methods and techniques for working with data served by IRSA.
They cover topics like querying IRSA, working with catalogs in Parquet format, and general parallelization techniques.


## Accessing IRSA archive holdings

### Image Thumbnailes

These notebooks show how to query IRSA's image services, inspect the results, and download and visualize images.

They use the tools: PyVO, Astropy.

```{toctree}
---
maxdepth: 1
---

tutorials/irsa-sia-examples/sia_2mass_allsky
tutorials/irsa-sia-examples/sia_allwise_atlas
tutorials/irsa-sia-examples/sia_cosmos
tutorials/irsa-sia-examples/siav2_seip

```

<!---
### Catalogs

```{toctree}
---
maxdepth: 1
---


```

## Visualizations
```{toctree}
---
maxdepth: 1
---

```
-->


## IRSA in the cloud

This notebook demonstrates basic access to the IRSA-curated datasets available in AWS S3 cloud storage buckets.

It uses the tools: Pandas, PyArrow, Astropy, Astroquery, PyVO, S3FS

```{toctree}
---
maxdepth: 1
---

tutorials/cloud_access/cloud-access-intro
```

### Catalogs

This notebook shows examples for the Parquet version of the AllWISE Source Catalog, located in AWS S3 cloud storage.
It uses the tools: Pandas, PyArrow, Astropy


```{toctree}
---
maxdepth: 1
---

tutorials/parquet-catalog-demos/wise-allwise-catalog-demo

```

### Explore OpenUniverse 2024 Data Preview

These notebooks explore simulared Roman observation stored in AWS S3 cloud storage.
They use the tools: Pandas, Astropy, S3FS, matplotlib, NumPy

```{toctree}
---
maxdepth: 1
---

tutorials/openuniversesims/openuniverse2024_roman_simulated_timedomainsurvey
tutorials/openuniversesims/openuniverse2024_roman_simulated_wideareasurvey

```

## Generally useful techniques

```{toctree}
---
maxdepth: 1
---

tutorials/parallelize/Parallelize_Convolution

```
