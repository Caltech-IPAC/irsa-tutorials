---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Euclid Q1 Catalogs in HATS Parquet

## Learning Goals

By the end of this tutorial, you will:

- Query the dataset to find and create figures for galaxies, QSOs, and stars with quality fluxes, redshifts, and morphology.
- Understand the format and schema of this dataset.
- Learn how to work with this HATS Parquet product using the LSDB and PyArrow python libraries.

+++

## Introduction

This notebook introduces the [Euclid Q1](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_q1.html) HATS Collection served by IPAC/IRSA and demonstrates access with python.
The Collection includes a HATS Catalog (main data product), Margin Cache (10 arcsec), and Index Table (OBJECT_ID).
The Catalog includes the twelve Euclid Q1 tables listed below, joined on the column 'OBJECT_ID' into a single Parquet dataset with 1,329 columns (one row per Euclid MER Object).
Among them, Euclid has provided several different redshift measurements, several flux measurements for each Euclid band, and flux measurements for bands from several ground-based observatories -- in addition to morphological and other measurements.
These were produced for different science goals using different algorithms and/or configurations.

Having all columns in the same dataset makes access convenient because the user doesn't have to make separate calls for data from different tables and/or join the results.
However, figuring out which, e.g., flux measurements to use amongst so many can be challenging.
In the sections below, we look at some of their distributions and reproduce figures from several papers in order to highlight some of the options and point out their differences.
The Appendix contains important information about the schema of this Parquet dataset, especially the syntax of the column names.
For more information about the meaning and provenance of a column, refer to the links provided with the list of tables below.

### Euclid Q1 tables and docs

The Euclid Q1 HATS Catalog includes the following twelve Q1 tables[*] which were products of three Euclid processing functions: MER, PHZ, and SPE.
Table names are linked to their original schemas and pointers to the Euclid papers describing them are also included.
Links to select other documentation follows.

- MER - [Romelli et al., 2025](https://arxiv.org/pdf/2503.15305) (hereafter, Romelli)
  - [MER](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#main-catalog-fits-file) - Sec. 6 & 8 (EUC_MER_FINAL-CAT)
  - [MORPH](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#morphology-catalog-fits-file) - Sec. 7 & 8 (EUC_MER_FINAL-MORPH-CAT)
  - [CUTOUTS](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#cutouts-catalog-fits-file) - Sec. 8 (EUC_MER_FINAL-CUTOUTS-CAT)
- PHZ - [Tucci et al., 2025](https://arxiv.org/pdf/2503.15306) (hereafter, Tucci)
  - [PHZ](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#photo-z-catalog) - Sec. 5 (phz_photo_z)
  - [CLASS](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#classification-catalog) - Sec. 4 (phz_classification)
  - [PHYSPARAM](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#physical-parameters-catalog) - Sec. 6 (6.1; phz_physical_parameters) _Notice that this is **galaxies** and uses a different algorithm._
  - [GALAXYSED](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#galaxy-sed-catalog) - App. B (B.1 phz_galaxy_sed)
  - [PHYSPARAMQSO](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#qso-physical-parameters-catalog) - Sec. 6 (6.2; phz_qso_physical_parameters)
  - [STARCLASS](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#star-template) - Sec. 6 (6.3; phz_star_template)
  - [STARSED](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#star-sed-catalog) - App. B (B.1 phz_star_sed)
  - [PHYSPARAMNIR](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#nir-physical-parameters-catalog) - Sec. 6 (6.4; phz_nir_physical_parameters)
- SPE - [Le Brun et al., 2025](https://arxiv.org/pdf/2503.15308) (hereafter, Le Brun)
  - [Z](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/spedpd/dpcards/spe_spepfoutputcatalog.html#redshift-catalog) - Sec. 2 (spectro_zcatalog_spe_quality, spectro_zcatalog_spe_classification, spectro_zcatalog_spe_galaxy_candidates, spectro_zcatalog_spe_star_candidates, and spectro_zcatalog_spe_qso_candidates)

See also:

- MER [Photometry](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/merphotometrycookbook.html) and [Morphology](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/mermorphologycookbook.html) Cookbooks
- [Frequently Asked Questions About Euclid Q1 data](https://euclid.caltech.edu/page/euclid-q1-data-faq) (hereafter, FAQ)
- [Q1 Explanatory Supplement](https://euclid.esac.esa.int/dr/q1/expsup/)

[*] Euclid typically calls these "catalogs", but this notebook uses "tables" to avoid any confusion with the HATS Catalog product.

### Parquet, HEALPix, and HATS

Parquet, HEALPix, and HATS are described in more detail at [https://irsadev.ipac.caltech.edu:9051/cloud_access/parquet/](https://irsadev.ipac.caltech.edu:9051/cloud_access/parquet/).
([FIXME] Currently requires IPAC VPN. Update url when the page is published to ops.)
In brief:

- [Apache Parquet](https://parquet.apache.org/docs/) is a file format that includes rich metadata and supports fast SQL-like queries on large datasets.
  - [PyArrow](https://arrow.apache.org/docs/python/) (python library) is a particularly efficient and flexible Parquet reader and we demonstrate its use below.
      Note that while PyArrow filters are more powerful than those offered by other libraries, they can be cumbersome to construct --
      For more information, we recommend IRSA's [AllWISE Parquet tutorial](https://caltech-ipac.github.io/irsa-tutorials/tutorials/parquet-catalog-demos/wise-allwise-catalog-demo.html).
- [HEALPix](https://healpix.jpl.nasa.gov/) (Hierarchical Equal Area isoLatitude Pixelization; Górski et al., 2005) is a tiling of the sky into equal-area pixels.
    - The HEALPix "order" determines the pixel area: higher order => smaller pixels, better resolution, and more total pixels required to cover the sky.
- [HATS](https://hats.readthedocs.io/) (Hierarchical Adaptive Tiling Scheme) is a HEALPix-based partitioning scheme (plus metadata) for Parquet datasets.
  - The HEALPix order at which data is partitioned is adaptive -- it varies along with the on-sky density of rows in a given catalog, with the aim of creating partitions (files) that have roughly equal numbers of rows (important for efficient access).
      The Euclid Q1 density and partitioning can be seen at [https://irsadev.ipac.caltech.edu/data/download/parquet/test/euclid/q1/catalogues/hats/skymap.png](https://irsadev.ipac.caltech.edu/data/download/parquet/test/euclid/q1/catalogues/hats/skymap.png).
      ([FIXME] Currently requires IPAC VPN. Update url when the page is published to ops. Alternately, could consider pulling this file and displaying it here. I have code to do it three different ways but all are a bit clunky and there's so much else going on in this notebook that currently it seems better not to add this.)
  - [LSDB](https://docs.lsdb.io/) is a python library specifically designed to support large-scale astronomy use cases with HATS datasets.
      It provides simple interfaces that allow users to perform (e.g.,) full-catalog cross matches with just a few lines of code.
      We demonstrate its use below.
      LSDB uses [Dask](https://docs.dask.org/) to run many core tasks.
      The Dask client configurations used in this notebook should work well for most use cases.
      But, if you run into trouble (e.g., worker memory issues) or just want to tinker and don't know where to start,
      we recommend LSDB's [Dask cluster configuration tips](https://docs.lsdb.io/en/stable/tutorials/dask-cluster-tips.html).

+++

## Installs and imports

+++

```{important}
We rely on ``lsdb>=0.5.2``, ``hpgeom>=1.4``, ``numpy>=2.0``, and ``pyerfa>=2.0.1.3`` for features that have been recently added. Please make sure that you have sufficiently recent versions installed.
```

[FIXME] Prune dependencies and imports once it's clear which libraries will be used.

```{code-cell}
# # Uncomment the next line to install dependencies if needed.
# %pip install 'lsdb>=0.5.2' 'hpgeom>=1.4' matplotlib 'numpy>=2.0' 'pyerfa>=2.0.1.3' s3fs
```

```{code-cell}
import os  # Determine number of CPUs (for parallelization)
import dask.distributed  # Parallelize catalog queries
import hpgeom
import lsdb  # Query the catalog
import matplotlib.colors  # Make figures look nice
import matplotlib.pyplot as plt  # Create figures
import numpy as np  # Math
import pandas as pd  # Manipulate query results
import pyarrow.compute as pc  # Filter dataset
import pyarrow.dataset  # Load the dataset
import pyarrow.parquet  # Load the schema
import pyarrow.fs  # Simple S3 filesystem pointer
```

```{tip}
If you run into an error that looks like,

> AttributeError: _ARRAY_API not found

or:

> A module that was compiled using NumPy 1.x cannot be run in NumPy 2.1.3 as it may crash.

make sure you have restarted the kernel since doing `pip install`. Then re-run the cell **twice**.
```

+++

## 1. Setup

+++

## 1.1 AWS S3 paths

```{code-cell}
s3_bucket = "irsa-fornax-testdata"
euclid_prefix = "EUCLID/q1/catalogues"

euclid_hats_collection_uri = f"s3://{s3_bucket}/{euclid_prefix}"  # for lsdb
euclid_parquet_metadata_path = f"{s3_bucket}/{euclid_prefix}/hats/dataset/_metadata"  # for pyarrow
euclid_parquet_schema_path = f"{s3_bucket}/{euclid_prefix}/hats/dataset/_common_metadata"  # for pyarrow

# Temporary try/except to handle credentials in different environments before public release.
try:
    # If running from within IPAC's network, your IP address acts as your credentials so this should work.
    lsdb.read_hats(euclid_hats_collection_uri)
except PermissionError:
    # If running from Fornax, credentials are provided automatically under the hood but
    # lsdb ignores them in the call above. Construct a UPath which will pick up the credentials.
    from upath import UPath

    euclid_hats_collection_uri = UPath(euclid_hats_collection_uri)
```

### 1.2 Helper functions

```{code-cell}
def magnitude_to_flux(magnitude: float) -> float:
    """Convert magnitude to flux following the MER Photometry Cookbook"""
    zeropoint = 23.9
    flux = 10 ** ((magnitude - zeropoint) / -2.5)
    return flux
```

In cases where we want magnitudes rather than fluxes, it will be convenient to have PyArrow do the conversion during the read so that we don't have to load and handle the flux columns ourselves.
The function below defines our instructions.

```{code-cell}
def flux_to_magnitude(flux_col_name: str, color_col_names: tuple[str, str] | None = None) -> pc.Expression:
    """Construct an expression for the magnitude of `flux_col_name` following the MER Photometry Cookbook.

    MAG = -2.5 * log10(TOTAL FLUX) + 23.9.
    If `color_col_names` is None, `flux_col_name` is taken as the TOTAL FLUX in the band.
    If not None, it should list the columns needed for the color correction such that
    TOTAL FLUX = flux_col_name * color_col_names[0] / color_col_names[1].

    Returns a `pyarrow.compute.Expression` which can be used in the `filter` (to filter based on it) and/or
    `columns` (to return it) keyword arguments when loading from the pyarrow dataset.
    """
    scale = pc.scalar(-2.5)
    zeropoint = pc.scalar(23.9)

    total_flux = pc.field(flux_col_name)
    if color_col_names is not None:
        color_scale = pc.divide(pc.field(color_col_names[0]), pc.field(color_col_names[1]))
        total_flux = pc.multiply(total_flux, color_scale)

    log10_flux = pc.log10(total_flux)
    mag_expression = pc.add(pc.multiply(scale, log10_flux), zeropoint)
    return mag_expression
```

### 1.3 PyArrow dataset

```{code-cell}
# Load the catalog as a PyArrow dataset to be used in several examples below.
dataset = pyarrow.dataset.parquet_dataset(euclid_parquet_metadata_path, partitioning="hive", filesystem=pyarrow.fs.S3FileSystem())
```

Define some columns to be used throughout the notebook.

```{code-cell}
# MER Object ID
OBJECT_ID = "OBJECT_ID"

# Whether the source was detected in the VIS mosaic (1) or only in the NIR mosaic (0).
VIS_DET = "MER_VIS_DET"

# Best estimate of the total flux in the detection band. From aperture photometry within a Kron radius.
# Detection band is VIS if MER_VIS_DET=1.
# Otherwise, this is a non-physical NIR-stack flux and there was no VIS detection (i.e., NIR-only).
FLUX_TOTAL = "MER_FLUX_DETECTION_TOTAL"
FLUXERR_TOTAL = "MER_FLUXERR_DETECTION_TOTAL"

# Whether the detection has a >50% probability of being spurious (1=Yes, 0=No).
SPURIOUS_FLAG = "MER_SPURIOUS_FLAG"

# Point-like morphology indicators.
POINTLIKE_PROB = "MER_POINT_LIKE_PROB"  # Always NaN for NIR-only (use MER_MUMAX_MINUS_MAG instead)
# Peak surface brightness minus magnitude in detection band.
MUMAX_MINUS_MAG = "MER_MUMAX_MINUS_MAG"  # < -2.5 => point-like; < -2.6 => compact (Tucci)

# PHZ classification, generated by a probabilistic random forest supervised ML algorithm.
PHZ_CLASS = "PHZ_PHZ_CLASSIFICATION"
PHZ_CLASS_MAP = {
    1: "Star",
    2: "Galaxy",
    4: "QSO",  # In Q1, this includes luminous AGN.
    # If multiple probability thresholds were exceeded, a combination of classes was reported.
    3: "Star and Galaxy",
    5: "Star and QSO",
    6: "Galaxy and QSO",
    7: "Star, Galaxy, and QSO",
    # Two other integers (-1 and 0) and nulls are present, indicating that the class was not determined.
    **{i: "Undefined" for i in [-1, 0, np.nan]},
}
```


Euclid Q1 includes data from three Euclid Deep Fields: EDF-N (North), EDF-S (South), EDF-F (Fornax).
There is also some data from a fourth field: LDN1641 (Lynds' Dark Nebula 1641), which was observed for technical reasons during Euclid's verification phase.
These regions are well separated, so we can distinguish them using a cone search without having to be too picky about the radius.
Rather than using the RA and Dec values directly, we'll find a set of HEALPix order 9 pixels that cover each area.
This will suffice for a simple and efficient cone search.
A column ('_healpix_9') of order 9 indexes was added to the catalog for this purpose.

[FIXME] The notebook does not currently use these but it might be good to do so.
Maybe in the Magnitudes section to show the differences as a function of class.
Anyway, either use it or remove it.

```{code-cell}
# Column name of HEALPix order 9 pixel indexes.
HEALPIX_9 = "_healpix_9"

# LDN1641 (Lynds' Dark Nebula 1641)
ra, dec, radius = 85.74, -8.39, 1.5  # need ~6 sq deg
ldn_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius)

# EDF-N (Euclid Deep Field - North)
ra, dec, radius = 269.733, 66.018, 4  # need ~20 sq deg
edfn_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius)

# EDF-S (Euclid Deep Field - South)
ra, dec, radius = 61.241, -48.423, 5  # need ~23 sq deg
edfs_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius)

# EDF-F (Euclid Deep Field - Fornax)
ra, dec, radius = 52.932, -28.088, 3  # need ~10 sq deg
edff_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius)


# ## Redshifts for cosmology

# Euclid is a cosmology mission focused on measuring the evolution of large-scale structures in order to study dark matter and dark energy.
# This means it must determine the redshifts for a large number of galaxies.
# In this section, we obtain and examine quality cosmological samples for three of the redshift point-estimates provided in Q1:
```

```{code-cell}
PHZ_Z = "PHZ_PHZ_MEDIAN"
PHYSPARAM_GAL_Z = "PHYSPARAM_PHZ_PP_MEDIAN_REDSHIFT"
SPE_GAL_Z = "Z_GALAXY_CANDIDATES_SPE_Z_RANK0"
```

"PHZ_PHZ_MEDIAN" is the median of the photometric redshift PDF that was produced for Euclid's core-science goals.
It was generated by Phosphoros, a fully Bayesian template-fitting code.
It should be better for cosmology than typical ML algorithms due to the scarcity of spectroscopic "truth" (needed to train ML) above z ~ 1.
The model grid was built spanning the parameters: redshift (z in [0, 6]), galaxy SED, intrinsic reddening curve, intrinsic attenuation.

"PHYSPARAM_PHZ_PP_MEDIAN_REDSHIFT" is the median of the photometric redshift PDF that was produced for galaxy-classed objects by the physical-properties branch of the PHZ pipeline (i.e., non-cosmology).
We include it here as a useful comparison.
It was generated by NNPZ, a k-nearest neighbors supervised learning algorithm.
The reference sample came from >1 million stellar population synthesis models spanning parameters: redshift (z in [0, 7]), age, star formation timescale, stellar metallicity, intrinsic attenuation, magnitude, and two dust laws.

"Z_GALAXY_CANDIDATES_SPE_Z_RANK0" is the top-ranked spectroscopic redshift estimate that was produced assuming the object is a galaxy (independent of the actual classification).
This is the most prominent peak of a PDF generated by a least-squares fitting algorithm (19 galaxy models).

+++

Load a quality PHZ sample. Cuts are from Tucci sec. 5.3.

```{code-cell}
# Photo-z flag: 0=good for core science, 10=NIR-only, 11=missing bands, 12=too faint.
PHZ_FLAG = "PHZ_PHZ_FLAGS"

# Columns we actually want to load.
phz_columns = [OBJECT_ID, PHZ_Z, FLUX_TOTAL, PHZ_FLAG]

# Partial filter for quality PHZ redshifts.
phz_filter = (
    (pc.field(VIS_DET) == 1)  # Detected in VIS
    & (pc.field(FLUX_TOTAL) > magnitude_to_flux(24.5))  # I < 24.5 (later, limit to 23)
    & (pc.divide(pc.field(FLUX_TOTAL), pc.field(FLUXERR_TOTAL)) > 5)  # I S/N > 5  # [CHECKME] Is this correct S/N?
    & ~(pc.field(PHZ_CLASS) == 1)  # Exclude objects classified as star.
    & (pc.field(SPURIOUS_FLAG) == 0)
)

# Execute the filter and load.
phz_df = dataset.to_table(columns=phz_columns, filter=phz_filter).to_pandas()
phz_df = phz_df.set_index(OBJECT_ID).sort_index()

# Final filter. To be applied later.
phz_final_filter = phz_df[FLUX_TOTAL] < magnitude_to_flux(23)
# 1m 6s
```

Load a quality PHYSPARAM sample. Cuts are from Tucci sec. 6.1.2.

```{code-cell}
# Properties to calculate specific star formation rate (sSFR).
PHYSPARAM_GAL_SFR = "PHYSPARAM_PHZ_PP_MEDIAN_SFR"
PHYSPARAM_GAL_MSTAR = "PHYSPARAM_PHZ_PP_MEDIAN_STELLARMASS"

# Columns we actually want to load.
pp_columns = [PHYSPARAM_GAL_Z, SPE_GAL_Z, MUMAX_MINUS_MAG, PHZ_FLAG, OBJECT_ID]

# Partial filter for quality PHYSPARAM redshifts.
# sSFR < 10^-8.2 /yr (excludes galaxies with unrealistically young ages and very high sSFR)
log_sfr_minus_log_mstar = pc.subtract(pc.field(PHYSPARAM_GAL_SFR), pc.field(PHYSPARAM_GAL_MSTAR))
# [FIXME] Does this also need to include a PHZ_CLASS cut to filter out non-galaxies?
pp_galaxy_filter = (log_sfr_minus_log_mstar < -8.2) & (pc.field(SPURIOUS_FLAG) == 0) & (pc.field("MER_DET_QUALITY_FLAG") < 4)

# Execute the filter and load.
pp_df = dataset.to_table(columns=pp_columns, filter=pp_galaxy_filter).to_pandas()
pp_df = pp_df.set_index(OBJECT_ID).sort_index()
# 1m 18s

# Final filter. To be applied later.
pp_final_filter = (pp_df[MUMAX_MINUS_MAG] > -2.6) & (pp_df[PHZ_FLAG] == 0)
```

Load a quality SPE sample. Cuts are from Le Brun sec. 3.3.

The NISP instrument was built to target Halpha emitting galaxies (z in [0.9, 1.8]).
SPE redshifts are reliable in that regime, but this represents <2% of the total spectra delivered by the SPE pipeline.
So it's crucial to make cuts in order to get it.

```{code-cell}
# SPE galaxy probability.
SPE_GAL_Z_PROB = "Z_GALAXY_CANDIDATES_SPE_Z_PROB_RANK0"

# Columns we actually want to load.
spe_columns = [SPE_GAL_Z, OBJECT_ID]

# Partial filter for quality SPE redshifts.
spe_filter = pc.field(SPE_GAL_Z_PROB) > 0.99
# [FIXME] & (pc.field(linewidth) < 680). Add when Halpha line is added to the catalog.
# Later, cut to z target range and
# sec. 6.2: Halpha flux > 2e-16 erg /s/cm2 and SN>3.5

# Execute the filter and load.
spe_df = dataset.to_table(columns=spe_columns, filter=spe_filter).to_pandas()
spe_df = spe_df.set_index(OBJECT_ID).sort_index()
# 27s

# Final filter. To be applied later.
spe_final_filter = (spe_df[SPE_GAL_Z] > 0.9) & (spe_df[SPE_GAL_Z] < 1.8)
# [FIXME] Add more when the columns are available.
```

Plot redshift distributions

```{code-cell}
tbl_colors = {"PHZ": "tab:green", "PHYSPARAM": "tab:orange", "SPE_GAL": "tab:purple"}

fig, ax = plt.subplots(1, 1, figsize=(9, 9))
# hist_kwargs = dict(bins=100, histtype="step")
hist_kwargs = dict(bins=np.linspace(0, 7.1, 100), histtype="step")

# PHZ
phz_kwargs = dict(label=PHZ_Z + " (partial quality)", color=tbl_colors["PHZ"], linestyle=":")
ax.hist(phz_df[PHZ_Z], **phz_kwargs, **hist_kwargs)
# Impose our final cuts.
phz_kwargs.update(dict(label=PHZ_Z + " (quality)"), linestyle="-")
ax.hist(phz_df.loc[phz_final_filter, PHZ_Z], **phz_kwargs, **hist_kwargs)

# PHYSPARAM
pp_kwargs = dict(label=PHYSPARAM_GAL_Z + " (partial quality)", color=tbl_colors["PHYSPARAM"], linestyle=":")
ax.hist(pp_df[PHYSPARAM_GAL_Z], **pp_kwargs, **hist_kwargs)
# Impose our final cuts.
pp_kwargs.update(label=PHYSPARAM_GAL_Z + " (quality)", linestyle="-")
ax.hist(pp_df.loc[pp_final_filter, PHYSPARAM_GAL_Z], **pp_kwargs, **hist_kwargs)

# SPE
spe_kwargs = dict(label=SPE_GAL_Z + " (partial quality)", color=tbl_colors["SPE_GAL"], linestyle=":")
ax.hist(spe_df[SPE_GAL_Z], **spe_kwargs, **hist_kwargs)
# Impose our final cuts.
spe_kwargs.update(label=SPE_GAL_Z + " (quality)", linestyle="-")
ax.hist(spe_df.loc[spe_final_filter, SPE_GAL_Z], **spe_kwargs, **hist_kwargs)

ax.set_xlabel("Redshift")
ax.set_ylabel("Counts")
plt.legend()
```

Notice that the maximum PHZ_Z redshift is z ~ 6 while PHYSPARAM_GAL_Z is z ~ 7, due to the input model parameters.
[TODO] More commenting...

+++

Next, make redshift comparison figures, treating PHZ as our ground truth.

```{code-cell}
def z_comparison_metrics(z_baseline: pd.Series, z_compare: pd.Series, ax: plt.Axes | None = None) -> dict[str, float]:
    """Calculate quality metrics for `z_baseline` vs `z_compare`. If `ax` is not None, annotate the axes."""
    residuals = (z_compare - z_baseline) / (1 + z_baseline)
    median_res = residuals.median()
    nmad = 1.4826 * (residuals - median_res).abs().median()
    outlier_frac = len(residuals.loc[residuals.abs() >= 0.15]) / len(z_baseline)
    metrics = {"Median residual": median_res, "NMAD": nmad, "Outlier fraction": outlier_frac}
    if ax:
        text = "\n".join(f"{k} = {v:.3f}" for k, v in metrics.items())
        ax.annotate(text, xy=(0.99, 0.99), xycoords="axes fraction", ha="right", va="top", bbox=dict(facecolor="w", alpha=0.8))
    return metrics
```

Compare PHZ to PHYSPARAM.
Here, we reproduce Tucci Figure 17 (left panel) except that we don't filter for EDF-F (# [FIXME] that would be easy to add).

```{code-cell}
# Get the common objects and set axes data x (PHZ) and y (PHYSPARAM).
phz_pp_df = phz_df.loc[phz_final_filter].join(pp_df.loc[pp_final_filter], how="inner", lsuffix="phz", rsuffix="pp")
x, y = phz_pp_df[PHZ_Z], phz_pp_df[PHYSPARAM_GAL_Z]
one_to_one_linspace = np.linspace(-0.01, 6, 100)

# Plot the redshift comparison and annotate with quality metrics.
fig, ax = plt.subplots(1, 1, figsize=(9, 9))
cb = ax.hexbin(x, y, bins="log")
z_comparison_metrics(x, y, ax=ax)
ax.plot(one_to_one_linspace, one_to_one_linspace, color="gray", linewidth=1)
ax.set_xlabel(PHZ_Z)
ax.set_ylabel(PHYSPARAM_GAL_Z)
plt.colorbar(cb)
```

The two photo-z estimates agree quite even though they were generated by very different different algorithms.
The two outlier clouds are roughly similar to those in Fig. 7 which were attributed to a confusion between the Balmer and Lyman breaks.

+++

Compare PHZ to SPE

```{code-cell}
# Get the common objects and set axes data x (PHZ) and y (SPE).
phz_spe_df = phz_df.loc[phz_final_filter].join(spe_df.loc[spe_final_filter], how="inner", lsuffix="phz", rsuffix="spe")
x, y = phz_spe_df[PHZ_Z], phz_spe_df[SPE_GAL_Z]
one_to_one_linspace = np.linspace(0.89, 1.81, 100)

fig, axes = plt.subplots(1, 2, figsize=(18, 9), sharey=True, width_ratios=(0.55, 0.45))
# Plot the redshift comparison and annotate with quality metrics.
ax = axes[0]
cb = ax.hexbin(x, y, bins="log")
z_comparison_metrics(x, y, ax=ax)
ax.plot(one_to_one_linspace, one_to_one_linspace, color="gray", linewidth=1)
ax.set_xlabel(PHZ_Z)
ax.set_ylabel(SPE_GAL_Z)
plt.colorbar(cb)
# Plot again, but also restrict photo-z to the cosmology target range 0.9 < z < 1.8.
ax = axes[1]
photo_z_in_range = x.loc[(x > 0.9) & (x < 1.8)].index
cb = ax.hexbin(x.loc[photo_z_in_range], y.loc[photo_z_in_range], bins="log", gridsize=25)
ax.plot(one_to_one_linspace, one_to_one_linspace, color="gray", linewidth=1)
ax.set_xlabel(PHZ_Z)
plt.colorbar(cb)
```

This should look better once the rest of the cuts (esp. for Halpha line) can be added.

+++

## Classification thresholds - purity vs completeness

+++

The PHZ_CLASS used above was determined using a threshold for galaxies which prioritized completeness over purity while the other thresholds (esp. star) gave purity more priority.
This was done in order meet certain requirements imposed by Euclid's main goals.
One way to see the affects of this are to look at the brightness and morphology as a function of class.
Here, we reproduce the first three panels (combining top and bottom) of Tucci Fig. 6.
Note that all objects in LDN1641 have a null (NaN) class designation, so our sample of stars, galaxies, and QSOs will come only from the EDF regions.

```{code-cell}
# Setup:

# Construct the filter.
class_filter = (
    (pc.field(PHZ_CLASS).isin([1, 2, 4]))  # Stars, galaxies, and QSOs (no mixed classes or Undefined).
    & (pc.divide(pc.field(FLUX_TOTAL), pc.field(FLUXERR_TOTAL)) > 5)  # S/N > 5
    & (pc.field(VIS_DET) == 1)  # Exclude NIR-only detections so that FLUX_TOTAL is VIS flux.
    # [FIXME] Are they making more cuts? Tried these, didn't help significantly.
    # & (pc.field(SPURIOUS_FLAG) == 0)
    # & (pc.field(PHZ_FLAG) == 0)
)

# Columns we actually want to load.
# Because we want PyArrow to construct a new column (magnitude), we must pass a dict mapping column names to expressions.
class_columns = {
    "I magnitude": flux_to_magnitude(FLUX_TOTAL),
    MUMAX_MINUS_MAG: pc.field(MUMAX_MINUS_MAG),
    PHZ_CLASS: pc.field(PHZ_CLASS),
}
```

```{code-cell}
# Load data.
classes_df = dataset.to_table(columns=class_columns, filter=class_filter).to_pandas()
# 30s
```

```{code-cell}
fig, axes = plt.subplots(1, 3, figsize=(18, 8))
for ax, (class_name, class_df) in zip(axes, classes_df.groupby(PHZ_CLASS)):
    ax.set_title(PHZ_CLASS_MAP[class_name])
    cb = ax.hexbin(class_df[MUMAX_MINUS_MAG], class_df["I magnitude"], bins="log")
    plt.colorbar(cb)
    ax.axvline(-2.5, color="k", linewidth=1)
    ax.set_xlabel(MUMAX_MINUS_MAG)
    ax.set_xlim(-5, 5)
    ax.set_ylabel("I magnitude")
    ax.set_ylim(15, 27)
```

Objects to the left of the vertical line are point-like.
Stars are highly concentrated there, especially those that are not faint (I < 24.5), which we should expect given Euclid's requirement for a pure sample.
Also as we should expect, most galaxies appear to the right of this line.
However, notice the strip of bright (e.g., I < 23) "galaxies" that are point-like.
Many of these are likely to be mis-classified stars or QSOs, and we could increase the purity of a galaxy sample by excluding them (at the expense of completeness, which was the original goal).
In the QSO panel, however, we see only a small cluster in the bright and point-like region where we would expect them to be.
These objects are likely to be correctly classified, and they are almost exclusively from EDF-N.
The remaining QSO classifications (the majority) should be considered doubtful.
Many QSOs are likely to be missing from the expected region due to the overlap of QSOs with galaxies in relevant color spaces and the relatively high probability thresholds imposed for QSOs.

+++

## Magnitudes

+++

Euclid Q1 contains two main types of flux measurements -- aperture and template-fit -- provided for both Euclid and external bands (e.g., DECam, PanSTARRS).
In this section, we look at the Q1 magnitude distributions in the four Euclid bands as a function of PHZ class and we compare the two types of flux measurements.
Additional flux measurements that will be of interest to some, but that we don't look at here, include: PSF-fit fluxes (VIS only); Sérsic-fit fluxes (computed for parametric morphology, discussed below); and fluxes that were corrected based on PHZ class.

```{code-cell}
# Construct a basic filter.
mag_filter = (
    (pc.field(VIS_DET) == 1)  # Exclude NIR-only detections for simplicity. (We'll inspect them later.)
    & (pc.field(PHZ_CLASS).isin([1, 2, 3, 4, 5, 6, 7]))  # Stars, Galaxies, QSOs, and mixed classes.
    & (pc.field(SPURIOUS_FLAG) == 0)
)
```

Template fluxes are expected to be more accurate than aperture fluxes for extended sources because the templates do a better job of excluding contamination from nearby sources.
Conversely, aperture fluxes were found to be more accurate for point-like objects (esp. bright stars) in the Q1 NIR photometry, likely due to better handling of PSF related issues.
In either case, Euclid recommends scaling the measured fluxes in the NIR bands with a color term in order to obtain the best estimate of the total flux in that band, which we demonstrate here.

```{code-cell}
# Columns we actually want to load. Dict instead of list because we're defining new columns (magnitudes).
# We'll have PyArrow return only magnitudes so that we don't have to handle all the flux columns in memory ourselves.
_mag_columns = {
    # FLUX_TOTAL is the best estimate for the total flux in the detection band (here, VIS) and comes from
    # aperture photometry. VIS provides the template for NIR bands. It has no unique templfit flux itself.
    "I total": flux_to_magnitude(FLUX_TOTAL),
    # Template-fit fluxes.
    "Y templfit total": flux_to_magnitude("MER_FLUX_Y_TEMPLFIT", (FLUX_TOTAL, "MER_FLUX_VIS_TO_Y_TEMPLFIT")),
    "J templfit total": flux_to_magnitude("MER_FLUX_J_TEMPLFIT", (FLUX_TOTAL, "MER_FLUX_VIS_TO_J_TEMPLFIT")),
    "H templfit total": flux_to_magnitude("MER_FLUX_H_TEMPLFIT", (FLUX_TOTAL, "MER_FLUX_VIS_TO_H_TEMPLFIT")),
    # Aperture fluxes.
    "Y aperture total": flux_to_magnitude("MER_FLUX_Y_2FWHM_APER", (FLUX_TOTAL, "MER_FLUX_VIS_2FWHM_APER")),
    "J aperture total": flux_to_magnitude("MER_FLUX_J_2FWHM_APER", (FLUX_TOTAL, "MER_FLUX_VIS_2FWHM_APER")),
    "H aperture total": flux_to_magnitude("MER_FLUX_H_2FWHM_APER", (FLUX_TOTAL, "MER_FLUX_VIS_2FWHM_APER")),
}
mag_columns = {**_mag_columns, PHZ_CLASS: pc.field(PHZ_CLASS), MUMAX_MINUS_MAG: pc.field(MUMAX_MINUS_MAG)}
```

Load data.

```{code-cell}
mags_df = dataset.to_table(columns=mag_columns, filter=mag_filter).to_pandas()
# 30s
```

Given Euclid's core science goals, we'll take the template fluxes as our baseline in this section.

Plot the magnitude distributions as a function of PHZ class.
Include multiply-classed objects and also separate point-likes, given what we learned above.

```{code-cell}
classes = {"Galaxy": (2, 3, 6, 7), "Star": (1, 3), "QSO": (4, 6)}  # , "Star and QSO": (5, 7)}
class_colors = ["tab:green", "tab:blue", "tab:orange"]  # , "tab:brown"]

bands = ["I total", "Y templfit total", "J templfit total", "H templfit total"]
mag_limits = (14, 28)  # Excluding all magnitudes outside this range.
hist_kwargs = dict(bins=20, range=mag_limits, histtype="step")

fig, axes = plt.subplots(3, 4, figsize=(18, 12), sharey="row", sharex=True)
for (class_name, class_ids), class_color in zip(classes.items(), class_colors):
    # Get the objects that are in this class only.
    class_df = mags_df.loc[mags_df[PHZ_CLASS] == class_ids[0]]
    # Plot histograms for each band. Galaxies on top row, then stars, then QSOs.
    axs = axes[0] if class_name == "Galaxy" else (axes[1] if class_name == "Star" else axes[2])
    for ax, band in zip(axs, bands):
        ax.hist(class_df[band], label=class_name, color=class_color, **hist_kwargs)

    # Get the objects that are in this class and possibly others.
    class_df = mags_df.loc[mags_df[PHZ_CLASS].isin(class_ids)]
    label = class_name + "+Galaxy" if class_name != "Galaxy" else "+any"
    # Of those objects, restrict to the ones that are point-like.
    classpt_df = class_df.loc[class_df[MUMAX_MINUS_MAG] < -2.5]
    # Plot histograms for both sets of objects.
    for ax, band in zip(axs, bands):
        ax.hist(class_df[band], color=class_color, label=label, linestyle=":", **hist_kwargs)
        ax.hist(classpt_df[band], color=class_color, linestyle="-.", label=label + " (point-like)", **hist_kwargs)

# Add axis labels, etc.
for ax in axes[:, 0]:
    ax.set_ylabel("Counts")
    ax.legend(framealpha=0.2, loc=2)
for axs, band in zip(axes.transpose(), bands):
    axs[0].set_title(band.split()[0])
    axs[-1].set_xlabel(f"{band} (mag)")
plt.tight_layout()
```

Euclid is tuned to detect galaxies for cosmology studies, so it's no surprise that there are many more galaxies than other objects in here.
We also notice that the star distributions are broader and peak at a brighter magnitude than the galaxy distributions, as we would expect.
showing more confusion when the sources are faint. # [FIXME] comment on the mixed classes and point-like.
The bottom panel shows very few point-like QSOs, reminding us that most QSO classifications are suspect.

Finally, the bottom row shows evidence of the different probability thresholds that were used for acceptance into each class.
When we add the objects that were also classified as galaxies to the star and QSO distributions

+++

Now, compare with the aperture fluxes.
Plot the difference between the reference (template) and aperture magnitudes.

```{code-cell}
fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharey=True, sharex=True)
bands = [("Y templfit total", "Y aperture total"), ("J templfit total", "J aperture total"), ("H templfit total", "H aperture total")]
hexbin_kwargs = dict(bins="log", extent=(*mag_limits, -5, 5), gridsize=50)
outlier_threshold = 0.25

for axs, (ref_band, sersic_band) in zip(axes.transpose(), bands):
    # band_df = mags_df.loc[(mags_df[ref_band] > mag_limits[0]) & (mags_df[ref_band] < mag_limits[1])]
    # Extended objects, top row.
    ax = axs[0]
    extended_df = mags_df.loc[mags_df[MUMAX_MINUS_MAG] >= -2.5]
    # mag_diff = extended_df[ref_band] - extended_df[sersic_band]
    # _not_faint_df = extended_df.loc[(extended_df["I total"] < 24.5)]
    # outlier_frac = len(_not_faint_df.loc[mag_diff.abs() > outlier_threshold]) / len(_not_faint_df)
    # cb = ax.hexbin(extended_df["I total"], mag_diff, label=f"Outliers: {outlier_frac:.3f}", **hexbin_kwargs)
    cb = ax.hexbin(extended_df["I total"], extended_df[ref_band] - extended_df[sersic_band], **hexbin_kwargs)
    plt.colorbar(cb)
    # ax.legend()
    ax.set_ylabel(f"{ref_band} - {sersic_band}")
    # Point-like objects, bottom row.
    ax = axs[1]
    pointlike_df = mags_df.loc[mags_df[MUMAX_MINUS_MAG] < -2.5]
    # mag_diff = pointlike_df[ref_band] - pointlike_df[sersic_band]
    # _not_faint_df = pointlike_df.loc[(pointlike_df["I total"] < 24.5)]
    # outlier_frac = len(_not_faint_df.loc[mag_diff.abs() > outlier_threshold]) / len(_not_faint_df)
    # cb = ax.hexbin(pointlike_df["I total"], mag_diff, label=f"Outliers: {outlier_frac:.3f}", **hexbin_kwargs)
    cb = ax.hexbin(pointlike_df["I total"], pointlike_df[ref_band] - pointlike_df[sersic_band], **hexbin_kwargs)
    plt.colorbar(cb)
    # ax.legend()
    ax.set_ylabel(f"{ref_band} - {sersic_band}")
# Add axis labels, etc.
for i, ax in enumerate(axes.flatten()):
    # ax.axhline(outlier_threshold, color="gray", linewidth=1)
    # ax.axhline(-outlier_threshold, color="gray", linewidth=1)
    ax.axhline(0, color="gray", linewidth=1)
    ax.axvline(23, color="k", linewidth=1)
    if i == 1:
        ax.set_title("Extended objects")
    if i == 4:
        ax.set_title("Point-like objects")
    if i > 2:
        ax.set_xlabel("I total")
plt.tight_layout()
```

For point-like objects, # [TODO] comment.

+++

## Galaxy morphology

+++

In addition to the two point-like

---- notes ----
morphology correlated to different large-scale structure environments.
early = elliptical and lenticular, smooth. late = spirals and irregulars.
color bimodality is related to morphology. red ~ early and blud ~ late
sersic index n measures curvature of light profile.
  - bigger and brighter ~ larger n.
  - smaller n = shallow at small radius and steeper at large radius.
  - n = 0.5 == gaussian. 1 == exponential. 4 == de Vaucouleurs'
bulges and ellipticals, n > 2
massive ellipticals have high n, high concentration. spiral disks have low n, low concentration.

+++

Reproduce figures 13 and 5 from https://arxiv.org/pdf/2503.15309.

```{code-cell}
# PHZ_CLASS defined above.
PHZ_GAL_PROB = "CLASS_PHZ_GAL_PROB"

# Morphology cols
morph_nonpara_cols = ["MORPH_CONCENTRATION", "MORPH_SMOOTHNESS"]
sersic_vis_index = "MORPH_SERSIC_SERSIC_VIS_INDEX"
sersic_vis_radius = "MORPH_SERSIC_SERSIC_VIS_RADIUS"
sersic_vis_q = "MORPH_SERSIC_SERSIC_VIS_AXIS_RATIO"
morph_sersic_cols = [
    "MORPH_SERSIC_SERSIC_VIS_INDEX",
    "MORPH_SERSIC_SERSIC_VIS_RADIUS",
    "MORPH_SERSIC_SERSIC_VIS_AXIS_RATIO",
]
morph_zoobot_cols = [
    "MORPH_SMOOTH_OR_FEATURED_FEATURED_OR_DISK",  # *
    "MORPH_SMOOTH_OR_FEATURED_ARTIFACT_STAR_ZOOM",
    "MORPH_SMOOTH_OR_FEATURED_SMOOTH",
    "MORPH_DISK_EDGE_ON_YES",
    "MORPH_DISK_EDGE_ON_NO",
    "MORPH_BAR_STRONG",
    "MORPH_BAR_WEAK",
    "MORPH_BAR_NO",
    "MORPH_BULGE_SIZE_DOMINANT",
    "MORPH_BULGE_SIZE_LARGE",
    "MORPH_BULGE_SIZE_MODERATE",
    "MORPH_BULGE_SIZE_SMALL",
    "MORPH_BULGE_SIZE_NONE",
]
# Phys param cols
pp_cols = [
    PHZ_GAL_PROB,
    "PHYSPARAM_PHZ_PP_MEDIAN_STELLARMASS",
    "PHYSPARAM_PHZ_PP_MEDIAN_STELLARMETALLICITY",  # [Fe/H]
    "PHYSPARAM_PHZ_PP_MEDIAN_SFR",  # log10(SFR/Msun/yr)
    "PHYSPARAM_PHZ_PP_MEDIAN_SFHAGE",  # (SFH) stellar age of galaxy
    "PHYSPARAM_PHYS_PARAM_FLAGS",
]
mag_and_z_cols = [
    FLUX_TOTAL,
    PHYSPARAM_GAL_Z,  # from SED fitting. use for plotting.
    "PHYSPARAM_PHZ_PP_MODE_REDSHIFT",
    PHZ_Z,  # for cosmology. compare with above for quality cuts.
    "PHYSPARAM_PHZ_PP_MEDIAN_MU",
    "PHYSPARAM_PHZ_PP_MEDIAN_MR",
]
extra_query_cols = [
    FLUX_TOTAL,
    "PHZ_PHZ_CLASSIFICATION",
    VIS_DET,
    SPURIOUS_FLAG,
    POINTLIKE_PROB,
]
columns = morph_nonpara_cols + morph_sersic_cols + morph_zoobot_cols + pp_cols + mag_and_z_cols + extra_query_cols
```

```{code-cell}
# Cuts described in sec. 2. Skipping EDF cut for now.
vis_mag_max = 23  # 24.5 but restrain to <23 (recommended for reliable single-Sérsic fits)
vis_flux_min = magnitude_to_flux(vis_mag_max)
paquery = (
    (pc.field(FLUX_TOTAL) > vis_flux_min)
    & (pc.field("PHZ_PHZ_CLASSIFICATION") == 2)
    & (pc.field(VIS_DET) == 1)
    & (pc.field(SPURIOUS_FLAG) == 0)
    & (pc.field("MER_POINT_LIKE_PROB") <= 0.1)
    # limit of param space was 5.5 and there is an artificial peak there.
    # n > 5.45 should be removed from any Sérsic-based analysis.
    & (pc.field("MORPH_SERSIC_SERSIC_VIS_INDEX") <= 5.45)
    # [FIXME] Cut galaxies with sSFR < 10^-8.2/yr, per https://arxiv.org/pdf/2503.15306.
)
galaxies = dataset.to_table(columns=columns, filter=paquery).to_pandas()
# 2 min 21 sec medium server 4cpu. 24 columns
```

```{code-cell}
# Fig. 5
# Following [MER Morphology Cookbook - Zoobot](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/mermorphologycookbook.html#zoobot-morphology)
fig, axes = plt.subplots(2, 2, figsize=(16, 8))
ax1, ax2, ax3, ax4 = axes.flatten()

# ax1
zoo_smooth_cols = [f"MORPH_SMOOTH_OR_FEATURED_{feature}" for feature in ["FEATURED_OR_DISK", "SMOOTH", "ARTIFACT_STAR_ZOOM"]]
galaxies["zoo_smooth_total"] = galaxies[zoo_smooth_cols].sum(axis=1)
gals = galaxies.loc[galaxies.zoo_smooth_total > 0]
zoo_feat_or_disk = gals[zoo_smooth_cols[0]] / gals.zoo_smooth_total
zoo_smooth = gals[zoo_smooth_cols[1]] / gals.zoo_smooth_total

# Disk galaxies typically have ~exponential profiles (n=1)
# Elliptical typically have high n (n=4; de Vaucouleurs)
labels = ["featured_or_disk > 0.5", "featured_or_disk > 0.9", "smooth > 0.5", "smooth > 0.9"]
colors = ["tab:cyan", "tab:blue", "tab:pink", "tab:red"]
data = [
    gals.loc[zoo_feat_or_disk > 0.5, sersic_vis_index],
    gals.loc[zoo_feat_or_disk > 0.9, sersic_vis_index],
    gals.loc[zoo_smooth > 0.5, sersic_vis_index],
    gals.loc[zoo_smooth > 0.9, sersic_vis_index],
]
for sersic_indexes, label, color in zip(data, labels, colors):
    weights = [1 / len(sersic_indexes)] * len(sersic_indexes)
    ax1.hist(sersic_indexes, label=label, color=color, histtype="step", bins=20, weights=weights)
ax1.legend()
ax1.set_xlabel("Sérsic index")
ax1.set_ylabel("Fraction of galaxies")

# ax2
edgeon_cols = ["MORPH_DISK_EDGE_ON_NO", "MORPH_DISK_EDGE_ON_YES"]
edgeon_total = galaxies[edgeon_cols].sum(axis=1)
edgeon_no = galaxies[edgeon_cols[0]] / edgeon_total
edgeon_yes = galaxies[edgeon_cols[1]] / edgeon_total

data = [
    galaxies.loc[edgeon_no > 0.5, sersic_vis_q],
    galaxies.loc[edgeon_no > 0.9, sersic_vis_q],
    galaxies.loc[edgeon_yes > 0.5, sersic_vis_q],
    galaxies.loc[edgeon_yes > 0.9, sersic_vis_q],
]
labels = ["edgeon no > 0.5", "edgeon no > 0.9", "edgeon yes > 0.5", "edgeon yes > 0.9"]
colors = ["tab:gray", "tab:cyan", "tab:blue", "tab:pink", "tab:red"]
for sersic_q, label, color in zip(data, labels, colors):
    weights = [1 / len(sersic_q)] * len(sersic_q)
    ax2.hist(sersic_q, label=label, color=color, histtype="step", bins=20, weights=weights)
ax2.legend()
ax2.set_xlabel("Sérsic q")
ax2.set_ylabel("Fraction of galaxies")

# ax3 (Fig. 4)
galaxies.plot.hexbin(sersic_vis_index, "MORPH_CONCENTRATION", ax=ax3, bins="log")
ax3.set_xlabel("Sérsic index")
ax3.set_ylabel("Concentration")
```

## NIR-only detections: high-redshift object or nearby ultra-cool star?

+++

Suppose we are looking for high-z objects which are far enough away that their light has been redshifted out of the VIS band.
https://arxiv.org/pdf/2503.15306 Tells us that >20% of MER objects are NIR-only detections (MER_VIS_DET=0).
  - says NIR-only => very likely either high-z & very luminous, brown dwarf, or else spurious (Andreas adds, "or very dust-obscured").
Given such a large fraction of NIR-only objects, suppose we want extract a sample that is relatively free of contamination from brown dwarfs and spurious detections prior to looking through the spectra.

```{code-cell}
OBJECT_ID = "OBJECT_ID"
PHZ_GAL_PROB = "CLASS_PHZ_GAL_PROB"
PHZ_QSO_PROB = "CLASS_PHZ_QSO_PROB"
PHZ_STAR_PROB = "CLASS_PHZ_STAR_PROB"
SPE_CLASS = "Z_SPE_CLASS"
SPE_STAR_PROB = "Z_SPE_STAR_PROB"
SPE_GAL_PROB = "Z_SPE_GAL_PROB"
SPE_QSO_PROB = "Z_SPE_QSO_PROB"
PHYSPARAMNIR_Z = "PHYSPARAMNIR_Z"
PHYSPARAMNIR_HIGH_Z_PROB = "PHYSPARAMNIR_HIGH_Z_PROB"  # want >0.8?
PHZ_BEST_CHI2 = "PHZ_BEST_CHI2"
PHYSPARAMNIR_Z_1D_PDF = "PHYSPARAMNIR_Z_1D_PDF"
PHYSPARAMNIR_SED = "PHYSPARAMNIR_SED"
PHZ_PHZ_PDF = "PHZ_PHZ_PDF"
PHZ_Z = "PHZ_PHZ_MEDIAN"
```

```{code-cell}
# Query for all nir-only. Plot zphos vs pz6, colored by chi2.
columns = [
    OBJECT_ID,
    PHYSPARAMNIR_Z,
    PHYSPARAMNIR_HIGH_Z_PROB,
    PHZ_BEST_CHI2,
    PHZ_Z,
    PHZ_CLASS,
    PHZ_GAL_PROB,
    PHZ_QSO_PROB,
    PHZ_STAR_PROB,
    SPE_CLASS,
    SPE_STAR_PROB,
    SPE_GAL_PROB,
    SPE_QSO_PROB,
    SPURIOUS_FLAG,
    MUMAX_MINUS_MAG,
    "MER_RIGHT_ASCENSION",
    "MER_DECLINATION",
]

# query = "MER_VIS_DET == 0 & PHYSPARAMNIR_Z >= 0 & PHYSPARAMNIR_HIGH_Z_PROB >= 0"
# nironly = lsdb.read_hats(euclid_hats_collection_uri, columns=columns).query(query).compute()
nironly = dataset.to_table(columns=columns, filter=(pc.field(VIS_DET) == 0) & (pc.field(SPURIOUS_FLAG) == 0)).to_pandas()
```

```{code-cell}
brown_dwarfs = [-523574860290315045, -600367386508373277]
highz_galaxies = [-531067351279302418]  # -531510650277773811, -531639279277732979

# nironly = nironly.dropna(subset=[physparam_nir_z, PHYSPARAMNIR_HIGH_Z_PROB], how="any")
fig, axes = plt.subplots(1, 2, figsize=(18, 6))
x, y = PHYSPARAMNIR_Z, PHYSPARAMNIR_HIGH_Z_PROB
class_groups = {"Stars": (1, 3, 5, 7), "Galaxies and QSOs": (2, 4, 3, 5, 6, 7)}
for ax, (class_name, class_ids) in zip(axes, class_groups.items()):
    nironly.loc[nironly[PHZ_CLASS].isin(class_ids)].plot.hexbin(x, y, bins="log", ax=ax, cmap="Greys")
    nironly.loc[nironly.OBJECT_ID.isin(brown_dwarfs)].plot.scatter(x, y, ax=ax, c="tab:cyan", s=PHZ_BEST_CHI2, label="Brown Dwarf")
    nironly.loc[nironly.OBJECT_ID.isin(highz_galaxies)].plot.scatter(x, y, ax=ax, c="tab:orange", s=PHZ_BEST_CHI2, label="High-z Galaxy")
    ax.scatter(6, 0.8, marker=(3, 0, -45), s=250, edgecolors="k", facecolors="none")
    ax.scatter(0, 0.8, marker="^", s=250, edgecolors="k", facecolors="none")
    ax.set_title(class_name)
#       - zphos > 6 and pz6 > 0.8 => good high-z candidate. but, disfavor if chi2 >> 10.
#       - zphos > 6 and low-mid pz6 => likely broad or multi-peaked zPDF
#       - zphos = 0 and large pz6 => model degeneracy between high-z galaxies and brown dwarfs
```

[TODO] comment

+++

Now look at PDFs for our three objects.

```{code-cell}
columns = [
    OBJECT_ID,
    PHYSPARAMNIR_Z_1D_PDF,
    PHYSPARAMNIR_SED,
    PHYSPARAMNIR_Z,
    PHYSPARAMNIR_HIGH_Z_PROB,
    PHZ_BEST_CHI2,
    PHZ_Z,
    PHZ_PHZ_PDF,
    "Z_GALAXY_CANDIDATES_SPE_PDF_RANK0",
    "Z_GALAXY_CANDIDATES_SPE_PDF_ZMIN_RANK0",
    "Z_GALAXY_CANDIDATES_SPE_PDF_ZMAX_RANK0",
    "Z_GALAXY_CANDIDATES_SPE_PDF_DELTAZ_RANK0",
    "Z_QSO_CANDIDATES_SPE_PDF_RANK0",
    "Z_QSO_CANDIDATES_SPE_PDF_ZMIN_RANK0",
    "Z_QSO_CANDIDATES_SPE_PDF_ZMAX_RANK0",
    "Z_QSO_CANDIDATES_SPE_PDF_DELTAZ_RANK0",
    PHZ_CLASS,
    PHZ_GAL_PROB,
    PHZ_QSO_PROB,
    PHZ_STAR_PROB,
    SPE_CLASS,
    SPE_STAR_PROB,
    SPE_GAL_PROB,
    SPE_QSO_PROB,
]
# collection = lsdb.read_hats(euclid_hats_collection_uri, columns=columns)
# index_df = pd.concat([collection.id_search({"OBJECT_ID": srcid}).compute() for srcid in brown_dwarfs + highz_galaxies])  # 7m 40s
index_df = dataset.to_table(columns=columns, filter=pc.field(OBJECT_ID).isin(brown_dwarfs + highz_galaxies)).to_pandas()  # 1m 7s
```

```{code-cell}
# objid = -523574860290315045  # from https://arxiv.org/pdf/2503.22442
# "PHYSPARAMNIR_Z_1D_PDF"  # 2503.15306: z<6 step is 0.01; z 6-10 step is 0.05
nir_zbins = np.concatenate([np.arange(0, 6, 0.01), np.arange(6, 10.01, 0.05)])
step_kwargs = dict(where="mid", alpha=0.8)

fig, axes = plt.subplots(1, 3, figsize=(24, 6), sharey=False)
for ax, (objid) in zip(axes, brown_dwarfs + highz_galaxies):
    obj_df = index_df.loc[index_df.OBJECT_ID == objid]
    ax.step(nir_zbins[:601], obj_df[PHZ_PHZ_PDF].squeeze(), label="PHZ zPDF", **step_kwargs)
    ax.step(nir_zbins, obj_df[PHYSPARAMNIR_Z_1D_PDF].squeeze(), label="PHYSPARAMNIR zPDF", **step_kwargs)
    # "Z_QSO_CANDIDATES_SPE_PDF_RANK0"  # how to plot? min, max, delta, and len don't make sense.
    ax.legend()
    ax.set_xlabel("Redshift")
    ax.set_title(objid)
```

## 4. Schema details

This data product contains the Euclid Q1 tables listed below, joined on 'OBJECT_ID' into a single dataset.
In addition, the Euclid 'TILEID' for each object has been added, along with a few HATS- and HEALPix-related columns.
All column names other than 'OBJECT_ID' and 'TILEID' have the table name prepended (e.g., 'DECLINATION' -> 'MER_DECLINATION').
In addition, all non-alphanumeric characters have been replaced with an underscore for compatibility with various libraries and services (e.g., 'E(B-V)' -> 'PHYSPARAMQSO_E_B_V_').
Finally, three of the spec-z tables required special handling and more information is given below.

- [MER](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#main-catalog-fits-file)
- [CUTOUTS](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#cutouts-catalog-fits-file)
- [MORPH](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#morphology-catalog-fits-file)
- [GALAXYSED](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#galaxy-sed-catalog)
- [PHZ](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#photo-z-catalog)
- [STARSED](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#star-sed-catalog)
- [CLASS](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#classification-catalog)
- [PHYSPARAM](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#physical-parameters-catalog)
- [PHYSPARAMNIR](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#nir-physical-parameters-catalog)
- [PHYSPARAMQSO](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#qso-physical-parameters-catalog)
- [STARCLASS](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#star-template)
- [Z](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/spedpd/dpcards/spe_spepfoutputcatalog.html#redshift-catalog)
    - The original FITS files for Z contain the spectroscopic redshift estimates for GALAXY_CANDIDATES, STAR_CANDIDATES, and QSO_CANDIDATES (in HDUs 3, 4, and 5 respectively) which required special handling to be included in this Parquet product.
      There are up to 5 redshift estimates per 'OBJECT_ID' (and per HDU).
      For the Parquet, the tables were pivoted so that there is one row per 'OBJECT_ID' to facilitate the table joins.
      The resulting columns are named by combining the table name (Z), the HDU name, the original column name, and the rank of the given redshift estimate (i.e., the value in the original 'SPE_RANK' column).
      For example, the 'SPE_PDF' column for the highest ranked redshift estimate in the 'GALAXY_CANDIDATES' table is called 'Z_GALAXY_CANDIDATES_SPE_PDF_RANK0'.

Here, we follow IRSA's
[Cloud Access notebook](https://caltech-ipac.github.io/irsa-tutorials/tutorials/cloud_access/cloud-access-intro.html#navigate-a-catalog-and-perform-a-basic-query)
to inspect the parquet schema.

```{code-cell}
s3_filesystem = pyarrow.fs.S3FileSystem()
schema = pyarrow.parquet.read_schema(euclid_parquet_schema_path, filesystem=s3_filesystem)
```

There are more than 1300 columns in this dataset.

```{code-cell}
print(f"{len(schema)} columns")
```

### 4.1 Euclid Q1 columns

+++

To find all columns from a given table, search for column names that start with the table name followed by an underscore.

```{code-cell}
# Find all column names from the photo-z table (PHZ).
[name for name in schema.names if name.startswith("PHZ_")]
```

Euclid Q1 offers many flux measurements, both from Euclid detections and from external ground-based observatories (see https://arxiv.org/pdf/2503.15305).
The latter are identified by the inclusion of "EXT" in the name, such as "MER_FLUX_U_EXT_DECAM_1FWHM_APER"
Even restricting to Euclid fluxes, there are 327 related columns, so it's helpful to narrow down further.
We'll look in the main MER table.

```{code-cell}
# Find all Euclid MER flux columns.
[name for name in schema.names if name.startswith("MER_") and "FLUX" in name and "EXT" not in name]
```

Column metadata includes unit and description.

```{code-cell}
schema.field("MER_FLUX_Y_4FWHM_APER").metadata
```

### 4.2 Additional columns

+++

In addition to the columns from Euclid Q1 tables, the following columns have been added to this dataset:

- 'TILEID' : Euclid MER tile index.
- '_healpix_19' : HEALPix order 19 pixel index. Useful for spatial queries.
- '_healpix_29' : (hats column) HEALPix order 29 pixel index. Useful for spatial queries.
- 'Norder' : (hats column) HEALPix order at which the data is partitioned.
- 'Npix' : (hats column) HEALPix pixel index at order Norder.
- 'Dir' : (hats column) Integer equal to 10_000 * floor[Npix / 10_000].

These columns are at either the beginning or the end of the schema.

```{code-cell}
schema.names[:5]
```

```{code-cell}
schema.names[-5:]
```

## About this notebook

**Authors:** Troy Raen (Developer; Caltech/IPAC-IRSA) and the IRSA Data Science Team.

**Updated:** 2025-06-04

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.
