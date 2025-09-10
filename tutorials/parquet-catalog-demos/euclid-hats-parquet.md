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

This notebook introduces the [Euclid Q1](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_q1.html) HATS Collection served by IPAC/IRSA and demonstrates access with python.

+++

## Learning Goals

By the end of this tutorial, you will:

- Query the dataset to find and create figures for galaxies, QSOs, and stars with quality fluxes, redshifts, and morphology.
- Understand the format and schema of this dataset.
- Learn how to work with this HATS Parquet product using the LSDB and PyArrow python libraries.

+++

## 1. Introduction

The Collection includes a HATS Catalog (main data product), Margin Cache (10 arcsec), and Index Table (object_id).
The Catalog includes the twelve Euclid Q1 tables listed below, joined on the column 'object_id' into a single Parquet dataset with 1,594 columns.
There are 29,953,430 rows, one per Euclid MER Object, and the total data volume is 400 GB.
The data includes several different redshift measurements, several flux measurements for each Euclid band, and flux measurements for bands from several ground-based observatories -- in addition to morphological and other measurements.
Each was produced for different science goals using different algorithms and/or configurations.

Having all columns in the same dataset makes access convenient because the user doesn't have to make separate calls for data from different tables and/or join the results.
However, figuring out which columns to use amongst so many can be challenging.
In the sections below, we look at some of their distributions and reproduce figures from several papers in order to highlight some of the options and point out their differences.
The Appendix explains how the columns in this Parquet dataset are named and organized.
For more information about the meaning and provenance of a column, refer to the links provided with the list of tables below.

### 1.1 Euclid Q1 tables and docs

The Euclid Q1 HATS Catalog includes the following 14 Q1 tables which are organized underneath the Euclid processing function (MER, PHZ, or SPE) that created it.
Links to the Euclid papers describing the processing functions are provided, as well as pointers for each table.
Table names are linked to their original schemas.

- MER - [Euclid Collaboration: Romelli et al., 2025](https://arxiv.org/pdf/2503.15305) (hereafter, Romelli)
  - [mer](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#main-catalog-fits-file) - Sec. 6 & 8 (EUC_MER_FINAL-CAT)
  - [morph](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#morphology-catalog-fits-file) - Sec. 7 & 8 (EUC_MER_FINAL-MORPH-CAT)
  - [cutouts](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#cutouts-catalog-fits-file) - Sec. 8 (EUC_MER_FINAL-CUTOUTS-CAT)
- PHZ - [Euclid Collaboration: Tucci et al., 2025](https://arxiv.org/pdf/2503.15306) (hereafter, Tucci)
  - [phz](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#photo-z-catalog) - Sec. 5 (phz_photo_z)
  - [class](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#classification-catalog) - Sec. 4 (phz_classification)
  - [physparam](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#physical-parameters-catalog) - Sec. 6 (6.1; phz_physical_parameters) _Notice this is for **galaxies**._
  - [galaxysed](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#galaxy-sed-catalog) - App. B (B.1 phz_galaxy_sed)
  - [physparamqso](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#qso-physical-parameters-catalog) - Sec. 6 (6.2; phz_qso_physical_parameters)
  - [starclass](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#star-template) - Sec. 6 (6.3; phz_star_template)
  - [starsed](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#star-sed-catalog) - App. B (B.1 phz_star_sed)
  - [physparamnir](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#nir-physical-parameters-catalog) - Sec. 6 (6.4; phz_nir_physical_parameters)
- SPE - [Euclid Collaboration: Le Brun et al., 2025](https://arxiv.org/pdf/2503.15308) (hereafter, Le Brun)
  - [z](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/spedpd/dpcards/spe_spepfoutputcatalog.html#redshift-catalog) - Sec. 2 (spectro_zcatalog_spe_quality, spectro_zcatalog_spe_classification, spectro_zcatalog_spe_galaxy_candidates, spectro_zcatalog_spe_star_candidates, and spectro_zcatalog_spe_qso_candidates)
  - [lines](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/spedpd/dpcards/spe_spepfoutputcatalog.html#lines-catalog) line measurements (HDU1) rows with SPE_LINE_NAME == Halpha only - Sec. 5 (spectro_line_features_catalog_spe_line_features_cat) _Notice that lines were identified assuming **galaxy** regardless of the classification._
  - [models](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/spedpd/dpcards/spe_spepfoutputcatalog.html#models-catalog) galaxy candidates (HDU2) only - Sec. 5 (spectro_model_catalog_spe_lines_catalog)

See also:

- MER [Photometry](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/merphotometrycookbook.html) and [Morphology](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/mermorphologycookbook.html) Cookbooks
- [Frequently Asked Questions About Euclid Q1 data](https://euclid.caltech.edu/page/euclid-q1-data-faq) (hereafter, FAQ)
- [Q1 Explanatory Supplement](https://euclid.esac.esa.int/dr/q1/expsup/)

### 1.2 Parquet, HEALPix, and HATS

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

## 2. Installs and imports

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
# import os  # Determine number of CPUs (for parallelization)
# import dask.distributed  # Parallelize catalog queries
# import lsdb  # Query the catalog
# import matplotlib.colors  # Make figures look nice
import hpgeom
import matplotlib.pyplot as plt  # Create figures
import numpy as np  # Math
import pandas as pd  # Manipulate query results
import pyarrow.compute as pc  # Filter dataset
import pyarrow.dataset  # Load the dataset
import pyarrow.parquet  # Load the schema
import pyarrow.fs  # Simple S3 filesystem pointer

# Copy-on-write will become the default in pandas 3.0 and is generally more performant
pd.options.mode.copy_on_write = True
```

```{tip}
If you run into an error that looks like,

> AttributeError: _ARRAY_API not found

or:

> A module that was compiled using NumPy 1.x cannot be run in NumPy 2.1.3 as it may crash.

make sure you have restarted the kernel since doing `pip install`. Then re-run the cell **twice**.
```

+++

## 3. Setup

+++

### 3.1 AWS S3 paths

```{code-cell}
s3_bucket = "nasa-irsa-euclid-q1"
euclid_prefix = "contributed/q1/merged_objects/hats"

euclid_hats_collection_uri = f"s3://{s3_bucket}/{euclid_prefix}"  # for lsdb
euclid_parquet_metadata_path = f"{s3_bucket}/{euclid_prefix}/euclid_q1_merged_objects-hats/dataset/_metadata"  # for pyarrow
euclid_parquet_schema_path = f"{s3_bucket}/{euclid_prefix}/euclid_q1_merged_objects-hats/dataset/_common_metadata"  # for pyarrow

s3_filesystem = pyarrow.fs.S3FileSystem(anonymous=True)
```

### 3.2 Helper functions

```{code-cell}
def magnitude_to_flux(magnitude: float) -> float:
    """Convert magnitude to flux following the MER Photometry Cookbook."""
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
    If `color_col_names` is None, `flux_col_name` is taken as the total flux in the band.
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

### 3.3 Load the catalog as a PyArrow dataset

```{code-cell}
# Load the catalog as a PyArrow dataset. This is used in many examples below.
dataset = pyarrow.dataset.parquet_dataset(euclid_parquet_metadata_path, partitioning="hive", filesystem=s3_filesystem)
```

### 3.4 Frequently used columns

+++

The following columns will be used throughout this notebook.
Many other columns are defined in the sections below where they are used.
Descriptors generally come from the respective paper (Romelli, Tucci, or Le Brun) unless noted.

```{code-cell}
# Object ID set by the MER pipeline.
OBJECT_ID = "object_id"
```

Flux and source detection columns.

```{code-cell}
# Whether the source was detected in the VIS mosaic (1) or only in the NIR-stack mosaic (0).
VIS_DET = "mer_vis_det"

# Best estimate of the total flux in the detection band. From aperture photometry within a Kron radius.
# Detection band is VIS if mer_vis_det=1.
# Otherwise, this is a non-physical NIR-stack flux and there was no VIS detection (aka, NIR-only).
FLUX_TOTAL = "mer_flux_detection_total"
FLUXERR_TOTAL = "mer_fluxerr_detection_total"
```

Point-like and spurious indicators.

```{code-cell}
# Peak surface brightness minus the magnitude used for mer_point_like_prob.
# Point-like: <-2.5. Compact: <-2.6. (Tucci)
MUMAX_MINUS_MAG = "mer_mumax_minus_mag"

# Probability from the star-galaxy classifier. Heavily biased toward high purity.
# This is always NaN for NIR-only objects (use mer_mumax_minus_mag instead).
POINTLIKE_PROB = "mer_point_like_prob"

# Whether the detection has a >50% probability of being spurious (1=Yes, 0=No).
SPURIOUS_FLAG = "mer_spurious_flag"
```

PHZ classifications. These were generated by a probabilistic random forest supervised ML algorithm.

```{code-cell}
PHZ_CLASS = "phz_phz_classification"
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

### 3.5 Euclid Deep Fields

+++

[FIXME] The notebook does not currently use these. Should either use them or remove them.

Euclid Q1 includes data from three Euclid Deep Fields: EDF-N (North), EDF-S (South), EDF-F (Fornax; also in the southern hemisphere).
There is also a small amount of data from a fourth field: LDN1641 (Lynds' Dark Nebula 1641), which was observed for technical reasons during Euclid's verification phase and mostly ignored here.
The regions are well separated, so we can distinguish them using a simple cone search without having to be too picky about the radius.
Rather than using the RA and Dec values directly, we'll find a set of HEALPix order 9 pixels that cover each area.
A column ('_healpix_9') of order 9 indexes was added to the catalog for this purpose.
These will suffice for a simple and efficient cone search.

```{code-cell}
# Column name of HEALPix order 9 pixel indexes.
HEALPIX_9 = "_healpix_9"

# EDF-N (Euclid Deep Field - North)
ra, dec, radius = 269.733, 66.018, 4  # need ~20 sq deg
edfn_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius)

# EDF-S (Euclid Deep Field - South)
ra, dec, radius = 61.241, -48.423, 5  # need ~23 sq deg
edfs_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius)

# EDF-F (Euclid Deep Field - Fornax)
ra, dec, radius = 52.932, -28.088, 3  # need ~10 sq deg
edff_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius)

# LDN1641 (Lynds' Dark Nebula 1641)
ra, dec, radius = 85.74, -8.39, 1.5  # need ~6 sq deg
ldn_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius)
```

## 4. Redshifts for cosmology

+++

Euclid is a cosmology mission focused on measuring the evolution of large-scale structures in order to study dark matter and dark energy.
This means it must determine the redshifts for a large number of galaxies.
In this section, we obtain and examine quality cosmological samples for three of the redshift point-estimates (PDFs deferred to sec. 6) provided in Q1:

```{code-cell}
PHZ_Z = "phz_phz_median"
PHYSPARAM_GAL_Z = "physparam_phz_pp_median_redshift"
SPE_GAL_Z = "z_galaxy_candidates_spe_z_rank0"
```

"phz_phz_median" is the median of the photometric redshift PDF that was produced for Euclid's core-science goals.
It was generated by Phosphoros, a fully Bayesian template-fitting code.
It was computed for all MER objects, but the input models assumed galaxy.
The model grid was built spanning the parameters: redshift (z in [0, 6]), galaxy SED, intrinsic reddening curve, intrinsic attenuation.
Phosphoros should be better for cosmology than ML algorithms (which are more typical) due to the scarcity of spectroscopic "truth" training data above z ~ 1.

"physparam_phz_pp_median_redshift" is the median of the photometric redshift PDF that was produced for galaxy-classed objects by the physical-properties branch of the PHZ pipeline (i.e., non-cosmology).
We include it here as a useful comparison.
It was generated by NNPZ, a k-nearest neighbors supervised learning algorithm.
The reference sample came from >1 million stellar population synthesis models spanning parameters: redshift (z in [0, 7]), age, star formation timescale, stellar metallicity, intrinsic attenuation, magnitude, and two dust laws.

"z_galaxy_candidates_spe_z_rank0" is the top-ranked spectroscopic redshift estimate that was produced assuming the object is a galaxy (independent of the actual classification).
This is the most prominent peak of a PDF generated by a least-squares fitting algorithm (19 galaxy models).

Note that the redshift estimates are based in part on external data, which is different in EDF-N vs EDF-S and EDF-F.
We mostly ignore that here for simplicity, but additional cuts could be added using the HEALPix order 9 indexes found above to separate the results.

+++

Load a quality PHZ sample. Cuts are from Tucci sec. 5.3.

```{code-cell}
# Photo-z flag: 0=good for core science, 10=NIR-only, 11=missing bands, 12=too faint.
PHZ_FLAG = "phz_phz_flags"

# Columns we actually want to load.
phz_columns = [OBJECT_ID, PHZ_Z]

# Filter for quality PHZ redshifts.
phz_filter = (
    (pc.field(VIS_DET) == 1)  # No NIR-only objects.
    & (pc.field(FLUX_TOTAL) > magnitude_to_flux(24.5))  # I < 24.5
    & (pc.divide(pc.field(FLUX_TOTAL), pc.field(FLUXERR_TOTAL)) > 5)  # I band S/N > 5
    & ~pc.field(PHZ_CLASS).isin([1, 3, 5, 7])  # Exclude objects classified as star.
    & (pc.field(SPURIOUS_FLAG) == 0)  # MER quality
)

# Execute the filter and load.
phz_df = dataset.to_table(columns=phz_columns, filter=phz_filter).to_pandas()
phz_df = phz_df.set_index(OBJECT_ID).sort_index()
# 1m 6s
```

Load a quality PHYSPARAM sample. Cuts are from Tucci sec. 6.1.2.

```{code-cell}
# Properties to calculate specific star formation rate (sSFR).
PHYSPARAM_GAL_SFR = "physparam_phz_pp_median_sfr"  # log10(SFR [Msun/yr])
PHYSPARAM_GAL_MSTAR = "physparam_phz_pp_median_stellarmass"  # log10(Stellar Mass [Msun])

# Columns we actually want to load.
# We'll have pyarrow construct and return the sSFR, so we must pass a dict mapping column names to expressions.
log10_ssfr = pc.subtract(pc.field(PHYSPARAM_GAL_SFR), pc.field(PHYSPARAM_GAL_MSTAR))
pp_columns = {PHYSPARAM_GAL_Z: pc.field(PHYSPARAM_GAL_Z), "log10_ssfr": log10_ssfr, OBJECT_ID: pc.field(OBJECT_ID)}

# Partial filter for quality PHYSPARAM redshifts.
pp_galaxy_filter = (
    (pc.field(MUMAX_MINUS_MAG) > -2.6)  # Non-compact objects
    & (pc.field(PHZ_FLAG) == 0)  # Good for core science
    & (pc.field(SPURIOUS_FLAG) == 0)  # MER quality
    & (pc.field("mer_det_quality_flag") < 4)  # MER quality
    & ~pc.is_null(pc.field(PHYSPARAM_GAL_Z))  # Galaxy class and redshift solution found
)

# Execute the filter and load.
pp_df = dataset.to_table(columns=pp_columns, filter=pp_galaxy_filter).to_pandas()
pp_df = pp_df.set_index(OBJECT_ID).sort_index()
# 1m 18s

# Final filter, to be applied later.
# sSFR < 10^-8.2 /yr. This excludes galaxies with unrealistically young ages and very high sSFR.
pp_final_filter = pp_df["log10_ssfr"] < -8.2
```

Load a quality SPE sample. Cuts are from Le Brun sec. 3.3 and 6.2.

The NISP instrument was built to target Halpha emitting galaxies, which effectively means 0.9 < z < 1.8.
SPE redshifts are reliable in that regime.
However, this represents <2% of the total delivered by the SPE pipeline, so it's crucial to make cuts in order to get it.

```{code-cell}
# SPE probability of the rank 0 (best) redshift estimate, assuming galaxy.
SPE_GAL_Z_PROB = "z_galaxy_candidates_spe_z_prob_rank0"
# [TODO] describe
HALPHA_LINE_FLUX = "lines_spe_line_flux_gf_rank0_halpha"
HALPHA_LINE_SNR = "lines_spe_line_snr_gf_rank0_halpha"
EMISSION_LINE_WIDTH = "models_galaxy_spe_vel_disp_e_rank0"

# Columns we actually want to load.
spe_columns = [SPE_GAL_Z, OBJECT_ID]

# Filter for quality SPE galaxy redshifts.
spe_filter = (
    # Euclid's target redshift range.
    (pc.field(SPE_GAL_Z) > 0.9)
    & (pc.field(SPE_GAL_Z) < 1.8)
    # MER quality
    & (pc.field(SPURIOUS_FLAG) == 0)
    & (pc.field("mer_det_quality_flag") < 4)
    # High quality SPE galaxies.
    & (pc.field(SPE_GAL_Z_PROB) > 0.99)  # [FIXME] Andreas says > 0.999. Also mentioned in sec. 6.2.
    & (pc.field(EMISSION_LINE_WIDTH) < 680)
    # Halpha emitters.
    & (pc.field(HALPHA_LINE_FLUX) > 2e-16)
    & (pc.field(HALPHA_LINE_SNR) > 3.5)  # [FIXME] Tiffany's notebook uses > 5.
    # These make no difference
    # & (pc.field("LINES_SPE_LINE_N_DITH_RANK0_Halpha") >= 3)  # all values in this column = 0
    # & (pc.field("Z_SPE_N_DITH_MED") >= 3)
    # & (pc.field("Z_SPE_ERROR_FLAG") == 0)
    # & (pc.field("Z_SPE_GAL_ERROR_FLAG") == 0)
)

# Execute the filter and load.
spe_df = dataset.to_table(columns=spe_columns, filter=spe_filter).to_pandas()
spe_df = spe_df.set_index(OBJECT_ID).sort_index()
# 1m 2s
```

Plot redshift distributions

```{code-cell}
tbl_colors = {"PHZ": "tab:orange", "PHYSPARAM": "tab:green", "SPE_GAL": "tab:purple"}

fig, ax = plt.subplots(1, 1, figsize=(12, 9))
hist_kwargs = dict(bins=np.linspace(0, 7.1, 100), histtype="step", linewidth=2, alpha=0.7)

# PHZ
phz_kwargs = dict(label=PHZ_Z + " (quality)", color=tbl_colors["PHZ"], linestyle="-", zorder=9)
ax.hist(phz_df[PHZ_Z], **phz_kwargs, **hist_kwargs)

# PHYSPARAM
hist_kwargs.update(linewidth=1)
pp_kwargs = dict(label=PHYSPARAM_GAL_Z + " (filtered)", color=tbl_colors["PHYSPARAM"], linestyle=":")
ax.hist(pp_df[PHYSPARAM_GAL_Z], **pp_kwargs, **hist_kwargs)
# Impose our final cuts.
pp_kwargs.update(label=PHYSPARAM_GAL_Z + " (quality)", linestyle="-")
ax.hist(pp_df.loc[pp_final_filter, PHYSPARAM_GAL_Z], **pp_kwargs, **hist_kwargs)

# SPE
spe_kwargs = dict(label=SPE_GAL_Z + " (filtered)", color=tbl_colors["SPE_GAL"], linestyle=":")
ax.hist(spe_df[SPE_GAL_Z], **spe_kwargs, **hist_kwargs)

ax.set_xlabel("Redshift")
ax.set_ylabel("Counts")
plt.legend()
```

The orange distribution is a quality sample of the redshifts (best point estimate) generated for cosmology by a Bayesian template-fitting code.
The maximum is z ~ 6, due to the model's input parameters.
Green represents redshifts that were generated to study galaxies' physical properties by a supervised learning, k-nearest neighbors algorithm.
The maximum is z ~ 7, again due to model inputs.
Several quality cuts were applied to produce the dotted-line sample, but this still includes a population of problematic galaxies for which the solutions pointed to unrealistically young ages and very high specific star formation rates.
The green solid line filters those out and represents a quality sample for this redshift estimate.
Purple represents the spectroscopic redshifts (best point estimate).
The dotted line has been filtered for reliable (SPE) galaxy solutions and the maximum is z ~ 5.
There is a clear bump between about 0.9 < z < 1.8 which results from a combination of the NISP instrument parameters (tuned to detect Halpha) and a model prior that strongly favored solutions in this regime.
However much more drastic cuts are needed to obtain a trustworthy sample.
The purple solid line filters down to the redshift target range, within which we have the most confidence, but still includes interlopers which were either noisy or had a different spectral line misidentified as Halpha.
Further cuts will need to be made once the Halpha line information is added to the parquet files # [FIXME].

+++

Next, compare the redshift estimates, treating PHZ (cosmology photo-z) as our ground truth.

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
Here, we reproduce Tucci Fig. 17 (left panel) except that we don't consider the problematic galaxies nor do we impose cuts on magnitude or region (EDF-F).

```{code-cell}
# Get the common objects and set axes data (PHZ on x, PHYSPARAM on y).
phz_pp_df = phz_df.join(pp_df.loc[pp_final_filter], how="inner", lsuffix="phz", rsuffix="pp")
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

Above is the redshift point-estimate comparison of quality samples of PHZ (cosmology photo-z) vs PHYSPARAM (physical-parameters photo-z).
The two agree quite well even though they were generated by very different different algorithms.
The two outlier clouds are very roughly similar to those in Tucci Fig. 7 which were attributed to a confusion between the Balmer and Lyman breaks.

+++

Compare PHZ to SPE

```{code-cell}
# Get the common objects and set axes data (PHZ on x, SPE on y).
phz_spe_df = phz_df.join(spe_df, how="inner", lsuffix="phz", rsuffix="spe")
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

[FIXME] This should look better once the rest of the cuts (esp. for Halpha line) can be added. Then need to comment.

+++

## 5. Classification thresholds - purity vs completeness

+++

The PHZ_CLASS used above was determined using probability thresholds (different in north vs south due to availability of external data) for galaxies which prioritized completeness over purity, while the others (esp. for stars) gave purity more priority.
This was done in order meet certain requirements imposed by Euclid's main goals, such as the need for a very pure sample of high S/N stars in order to model the PSF.
One way to see the effects of these differing thresholds are to look at the brightness and morphology as a function of class.

```{code-cell}
# Construct the filter.
class_filter = (
    (pc.field(PHZ_CLASS).isin([1, 2, 4]))  # Stars, galaxies, and QSOs (no mixed classes or Undefined).
    & (pc.divide(pc.field(FLUX_TOTAL), pc.field(FLUXERR_TOTAL)) > 5)  # S/N > 5
    & (pc.field(VIS_DET) == 1)  # No NIR-only objects, so FLUX_TOTAL == VIS flux. For simplicity; not required.
)

# Columns we actually want to load.
# Because we want PyArrow to construct a new column (magnitude), we must pass a dict mapping column names to expressions.
class_columns = {"I magnitude": flux_to_magnitude(FLUX_TOTAL), MUMAX_MINUS_MAG: pc.field(MUMAX_MINUS_MAG), PHZ_CLASS: pc.field(PHZ_CLASS)}
```

```{code-cell}
# Load data.
classes_df = dataset.to_table(columns=class_columns, filter=class_filter).to_pandas()
```

Plot point-like morphology vs brightness as a function of class.
Here, we reproduce the first three panels of Tucci Fig. 6, combining top and bottom.

```{code-cell}
fig, axes = plt.subplots(1, 3, figsize=(20, 6))
for ax, (class_name, class_df) in zip(axes, classes_df.groupby(PHZ_CLASS)):
    ax.set_title(PHZ_CLASS_MAP[class_name])
    cb = ax.hexbin(class_df[MUMAX_MINUS_MAG], class_df["I magnitude"], bins="log")
    ax.annotate(f"n={len(class_df):,}", xy=(0.99, 0.99), xycoords="axes fraction", ha="right", va="top")
    plt.colorbar(cb)
    ax.axvline(-2.5, color="k", linewidth=1)
    ax.set_xlabel(MUMAX_MINUS_MAG)
    ax.set_xlim(-5, 5)
    ax.set_ylabel("I magnitude")
    ax.set_ylim(15, 27)
```

mer_mumax_minus_mag is the peak surface brightness above the background minus the magnitude that was used to compute mer_point_like_prob.
Objects to the left of the vertical line (<-2.5) are point-like.
Stars are highly concentrated there, especially those that are not faint (I < 24.5), which we should expect given Euclid's requirement for a pure sample.
Also as we should expect, most galaxies appear to the right of this line.
However, notice the strip of bright (e.g., I < 23) "galaxies" that are point-like.
Many of these are likely to be mis-classified stars or QSOs, and we could increase the purity of a galaxy sample by excluding them (at the expense of completeness, which was the original goal).
In the QSO panel, however, we see only a small cluster in the bright and point-like region where we would expect them to be.
These objects are likely to be correctly classified and are located almost exclusively in EDF-N -- likely due to the additional availability of *u* band data from UNIONS in the north.
The remaining QSO classifications (the majority) should be considered doubtful.
Many QSOs are likely to be missing from the expected region due to the overlap of QSOs with galaxies in relevant color spaces and the relatively high probability thresholds imposed for QSOs.

+++

## 6. Magnitudes

+++

Euclid Q1 contains two main types of flux measurements -- aperture and template-fit -- provided for both Euclid and external survey bands.
In this section, we look at the Q1 magnitude distributions in the four Euclid bands as a function of PHZ class and we compare the two types of flux measurements.
Additional flux measurements that will be of interest to some, but that we don't look at here, include: PSF-fit fluxes (VIS only); Sérsic-fit fluxes (computed for parametric morphology, discussed below); and fluxes that were corrected based on PHZ class.

```{code-cell}
# Construct a basic filter.
mag_filter = (
    (pc.field(VIS_DET) == 1)  # No NIR-only objects. For simpler total-flux calculations; not required.
    & (pc.field(PHZ_CLASS).isin([1, 2, 3, 4, 5, 6, 7]))  # Stars, Galaxies, QSOs, and mixed classes.
    & (pc.field(SPURIOUS_FLAG) == 0)  # Basic quality cut.
)
```

Template fluxes are expected to be more accurate than aperture fluxes for extended sources because the templates do a better job of excluding contamination from nearby sources.
Conversely, aperture fluxes were found to be more accurate for point-like objects (esp. bright stars) in the Q1 NIR photometry, likely due to better handling of PSF related issues.
In either case, Euclid recommends scaling the measured fluxes in the NIR bands with a color term in order to obtain the best estimate of the total flux in that band, which we demonstrate here (see also: Romelli; MER Photometry Cookbook).

```{code-cell}
# Columns we actually want to load. Dict instead of list because we're defining new columns (magnitudes).
# We'll have PyArrow return only magnitudes so that we don't have to handle all the flux columns in memory ourselves.
_mag_columns = {
    # FLUX_TOTAL is the best estimate for the total flux in the detection band (here, VIS) and comes from
    # aperture photometry. VIS provides the template for NIR bands. It has no unique templfit flux itself.
    "I total": flux_to_magnitude(FLUX_TOTAL),
    # Template-fit fluxes.
    "Y templfit total": flux_to_magnitude("mer_flux_y_templfit", color_col_names=(FLUX_TOTAL, "mer_flux_vis_to_y_templfit")),
    "J templfit total": flux_to_magnitude("mer_flux_j_templfit", color_col_names=(FLUX_TOTAL, "mer_flux_vis_to_j_templfit")),
    "H templfit total": flux_to_magnitude("mer_flux_h_templfit", color_col_names=(FLUX_TOTAL, "mer_flux_vis_to_h_templfit")),
    # Aperture fluxes.
    "Y aperture total": flux_to_magnitude("mer_flux_y_2fwhm_aper", color_col_names=(FLUX_TOTAL, "mer_flux_vis_2fwhm_aper")),
    "J aperture total": flux_to_magnitude("mer_flux_j_2fwhm_aper", color_col_names=(FLUX_TOTAL, "mer_flux_vis_2fwhm_aper")),
    "H aperture total": flux_to_magnitude("mer_flux_h_2fwhm_aper", color_col_names=(FLUX_TOTAL, "mer_flux_vis_2fwhm_aper")),
}
mag_columns = {**_mag_columns, PHZ_CLASS: pc.field(PHZ_CLASS), MUMAX_MINUS_MAG: pc.field(MUMAX_MINUS_MAG)}
```

Load data.

```{code-cell}
mags_df = dataset.to_table(columns=mag_columns, filter=mag_filter).to_pandas()
```

Given Euclid's core science goals, we'll take the template fluxes as our baseline in this section.

Plot the magnitude distributions as a function of PHZ class.
Include multiply-classed objects and separate point-likes, given what we learned above.

```{code-cell}
# Galaxy + any. Star + galaxy. QSO + galaxy.
classes = {"Galaxy": (2, 3, 6, 7), "Star": (1, 3), "QSO": (4, 6)}
class_colors = ["tab:green", "tab:blue", "tab:orange"]

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

    # Get the objects that were accepted as multiple classes.
    class_df = mags_df.loc[mags_df[PHZ_CLASS].isin(class_ids)]
    label = "+Galaxy" if class_name != "Galaxy" else "+any"
    # Of those objects, restrict to the ones that are point-like.
    classpt_df = class_df.loc[class_df[MUMAX_MINUS_MAG] < -2.5]
    # Plot histograms for both sets of objects.
    for ax, band in zip(axs, bands):
        ax.hist(class_df[band], color=class_color, label=label, linestyle=":", **hist_kwargs)
        ax.hist(classpt_df[band], color=class_color, linestyle="-.", label=label + " and point-like", **hist_kwargs)

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
In the green dash-dot line we see the population of point-like "galaxies", probably misclassified but mostly at faint magnitudes.
The star distributions are broader and peak at a brighter magnitude than the galaxy distributions, as we expect.
Adding to the star-only class the objects that were classified as both star and galaxy adds a significant number, especially to the peak of the distribution and towards the right where we might expect more confusion due to faint sources.
Restricting these to point-like objects, we see that many of the bright objects that surpassed both the star and galaxy probability thresholds are likely to be, in fact, stars.
(Evidence of this can also be seen in Tucci Fig. 4.)
However, this don't hold at the faint end where even some star-only classified objects are removed by the point-like cut.
The bottom panel shows very few point-like QSOs, reminding us that most QSO classifications are suspect.

+++

Now, compare with the aperture fluxes by plotting the difference between the reference (template) and aperture magnitudes.
This figure is inspired by Romelli Fig. 6 (top panel).

```{code-cell}
# Only consider objects within these mag and mag difference limits.
mag_limits, mag_diff_limits = (16, 24), (-1, 1)
mag_limited_df = mags_df.loc[(mags_df["I total"] > mag_limits[0]) & (mags_df["I total"] < mag_limits[1])]

fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharey=True, sharex=True)
bands = [("Y templfit total", "Y aperture total"), ("J templfit total", "J aperture total"), ("H templfit total", "H aperture total")]
hexbin_kwargs = dict(bins="log", extent=(*mag_limits, *mag_diff_limits), gridsize=25)
annotate_kwargs = dict(xycoords="axes fraction", ha="left", fontweight="bold", bbox=dict(facecolor="white", alpha=0.8))

# Plot
for axs, (ref_band, sersic_band) in zip(axes.transpose(), bands):
    # Extended objects, top row.
    ax = axs[0]
    extended_df = mags_df.loc[mags_df[MUMAX_MINUS_MAG] >= -2.5]
    extended_mag_diff = extended_df[ref_band] - extended_df[sersic_band]
    cb = ax.hexbin(extended_df["I total"], extended_mag_diff, **hexbin_kwargs)
    plt.colorbar(cb)
    ax.set_ylabel(f"{ref_band} - {sersic_band}")
    # Annotate top (bottom) with the fraction of objects having a magnitude difference greater (less) than 0.
    frac_tmpl_greater = len(extended_mag_diff.loc[extended_mag_diff > 0]) / len(extended_mag_diff)
    ax.annotate(f"{frac_tmpl_greater:.3f}", xy=(0.01, 0.99), va="top", **annotate_kwargs)
    frac_tmpl_less = len(extended_mag_diff.loc[extended_mag_diff < 0]) / len(extended_mag_diff)
    ax.annotate(f"{frac_tmpl_less:.3f}", xy=(0.01, 0.01), va="bottom", **annotate_kwargs)

    # Point-like objects, bottom row.
    ax = axs[1]
    pointlike_df = mags_df.loc[mags_df[MUMAX_MINUS_MAG] < -2.5]
    pointlike_mag_diff = pointlike_df[ref_band] - pointlike_df[sersic_band]
    cb = ax.hexbin(pointlike_df["I total"], pointlike_mag_diff, **hexbin_kwargs)
    plt.colorbar(cb)
    ax.set_ylabel(f"{ref_band} - {sersic_band}")
    # Annotate top (bottom) with the fraction of objects having a magnitude difference greater (less) than 0.
    frac_tmpl_greater = len(pointlike_mag_diff.loc[pointlike_mag_diff > 0]) / len(pointlike_mag_diff)
    ax.annotate(f"{frac_tmpl_greater:.3f}", xy=(0.01, 0.99), va="top", **annotate_kwargs)
    frac_tmpl_less = len(pointlike_mag_diff.loc[pointlike_mag_diff < 0]) / len(pointlike_mag_diff)
    ax.annotate(f"{frac_tmpl_less:.3f}", xy=(0.01, 0.01), va="bottom", **annotate_kwargs)

# Add axis labels, etc.
for i, ax in enumerate(axes.flatten()):
    ax.axhline(0, color="gray", linewidth=1)
    if i == 1:
        ax.set_title("Extended objects")
    if i == 4:
        ax.set_title("Point-like objects")
    if i > 2:
        ax.set_xlabel("I total")
plt.tight_layout()
```

The template - aperture magnitude difference is fairly tightly clustered around 0 for extended objects (top row) but the outliers are asymmetric (fractions above and below zero are noted).
We see a positive offset which indicates a fainter template-fit magnitude, as we should expect given that the templates do a better job of excluding contaminating light from nearby sources.
The offset is more pronounced for point-like objects (bottom row), likely due to the PSF handling mentioned above, and we are reminded that aperture magnitudes are more reliable here.

+++

## 7. Galaxy morphology

+++

The MORPH table includes three kinds of galaxy morphology indicators:
non-parametric parameters like Smoothness and Concentration,
parametric parameters like Sérsic index and axis ratio,
and Zoobot responses to questions like "does it have spiral arms?".
Zoobot is a deep learning model that was trained on Galaxy Zoo classifications.

Here, we compare one parameter of each type.
This is based on [Euclid Collaboration: Quilley et al., 2025](https://arxiv.org/pdf/2503.15309) (hereafter, Quilley).

```{code-cell}
# Non-parametric column
CONCENTRATION = "morph_concentration"
# Sérsic
SERSIC_VIS_INDEX = "morph_sersic_sersic_vis_index"
# Zoobot columns use the syntax {question}_{answer}.
zoo_question, zoo_answers = "morph_smooth_or_featured", ["smooth", "featured_or_disk", "artifact_star_zoom"]
zoo_smooth_columns = [f"{zoo_question}_{answer}" for answer in zoo_answers]

morph_filter = (
    # Cuts described in Quilley sec. 2.
    (pc.field(PHZ_CLASS) == 2)  # Galaxies
    & (pc.field(VIS_DET) == 1)
    & (pc.field(FLUX_TOTAL) > magnitude_to_flux(23))  # I<23 recommended for reliable Sérsic fits.
    & (pc.field(SPURIOUS_FLAG) == 0)
    & (pc.field(POINTLIKE_PROB) <= 0.1)
    # Sec. 4. Remove an artificial peak at the limit of the param space. Recommended for any Sérsic-based analysis.
    & (pc.field(SERSIC_VIS_INDEX) <= 5.45)
    # Secs. 4 & 5 make additional quality cuts that we skip for simplicity.
)

# Columns we'll load.
morph_columns = [CONCENTRATION, SERSIC_VIS_INDEX, *zoo_smooth_columns, FLUX_TOTAL]
```

```{code-cell}
# Load data.
galaxies = dataset.to_table(columns=morph_columns, filter=morph_filter).to_pandas()
```

Transform the Zoobot columns to “the fraction of volunteers expected to select this morphology answer”, following the MER Morphology Cookbook.

```{code-cell}
galaxies["zoo_smooth_total"] = galaxies[zoo_smooth_columns].sum(axis=1)
zoo_galaxies = galaxies.loc[galaxies.zoo_smooth_total > 0]
zoo_galaxies["featured_or_disk"] = zoo_galaxies[f"{zoo_question}_featured_or_disk"] / zoo_galaxies.zoo_smooth_total
zoo_galaxies["smooth"] = zoo_galaxies[f"{zoo_question}_smooth"] / zoo_galaxies.zoo_smooth_total
```

Plot the comparisons.
This is a combination of Quilley fig. 4 and the top left panel of fig. 5.

```{code-cell}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# ax1, Sérsic index vs Concentration
bright_galaxies = galaxies.loc[galaxies[FLUX_TOTAL] > magnitude_to_flux(21)]
bright_galaxies.plot.hexbin(SERSIC_VIS_INDEX, CONCENTRATION, ax=ax1, bins="log", extent=(0, 5, 1, 5))

# ax2, Sérsic index vs Zoobot separated by probability for each answer.
labels = ["featured_or_disk > 0.5", "featured_or_disk > 0.9", "smooth > 0.5", "smooth > 0.9"]
colors = ["tab:cyan", "tab:blue", "tab:pink", "tab:red"]
data = [
    zoo_galaxies.loc[zoo_galaxies["featured_or_disk"] > 0.5, SERSIC_VIS_INDEX],
    zoo_galaxies.loc[zoo_galaxies["featured_or_disk"] > 0.9, SERSIC_VIS_INDEX],
    zoo_galaxies.loc[zoo_galaxies["smooth"] > 0.5, SERSIC_VIS_INDEX],
    zoo_galaxies.loc[zoo_galaxies["smooth"] > 0.9, SERSIC_VIS_INDEX],
]
for sersic_indexes, label, color in zip(data, labels, colors):
    # Histogram of the fraction of galaxies in this category.
    weights = [1 / len(sersic_indexes)] * len(sersic_indexes)
    ax2.hist(sersic_indexes, label=label, color=color, histtype="step", bins=20, weights=weights)
ax2.legend()
ax2.set_xlabel(SERSIC_VIS_INDEX)
ax2.set_ylabel("Fraction of galaxies")

plt.tight_layout()
```

The left panel shows the expected relationship, since the concentration is known to correlate strongly with Sérsic index (Romelli).
The differences with Quilley fig. 4 are due to the cuts that we skipped for simplicity.
The right panel also largely agrees with expectations.
"featured_or_disk" galaxies are likely to be spirals (late-type) and these typically have a Sérsic index around 1.
"smooth" galaxies are likely to be ellipticals (early-type) which typically have a Sérsic index of 4 or more, though here we need a higher threshold on the Zoobot reponse to see it.

+++

## 8. NIR-only detections: high-redshift galaxy or nearby brown dwarf?

+++

More than 20% of Q1 MER Objects were detected only in the NIR-stack mosaic, not in VIS.
We mostly ignored them above, either for quality control or simplicity.
These NIR-only detections received special handling in many areas of the Euclid processing pipelines.
They may also require careful handling by the user, some of which we pointed to above.
Here, we explore the population through data in the PHYSPARAMNIR table, produced by the PHZ branch tasked with determining physical properties of NIR-only objects, and then take the opportunity to look at redshift PDFs from PHZ, PHYSPARAMNIR, and SPE.

```{code-cell}
# Simple filter for non-spurious NIR-only objects.
nironly_nonspurious_filter = (pc.field(VIS_DET) == 0) & (pc.field(SPURIOUS_FLAG) == 0)
```

NIR-only objects are a mix of nearby brown dwarfs and high-redshift galaxies & quasars.
These two broad but very different types of objects overlap in relevant color spaces and can be difficult to separate.
(Weaver et al., 2024).
Spectra will often be required to confirm membership, but we can use photometric properties produced by PHZ to make some useful cuts first.
We'll track the following three objects to illustrate:

- object_id: -523574860290315045. **T4 dwarf**, discovered spectroscopically ([Dominguez-Tagle et al., 2025](https://arxiv.org/abs/2503.22442)).
- object_id: -600367386508373277. **L-type dwarf**, spectroscopically confirmed ([Zhang, Lodieu, and Martín, 2024](https://arxiv.org/abs/2403.15288) Table C.2. '04:00:08.99 −50:50:14.4'. Found in Q1 via cone search; separation = 1.6 arcsec).
- object_id: -531067351279302418. **Star-forming galaxy at z=5.78**, spectroscopically confirmed ([Bunker et al., 2003](https://arxiv.org/abs/astro-ph/0302401). Found in Q1 via cone search; separation = 0.59 arcsec).

```{code-cell}
targets = {
    -523574860290315045: ("T dwarf", "tab:orange"),
    -600367386508373277: ("L dwarf", "tab:brown"),
    -531067351279302418: ("Galaxy w/ spec-z=5.78", "tab:cyan"),
}
```

The following three columns will be especially helpful:

- "physparamnir_z" is the best fitted redshift point estimate that was produced by Phosphoros (same as phz_phz_median) using a configuration specific to NIR-only sources.
  Of note: the redshift range was extended up to z=10 and the input models included stars, galaxies, and QSOs.
- "physparamnir_high_z_prob" gives the probability that the redshift is greater than 6.
  It is the integral of the PDF at z>6.
- "phz_best_chi2" is the Chi^2 associated with the best fit model.

```{code-cell}
PPNIR_Z = "physparamnir_z"
PPNIR_Z_GT6_PROB = "physparamnir_high_z_prob"
PHZ_BEST_CHI2 = "phz_best_chi2"
```

Load all non-spurious NIR-only objects.

```{code-cell}
nironly_columns = [OBJECT_ID, PPNIR_Z, PPNIR_Z_GT6_PROB, PHZ_BEST_CHI2]
nironly_df = dataset.to_table(columns=nironly_columns, filter=nironly_nonspurious_filter).to_pandas()
```

Plot the best fitted redshift vs the probability of z>6 to get a sense of which regions indicate confidence vs confusion.
This is based on Tucci sec. 6.4.

```{code-cell}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6), sharex=True, sharey=True)
x, y = PPNIR_Z, PPNIR_Z_GT6_PROB

# Plot z vs prob(z>6), colored by density.
nironly_df.plot.hexbin(x, y, bins="log", ax=ax1, cmap="Greys")
ax1.set_title("Number of NIR-only objects")

# Remove top 5% of objects with worst Chi^2. Plot the rest, colored by mean Chi^2.
nironly_trim_chi2 = nironly_df.loc[nironly_df[PHZ_BEST_CHI2] < 200]
nironly_trim_chi2.plot.hexbin(x, y, ax=ax2, cmap="viridis_r", C=PHZ_BEST_CHI2, reduce_C_function=np.mean)
ax2.set_title("NIR-only objects (trimmed) colored by mean Chi^2")

# Add some indicator lines and symbols.
for ax in (ax1, ax2):
    for target_id, (target_name, target_color) in targets.items():
        target = nironly_df.loc[nironly_df[OBJECT_ID] == target_id]
        ax.scatter(target[x], target[y], marker="X", s=100, c=target_color, label=target_name)
        if ax == ax1:
            ax.legend(loc=10)
    ax.axvline(6, color="k", linewidth=1)
    ax.axhline(0.8, color="k", linewidth=1)
    ax.scatter(6, 0.8, marker="*", s=350, edgecolors="k", facecolors="none")
    ax.scatter(0, 0.78, marker="^", s=350, edgecolors="k", facecolors="none")
```

First, we notice that the redshifts for two of our target objects are close to what we'd expect (T dwarf at 0 and galaxy at 5.75), but the L dwarf's best fit redshift is much too high and, worse, the probability is relatively high as well.
However, according to Tucci, the L dwarf's position in this plane (below and right of the black star) actually indicates that its PDF is likely either broad or multi-peaked, so we should be skeptical of the redshift point estimate.
Good candidates for high-redshift objects sit in the upper right quadrant, z>6 and prob(z>6) > 0.8, and have Chi^2 values that are not much larger than 10.
(Good Chi^2 is yellow in the right panel, but note that the color represents a mean for the bin and not the Chi^2 of individual objects.)
A third region identified by Tucci is near or above the black triangle: z=0 and a high probability for z>6.
This combination indicates model degeneracy between brown dwarfs and high-redshift galaxies.

We want to look at the redshift PDFs for our objects, but first we should inspect the classifications a little more, so we'll load those too.

```{code-cell}
# Redshift PDFs that we'll look at.
PHZ_ZPDF = "phz_phz_pdf"
PPNIR_ZPDF = "physparamnir_z_1d_pdf"
SPE_QSO_ZPDF = "z_qso_candidates_spe_pdf_rank0"

# Columns to load.
targets_columns = [
    OBJECT_ID,
    # PHZ classifications
    PHZ_CLASS,
    "class_phz_gal_prob",
    "class_phz_qso_prob",
    "class_phz_star_prob",
    # PHZ redshifts
    PHZ_Z,
    PHYSPARAM_GAL_Z,
    PPNIR_Z,
    PHZ_ZPDF,
    PPNIR_ZPDF,
    # SPE classifications
    "z_spe_class",
    "z_spe_star_prob",
    "z_spe_gal_prob",
    "z_spe_qso_prob",
    # SPE redshifts
    SPE_GAL_Z,
    SPE_GAL_Z_PROB,
    "z_qso_candidates_spe_z_rank0",
    "z_qso_candidates_spe_z_prob_rank0",
    SPE_QSO_ZPDF,
]
```

```{code-cell}
# Load data.
targets_filter = pc.field(OBJECT_ID).isin(targets.keys())
targets_df = dataset.to_table(columns=targets_columns, filter=targets_filter).to_pandas()
```

```{code-cell}
# Add our target names, then look at the data.
targets_df.insert(1, "target_name", targets_df[OBJECT_ID].apply(lambda oid: targets[oid][0]))
targets_df
```

PHZ classified our galaxy as a galaxy, L dwarf as a QSO, and was unable to classify the T dwarf but gave the highest probability for a star.
For the galaxy, the main photo-z estimate (phz_phz_median) is roughly close to the spectroscopic redshift but the one produced along with galaxy physical properties (physparam_phz_pp_median_redshift) is much closer, and the one produced by the NIR-only branch (physparamnir_z) is closer still.
SPE classified the T dwarf as a QSO with probability > 0.999, but we should still be skeptical (even if we didn't know the truth already) because both the object type and the redshift solution are outside of its target range.
It was unable to produce classifications or even probabilities for our other two targets.
SPE will have produced zPDFs assuming both galaxy and QSO for all objects, regardless of the actual classification, but we'll ignore all of them except the QSO solution for the T dwarf.

Plot the redshift PDFs.

```{code-cell}
# PHZ z<6 step is 0.01 and z in [6, 10] step is 0.05.
phz_zbins = np.concatenate([np.arange(0, 6, 0.01), np.arange(6, 10.01, 0.05)])
step_kwargs = dict(where="mid", alpha=0.7)
ppnir_color, phz_color, spe_color = "tab:orange", "tab:cyan", "tab:purple"

fig, axes = plt.subplots(1, 3, figsize=(24, 6))
for ax, (target_id, (target_name, target_color)) in zip(axes, targets.items()):
    obj_df = targets_df.loc[targets_df[OBJECT_ID] == target_id]
    # PHYSPARAMNIR range extends to z=10.
    ax.step(phz_zbins, obj_df[PPNIR_ZPDF].squeeze(), label=PPNIR_ZPDF, zorder=9, color=ppnir_color, **step_kwargs)
    ax.scatter(obj_df[PPNIR_Z], 0, label=PPNIR_Z, marker="*", color=ppnir_color)
    # PHZ range stops at z=6.
    ax.step(phz_zbins[:601], obj_df[PHZ_ZPDF].squeeze(), label=PHZ_ZPDF, color=phz_color, **step_kwargs)
    ax.scatter(obj_df[PHZ_Z], 0, label=PHZ_Z, marker="*", color=phz_color)

    ax.legend()
    ax.set_xlabel("Redshift")
    ax.set_title(target_name)

    # [FIXME] SPE... how to plot? 'Z_QSO_CANDIDATES_SPE_PDF_RANK0' is a PDF array of length 80, but what is its redshift range?
    # There are also columns for 'ZMIN', 'ZMAX', and 'DELTAZ' but naively applying them (commented out below)
    # produces nonsense and I can't find documentation that gives instructions.
    # if target_name == "T dwarf":
    #     start, step = obj_df["Z_QSO_CANDIDATES_SPE_PDF_ZMIN_RANK0"].squeeze(), obj_df["Z_QSO_CANDIDATES_SPE_PDF_DELTAZ_RANK0"].squeeze()
    #     spe_zbins = start + step * np.arange(80)
    #     ax.step(spe_zbins, obj_df["Z_QSO_CANDIDATES_SPE_PDF_RANK0"].squeeze(), label="SPE QSO zPDF", **step_kwargs)
```

In the left panel (T dwarf), we see that the photo-z PDF produced by the NIR-only branch is very strongly peaked at z=0 with a small secondary bump near z=7, consistent with its placement in the previous figure.
Recall that PHZ_PHZ_PDF was produced using galaxy models, regardless of the object's class.
In the middle panel (L dwarf), we see the broad and multi-peaked NIR PDF that was guessed at based on the previous figure.
While the strongest peak is near z=8 (a QSO solution, perhaps?), there are also peaks at z=0 (consistent with a star) and near z=1.75 (consistent with the PHZ (galaxy) solution) which are prominent enough to reduce the probability of z>6 below the 0.8 threshold.
In the right panel (Galaxy), we see good agreement between the PDFs except at z=0.

[FIXME] Once there is spectra in this hats product we can look at those too.

+++

## Appendix: Schema details

This Euclid Q1 HATS Catalog contains the 14 Euclid tables listed in the introduction, joined on 'object_id' into a single parquet dataset.
In addition, the Euclid 'tileid' for each object has been added, as well as a few HATS- and HEALPix-related columns.
All Euclid column names have been lower-cased and the table name has been prepended (e.g., 'FLUX_H_TEMPLFIT' -> 'mer_flux_h_templfit'), except for the following:

- object_id : Euclid MER Object ID. Unique identifier of a row in this dataset.
- tileid : ID of the Euclid Tile the object was detected in.
- ra : 'RIGHT_ASCENSION' from the 'mer' table. Named shortened to match other IRSA services.
- dec : 'DECLINATION' from the 'mer' table. Named shortened to match other IRSA services.

In addition to the above changes, all non-alphanumeric characters in column names have been replaced with an underscore for compatibility with various libraries and services (e.g., 'E(B-V)' -> 'physparamqso_e_b_v_').
Finally, the SPE tables 'z', 'lines', and 'models' required special handling as follows:

- z : The original FITS files contain the spectroscopic redshift estimates for GALAXY_CANDIDATES, STAR_CANDIDATES, and QSO_CANDIDATES (HDUs 3, 4, and 5 respectively) with up to 5 estimates per 'object_id', per HDU.
  For the parquet dataset, these were pivoted so that there is one row per 'object_id' in order to facilitate the table joins.
  The resulting columns were named by combining the table name (z), the HDU name, the original column name, and the rank of the given redshift estimate (i.e., the value in the original 'SPE_RANK' column).
  For example, the 'SPE_PDF' column for the highest ranked redshift estimate in the 'GALAXY_CANDIDATES' table is called 'z_galaxy_candidates_spe_pdf_rank0'.
- lines : The parquet dataset only includes the rows from HDU1 with 'SPE_LINE_NAME' == 'Halpha'.
  Similar to above, there are up to 5 sets of columns per 'object_id', one per redshift estimate.
  Column names have been appended with both the rank and the line name.
  For example, the column originally called 'SPE_LINE_FLUX_GF' is named 'lines_spe_line_flux_gf_rank0_halpha' for the Halpha line identified with the highest ranked redshift estimate.
- models : The parquet dataset only includes HDU2 -- the model parameters for the galaxy solutions.
  This table has the same structure as 'z'.
  In addition to the table name, 'galaxy' has been appended to the column names.
  For example, the column originally called 'SPE_VEL_DISP_E' is named 'models_galaxy_spe_vel_disp_e_rank0' for the velocity dispersion of emission lines needed to fit the highest ranked galaxy redshift estimate.

Below, we follow IRSA's
[Cloud Access notebook](https://caltech-ipac.github.io/irsa-tutorials/tutorials/cloud_access/cloud-access-intro.html#navigate-a-catalog-and-perform-a-basic-query)
to inspect the parquet schema.
The schema is accessible from the PyArrow dataset object we already loaded (`dataset.schema`) but we load it here from the `_common_metadata` file because it includes column metadata (units and descriptions) that is not present in `dataset.schema`.

```{code-cell}
schema = pyarrow.parquet.read_schema(euclid_parquet_schema_path, filesystem=s3_filesystem)
```

There are almost 1600 columns in this dataset.

```{code-cell}
print(f"{len(schema)} columns total")
```

### A.1 Euclid Q1 columns

+++

To find all columns from a given table, search for column names that start with the table name followed by an underscore.
Table names are given in section 1.1.

```{code-cell}
# Find all column names from the phz table.
phz_columns = [name for name in schema.names if name.startswith("phz_")]

print(f"{len(phz_columns)} columns from the PHZ table. First four are:")
phz_columns[:4]
```

Column metadata includes unit and description.

```{code-cell}
schema.field("mer_flux_y_2fwhm_aper").metadata
```

Euclid Q1 offers many flux measurements, both from Euclid detections and from external ground-based surveys.
They are given in microjanskys, so all flux columns can be found by searching the metadata for this unit.

```{code-cell}
# Find all flux columns.
flux_columns = [field.name for field in schema if field.metadata[b"unit"] == b"uJy"]

print(f"{len(flux_columns)} flux columns. First four are:")
flux_columns[:4]
```

Columns associated with external surveys are identified by the inclusion of "ext" in the name.

```{code-cell}
external_flux_columns = [name for name in flux_columns if "ext" in name]
print(f"{len(external_flux_columns)} flux columns from external surveys. First four are:")
external_flux_columns[:4]
```

### A.2 Additional columns

+++

Several columns have been added to this dataset that are not in the original Euclid Q1 tables:

- 'tileid' : Euclid MER tile index.
- '_healpix_9' : HEALPix order 9 pixel index. Useful for spatial queries.
- '_healpix_19' : HEALPix order 19 pixel index. Useful for spatial queries.
- '_healpix_29' : (hats column) HEALPix order 29 pixel index. Useful for spatial queries.

```{code-cell}
schema.names[:5]
```

In addition, the following columns are used for the HATS partitioning:

- 'Norder' : (hats column) HEALPix order at which the data is partitioned.
- 'Npix' : (hats column) HEALPix pixel index at order Norder.
- 'Dir' : (hats column) Integer equal to 10_000 * floor[Npix / 10_000].

They appear at the end of the schema and are also accessible from the PyArrow dataset partitioning object.

```{code-cell}
dataset.partitioning.schema
```

## About this notebook

**Authors:** Troy Raen (Developer; Caltech/IPAC-IRSA) and the IRSA Data Science Team.

**Updated:** 2025-06-30

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.
