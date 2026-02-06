---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Euclid Q1 Merged Objects HATS Catalog: Introduction

+++

This tutorial is an introduction to the content and format of the Euclid Q1 Merged Objects HATS Catalog.
Later tutorials in this series will show how to load quality samples.
See [Euclid Tutorial Notebooks](../../euclid_access/euclid.md) for a complete list of Euclid tutorials.

+++

## Learning Goals

In this tutorial, we will:

- Learn about the Euclid Merged Objects catalog that IRSA created by combining information from multiple Euclid Quick Release 1 (Q1) catalogs.
- Find columns of interest.
- Perform a basic query using the Python library PyArrow.

+++

## 1. Introduction

The [Euclid Q1](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_q1.html) catalogs were derived from Euclid photometry and spectroscopy, taken by the Visible Camera (VIS) and the Near-Infrared Spectrometer and Photometer (NISP), and from photometry taken by other ground-based instruments.
The data include several flux measurements per band, several redshift estimates, several morphology parameters, etc.
Each was derived for different science goals using different algorithms or configurations.

The Euclid Q1 Merged Objects HATS Catalog was produced by IRSA by joining 14 of the original catalogs on object ID (column: `object_id`).
Following the Hierarchical Adaptive Tiling Scheme [HATS](https://hats.readthedocs.io/) framework, the data were then partitioned spatially (by right ascension and declination) and written as an Apache Parquet dataset.

- Columns: 1,594
- Rows: 29,953,430 (one per Euclid Q1 object)
- Size: 400 GB

The catalog is served from an AWS S3 cloud storage bucket.
Access is free and no credentials are required.

+++

## 2. Imports

```{code-cell} ipython3
# # Uncomment the next line to install dependencies if needed.
# %pip install hpgeom pandas pyarrow
```

```{code-cell} ipython3
import hpgeom  # Find HEALPix indexes from RA and Dec
import pyarrow.compute as pc  # Filter the catalog
import pyarrow.dataset  # Load the catalog
import pyarrow.fs  # Simple S3 filesystem pointer
import pyarrow.parquet  # Load the schema
```

## 3. Load Parquet Metadata

First we'll load the Parquet schema (column information) of the Merged Objects catalog so we can use it in later sections.
The Parquet schema is accessible from a few locations, all of which include the column names and types.
Here, we load it from the `_common_metadata` file because it also includes the column units and descriptions.

```{code-cell} ipython3
# AWS S3 paths.
s3_bucket = "nasa-irsa-euclid-q1"
dataset_prefix = "contributed/q1/merged_objects/hats/euclid_q1_merged_objects-hats/dataset"

dataset_path = f"{s3_bucket}/{dataset_prefix}"
schema_path = f"{dataset_path}/_common_metadata"

# S3 pointer. Use `anonymous=True` to access without credentials.
s3 = pyarrow.fs.S3FileSystem(anonymous=True)
```

```{code-cell} ipython3
# Load the Parquet schema.
schema = pyarrow.parquet.read_schema(schema_path, filesystem=s3)

# There are almost 1600 columns in this dataset.
print(f"{len(schema)} columns in the Euclid Q1 Merged Objects catalog")
```

## 4. Merged Objects Catalog Contents

+++

The Merged Objects catalog contains data from 14 Euclid Q1 tables, joined on the column `object_id`.
The tables were produced by three Euclid processing functions: MER (multi-wavelength mosaics on common spatial and pixel scales), PHZ (photometric redshifts), and SPE (spectroscopy).
The subsections below include the table names, links to reference papers, URLs to the original table schemas, and examples of how the column names were transformed for the Merged Objects catalog.

The original tables' column names are mostly in all caps.
In the Merged Objects catalog and the catalogs available through IRSA's TAP service, all column names have been lower-cased.
In addition, all non-alphanumeric characters have been replaced with an underscore for compatibility with various libraries and services.
Finally, the original table name has been prepended to column names in the Merged Objects catalog, both for provenance and to avoid duplicates.
An example that includes all of these transformations is: `E(B-V)` -> `physparamqso_e_b_v_`.

Three columns have special names that differ from the standard naming convention described above:

- `object_id` : Euclid MER Object ID. Unique identifier of a row in this dataset. No table name prepended.
- `ra` : 'RIGHT_ASCENSION' from the 'mer' table. Name shortened to match other IRSA services. No table name prepended.
- `dec` : 'DECLINATION' from the 'mer' table. Name shortened to match other IRSA services. No table name prepended.

Seven additional columns have been added to the Merged Objects catalog that are not in the original Euclid tables.
They are described below, after the Euclid tables.

### 4.1 MER tables

The Euclid MER processing function produced three tables.
The reference paper is [Euclid Collaboration: Romelli et al., 2025](https://arxiv.org/pdf/2503.15305) (hereafter, Romelli).
The tables are:

**Main table (mer)**
- Described in Romelli sections 6 & 8 ("EUC_MER_FINAL-CAT")
- Original schema: [Main catalog FITS file](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#main-catalog-fits-file)
- Example column name transform: `FLUX_DETECTION_TOTAL` --> `mer_flux_detection_total`

**Morphology (morph)**
- Described in Romelli sections 7 & 8 ("EUC_MER_FINAL-MORPH-CAT")
- Original schema: [Morphology catalog FITS file](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#morphology-catalog-fits-file)
- Example column name transform: `CONCENTRATION` --> `morph_concentration`

**Cutouts (cutouts)**
- Described in Romelli section 8 ("EUC_MER_FINAL-CUTOUTS-CAT")
- Original schema: [Cutouts catalog FITS file](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html#cutouts-catalog-fits-file)
- Example column name transform: `CORNER_0_RA` --> `cutouts_corner_0_ra`

Find all columns from these tables in the Parquet schema:

```{code-cell} ipython3
mer_prefixes = ["mer_", "morph_", "cutouts_"]
mer_col_counts = {p: len([n for n in schema.names if n.startswith(p)]) for p in mer_prefixes}

print(f"MER tables: {sum(mer_col_counts.values())} columns total")
for prefix, count in mer_col_counts.items():
    print(f"  {prefix}: {count}")
```

### 4.2 PHZ tables

The Euclid PHZ processing function produced eight tables.
The reference paper is [Euclid Collaboration: Tucci et al., 2025](https://arxiv.org/pdf/2503.15306) (hereafter, Tucci).
The tables are:

**Photometric Redshifts (phz)**
- Described in Tucci section 5 ("phz_photo_z")
- Original schema: [Photo Z catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#photo-z-catalog)
- Example column name transform: `PHZ_PDF` --> `phz_phz_pdf`

**Classifications (class)**
- Described in Tucci section 4 ("phz_classification")
- Original schema: [Classification catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#classification-catalog)
- Example column name transform: `PHZ_CLASSIFICATION` --> `class_phz_classification`

**Galaxy Physical Parameters (physparam)**
- Described in Tucci section 6 (6.1; "phz_physical_parameters")
- Original schema: [Physical Parameters catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#physical-parameters-catalog)
- Example column name transform: `PHZ_PP_MEDIAN_REDSHIFT` --> `physparam_phz_pp_median_redshift`

**Galaxy SEDs (galaxysed)**
- Described in Tucci appendix B (B.1 "phz_galaxy_sed")
- Original schema: [Galaxy SED catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#galaxy-sed-catalog)
- Example column name transform: `FLUX_4900_5000` --> `galaxysed_flux_4900_5000`

**QSO Physical Parameters (physparamqso)**
- Described in Tucci section 6 (6.2; "phz_qso_physical_parameters")
- Original schema: [QSO Physical Parameters catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#qso-physical-parameters-catalog)
- Example column name transform: `E(B-V)` --> `physparamqso_e_b_v_`

**Star Parameters (starclass)**
- Described in Tucci section 6 (6.3; "phz_star_template")
- Original schema: [Star Template](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#star-template)
- Example column name transform: `FLUX_VIS_Total_Corrected` --> `starclass_flux_vis_total_corrected`

**Star SEDs (starsed)**
- Described in Tucci appendix B (B.1 "phz_star_sed")
- Original schema: [Star SED catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputcatalog.html#star-sed-catalog)
- Example column name transform: `FLUX_4900_5000` --> `starsed_flux_4900_5000`

**NIR Physical Parameters (physparamnir)**
- Described in Tucci section 6 (6.4; "phz_nir_physical_parameters")
- Original schema: [NIR Physical Parameters catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/phzdpd/dpcards/phz_phzpfoutputforl3.html#nir-physical-parameters-catalog)
- Example column name transform: `E(B-V)` --> `physparamnir_e_b_v_`

Find all columns from these tables in the Parquet schema:

```{code-cell} ipython3
phz_prefixes = ["phz_", "class_", "physparam_", "galaxysed_", "physparamqso_",
                "starclass_", "starsed_", "physparamnir_"]
phz_col_counts = {p: len([n for n in schema.names if n.startswith(p)]) for p in phz_prefixes}

print(f"PHZ tables: {sum(phz_col_counts.values())} columns total")
for prefix, count in phz_col_counts.items():
    print(f"  {prefix}: {count}")
```

### 4.3 SPE tables

The Euclid SPE processing function produced three tables from which data are included in the Merged Objects catalog.
The reference paper is [Euclid Collaboration: Le Brun et al., 2025](https://arxiv.org/pdf/2503.15308) (hereafter, Le Brun).

These tables required special handling because they contain multiple rows per object (identified by column `object_id`).
The tables were pivoted before being joined so that the Merged Objects catalog contains one row per object.
The pivoted columns were named by combining at least the table name, the original column name, and the rank of the redshift estimate (i.e., the value in the original 'SPE_RANK' column).

The tables are:

**Spectroscopic Redshifts (z)**
- Described in Le Brun section 2 ("spectro_zcatalog_spe_quality", "spectro_zcatalog_spe_classification", "spectro_zcatalog_spe_galaxy_candidates", "spectro_zcatalog_spe_star_candidates", and "spectro_zcatalog_spe_qso_candidates")
- Original schema: [Redshift catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/spedpd/dpcards/spe_spepfoutputcatalog.html#redshift-catalog)
- Pivot and naming: For each object, the original table contains up to 5 redshift estimates (ranked by confidence) produced by assuming the object is a galaxy, star, and QSO -- for a total of up to 15 rows per object.
  The table was pivoted to one row per object and the resulting columns were named by prepending the table name (`z`) and the assumed type (`galaxy_candidates`, `star_candidates`, and `qso_candidates`), and appending the rank (`rank[0-4]`).
- Example column name transform: `SPE_PDF` --> `z_galaxy_candidates_spe_pdf_rank0` (top-ranked redshift PDF, assuming galaxy)

**Spectral Line Measurements (lines)**
- Described in Le Brun section 5 ("spectro_line_features_catalog_spe_line_features_cat").
  _Notice that lines were identified by assuming the object is a **galaxy**._
- Original schema: [Lines catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/spedpd/dpcards/spe_spepfoutputcatalog.html#lines-catalog) (HDU#1 only)
- Pivot and naming: For each object and line, the original table contains up to 5 rows, one per galaxy redshift estimate from the **z** table.
  The Merged Objects catalog only contains the Halpha line measurements.
  The table was pivoted to one row per object and the resulting columns were named by prepending the table name (`lines`) and appending both the redshift rank (`rank[0-4]`) and the name of the line (`halpha`).
- Example column name transform: `SPE_LINE_FLUX_GF` --> `lines_spe_line_flux_gf_rank0_halpha` (Halpha line flux of the top-ranked redshift, assuming galaxy)

**Models (models)**
- Described in Le Brun section 5 ("spectro_model_catalog_spe_lines_catalog")
- Original schema: [Models catalog](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/spedpd/dpcards/spe_spepfoutputcatalog.html#models-catalog) (HDU#2 only)
- Pivot and naming: The original table has the same structure as the **z** table.
  The Merged Objects catalog only contains the galaxy solutions.
  The table was pivoted to one row per object and the resulting columns were named by prepending the table name (`models`) and the assumed type (`galaxy`), and appending the redshift rank (`rank[0-4]`).
- Example column name transform: `SPE_VEL_DISP_E` --> `models_galaxy_spe_vel_disp_e_rank0` (velocity dispersion of the top-ranked redshift, assuming galaxy)

Find all columns from these tables in the Parquet schema:

```{code-cell} ipython3
spe_prefixes = ["z_", "lines_", "models_"]
spe_col_counts = {p: len([n for n in schema.names if n.startswith(p)]) for p in spe_prefixes}

print(f"SPE tables: {sum(spe_col_counts.values())} columns total")
for prefix, count in spe_col_counts.items():
    print(f"  {prefix[:-1]}: {count}")
```

### 4.4 Additional columns

The following columns were added to the Merged Objects catalog but do not appear in the original Euclid tables.

**Euclid columns:**

- `tileid` : ID of the Euclid Tile the object was detected in.
  The Euclid tiling is described in Romelli section 3.1.

**HEALPix columns:**

These HEALPix indexes correspond to the object's RA and Dec coordinates.
They are useful for spatial queries, as demonstrated in the Euclid Deep Fields section below.

- `_healpix_9` : HEALPix order 9 pixel index.
  Order 9 pixels have a resolution (square root of area) of ~400 arcseconds or ~0.1 degrees.
- `_healpix_19` : HEALPix order 19 pixel index.
  Order 19 pixels have a resolution of ~0.4 arcseconds.
- `_healpix_29` : HEALPix order 29 pixel index.
  Order 29 pixels have a resolution of ~4e-4 arcseconds.

The HEALPix, Euclid object ID, and Euclid tile ID columns appear first:

```{code-cell} ipython3
schema.names[:5]
```

**HATS columns:**

These are the HATS partitioning columns.
They appear in the Parquet file names but are not included inside the files.
However, PyArrow automatically makes them available as regular columns when the dataset is loaded as demonstrated in these tutorials.

- `Norder` : (hats column) HEALPix order at which the data is partitioned.
- `Npix` : (hats column) HEALPix pixel index at order Norder.
- `Dir` : (hats column) Integer equal to 10_000 * floor[Npix / 10_000].

The HATS columns appear at the end:

```{code-cell} ipython3
schema.names[-3:]
```

### 4.5 Find columns of interest

The subsections above show how to find all columns from a given Euclid table as well as the additional columns.
Here we show some additional techniques for finding columns.

```{code-cell} ipython3
# Access the data type using the `field` method.
schema.field("mer_flux_y_2fwhm_aper")
```

```{code-cell} ipython3
# The column metadata includes unit and description.
# Parquet metadata is always stored as bytestrings, which are denoted by a leading 'b'.
schema.field("mer_flux_y_2fwhm_aper").metadata
```

Euclid Q1 offers many flux measurements, both from Euclid detections and from external ground-based surveys.
They are given in microjanskys, so all flux columns can be found by searching the metadata for this unit.

```{code-cell} ipython3
# Find all flux columns.
flux_columns = [field.name for field in schema if field.metadata[b"unit"] == b"uJy"]

print(f"{len(flux_columns)} flux columns. First four are:")
flux_columns[:4]
```

Columns associated with external surveys are identified by the inclusion of "ext" in the name.

```{code-cell} ipython3
external_flux_columns = [name for name in flux_columns if "ext" in name]
print(f"{len(external_flux_columns)} flux columns from external surveys. First four are:")
external_flux_columns[:4]
```

## 5. Euclid Deep Fields

+++

Euclid Q1 includes data from three Euclid Deep Fields: EDF-N (North), EDF-S (South), EDF-F (Fornax; also in the southern hemisphere).
There is also a small amount of data from a fourth field: LDN1641 (Lynds' Dark Nebula 1641), which was observed for technical reasons during Euclid's verification phase.
The fields are described in [Euclid Collaboration: Aussel et al., 2025](https://arxiv.org/pdf/2503.15302) and can be seen on this [skymap](https://irsa.ipac.caltech.edu/data/download/parquet/euclid/q1/merged_objects/hats/euclid_q1_merged_objects-hats/skymap.png).

The regions are well separated, so we can distinguish them using a simple cone search without having to be too picky about the radius.
We can load data more efficiently using the HEALPix order 9 pixels that cover each area rather than using RA and Dec values directly.
These will be used in later tutorials.

```{code-cell} ipython3
# EDF-N (Euclid Deep Field - North)
ra, dec, radius = 269.733, 66.018, 4  # 20 sq deg
edfn_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius, inclusive=True)

# EDF-S (Euclid Deep Field - South)
ra, dec, radius = 61.241, -48.423, 5  # 23 sq deg
edfs_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius, inclusive=True)

# EDF-F (Euclid Deep Field - Fornax)
ra, dec, radius = 52.932, -28.088, 3  # 10 sq deg
edff_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius, inclusive=True)

# LDN1641 (Lynds' Dark Nebula 1641)
ra, dec, radius = 85.74, -8.39, 1.5  # 6 sq deg
ldn_k9_pixels = hpgeom.query_circle(hpgeom.order_to_nside(9), ra, dec, radius, inclusive=True)
```

## 6. Basic Query

To demonstrate a basic query, we'll search for objects with a galaxy photometric redshift estimate of 6.0 (largest possible).
Other tutorials in this series will show more complex queries, and describe the redshifts and other data in more detail.
PyArrow dataset filters are described at [Filtering by Expressions](https://arrow.apache.org/docs/python/compute.html#filtering-by-expressions), and the list of available functions is at [Compute Functions](https://arrow.apache.org/docs/python/api/compute.html).

```{code-cell} ipython3
dataset = pyarrow.dataset.dataset(dataset_path, partitioning="hive", filesystem=s3, schema=schema)

highz_objects = dataset.to_table(
    columns=["object_id", "phz_phz_median"], filter=pc.field("phz_phz_median") == 6
).to_pandas()
highz_objects
```

## About this notebook

**Authors:** Troy Raen, Vandana Desai, Andreas Faisst, Shoubaneh Hemmati, Jaladh Singhal, Brigitta Sipőcz, Jessica Krick, the IRSA Data Science Team, and the Euclid NASA Science Center at IPAC (ENSCI).

**Updated:** 2025-12-23

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html)
