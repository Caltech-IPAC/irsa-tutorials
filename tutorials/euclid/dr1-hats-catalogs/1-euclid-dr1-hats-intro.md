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
authors:
  - name: Deeksha Mohanty
  - name: Troy Raen
---

# Euclid DR1 HATS Catalogs: Introduction

+++

This tutorial is an introduction to the content and format of the Euclid DR1 HATS Catalogs, a set of catalogs that IRSA created from the Euclid Data Release 1 (DR1) catalog tables.

+++

## Learning Goals

In this tutorial, we will:

- Learn which Euclid DR1 catalogs are available in HATS format and how they are named and organized.
- Load the schema (column information) of a catalog and find columns of interest.
- Perform a basic query using the Python library PyArrow.
- Perform a cone search using the Python library LSDB.

+++

## 1. Introduction

Euclid DR1 contains data from the first year of survey operations, covering about 1,900 square degrees.
The DR1 catalogs were derived from Euclid photometry and spectroscopy, taken by the Visible Camera (VIS) and the Near-Infrared Spectrometer and Photometer (NISP), and from photometry taken by other ground-based instruments.

IRSA converted a subset of the Euclid DR1 catalog tables into catalogs in the Hierarchical Adaptive Tiling Scheme ([HATS](https://hats.readthedocs.io/)) format, one HATS catalog per table.
The HATS catalogs cover the Wide Survey, Wide, Deep Survey, and Deep catalog tables, and the combined spectra.
Other DR1 catalog tables, such as the Deep mode, Deep vsiits, Wide Mode tables and `phz_photo_z_wide`, are not available in HATS format.
In HATS, the data are partitioned spatially (by right ascension and declination) and written as an Apache Parquet dataset.
This is different from the Euclid Q1 release, for which IRSA joined 14 tables into a single Merged Objects HATS Catalog (see [Euclid Q1 Merged Objects HATS Catalog: Introduction](../merged-objects-hats-catalog/1-euclid-q1-hats-intro.md)).
Each DR1 HATS catalog has the same rows and columns as the original Euclid table, with a few changes described in Section 4 below.

The catalogs are served from an AWS S3 cloud storage bucket.
Access is free and no credentials are required.

+++

## 2. Imports

```{code-cell} ipython3
# # Uncomment the next line to install dependencies if needed.
# %pip install "lsdb>=0.8.1" pandas pyarrow s3fs
```

```{code-cell} ipython3
import lsdb  # Cone search with HATS catalogs
import pyarrow.compute as pc  # Filter the catalog
import pyarrow.dataset  # Load the catalog
import pyarrow.fs  # Simple S3 filesystem pointer
import pyarrow.parquet  # Load the schema
```

## 3. Available Catalogs

+++

The HATS catalogs are from the Wide Survey, Wide, Deep Survey, and Deep tables, plus the combined spectra.
The tables were produced by four Euclid processing functions: MER (multi-wavelength mosaics on common spatial and pixel scales), PHZ (photometric redshifts), SIR (NISP spectra), and SPE (spectroscopy).

The suffix in each catalog name is the suffix of the original Euclid table, and the suffixes differ between processing functions:

- MER: `_deep_survey` and `_wide_survey`
- PHZ: `_deep_survey` and `_wide`
- SPE: `_deep`, `_wide_survey`, and `_wide`

[FIXME] Explain what the suffixes mean. `_deep` and `_wide` may each combine more than one type of observation (e.g., `_deep` may include Deep Survey, Deep mode, and deep visits), but this is not confirmed. Also explain how the SPE `_wide_survey` and `_wide` tables differ.

The tables below list each catalog and its number of rows.
A dash means there is no HATS catalog for that table.

### 3.1 MER catalogs

Original schemas: [MER final catalog data product](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dm10/merdpd/dpcards/mer_finalcatalog.html)

| Table | `_deep_survey` | `_wide_survey` |
|---|---|---|
| `mer_catalogue` | 51,820,952 | 888,152,165 |
| `mer_morphology` | 51,820,952 | 888,152,165 |
| `mer_cutouts` | 51,820,952 | 888,152,165 |
| `mer_catalogue_deep_survey_photometry` (full name) | 26,008,743 | |

### 3.2 PHZ catalogs

Original schemas: [PHZ output catalog data product](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dm10/phzdpd/dpcards/phz_phzpfoutputcatalog.html) (`phz_photo_z`, `phz_galaxy_sed`, `phz_star_sed`) and [PHZ output for L3 data product](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dm10/phzdpd/dpcards/phz_phzpfoutputforl3.html) (the other PHZ catalogs)

| Table | `_deep_survey` | `_wide` |
|---|---|---|
| `phz_photo_z` | 51,550,740 | - |
| `phz_classification` | 51,550,740 | 847,313,811 |
| `phz_physical_parameters` | 51,550,740 | 847,313,811 |
| `phz_nir_physical_parameters` | 6,415,078 | 187,662,926 |
| `phz_qso_physical_parameters` | 6,154,637 | 116,564,738 |
| `phz_star_template` | 7,275,888 | 189,466,498 |
| `phz_galaxy_sed` | - | 846,826,349 |
| `phz_star_sed` | - | 12,275,216 |

### 3.3 SIR catalog

Original schema: [SIR combined spectra data product](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dm10/sirdpd/dpcards/sir_combinedspectra.html)

| Table | Rows |
|---|---|
| `sir_combined_spectra` | 309,078,775 |

The `sir_combined_spectra` catalog contains the combined 1D spectra, with one row per spectrum.
An object can have more than one spectrum, so it can appear in more than one row.
For example, an object can have both a red grism spectrum and a blue grism spectrum, which are identified by the column `lrange` (values `RGS` and `BGS`).
An object can also have more than one spectrum with the same `lrange`, so `object_id` and `lrange` together do not uniquely identify a row.
The spectrum itself is stored in a single nested column, `spectrum`.

### 3.4 SPE catalogs

Original schemas: [SPE output catalog data product](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dm10/spedpd/dpcards/spe_spepfoutputcatalog.html)

| Table | `_deep` | `_wide_survey` | `_wide` |
|---|---|---|---|
| `spectro_z_spe_quality` | 51,930,919 | 875,837,267 | 933,172,778 |
| `spectro_z_spe_classification` | 7,456,278 | 161,471,115 | 171,046,741 |
| `spectro_z_spe_galaxy_candidates` | 36,910,593 | 800,614,053 | 848,016,459 |
| `spectro_z_spe_star_candidates` | 13,622,726 | 277,528,395 | 293,737,969 |
| `spectro_z_spe_qso_candidates` | 36,063,712 | 735,013,630 | 778,703,735 |
| `spectro_line_features_spe_line_features` | 943,908,942 | 9,524,030,381 | 10,275,874,115 |
| `spectro_line_features_spe_atomic_indices` | 402,756,480 | 3,033,760,428 | 3,329,048,376 |
| `spectro_line_features_spe_molecular_indices` | 67,126,080 | 505,626,738 | 554,841,396 |
| `spectro_line_features_spe_continuum_features` | 22,375,360 | 168,542,246 | 184,947,132 |
| `spectro_model_spe_galaxy_models` | 36,936,097 | 800,701,396 | 848,103,802 |
| `spectro_model_spe_star_models` | 13,629,798 | 277,559,167 | 293,768,741 |
| `spectro_model_spe_qso_models` | 36,082,862 | 735,094,834 | 778,784,939 |

Unlike the MER and PHZ catalogs, most SPE catalogs contain multiple rows per object (identified by column `object_id`).
For example, the galaxy redshift candidate tables contain several redshift estimates per object, ranked by the column `spe_rank` (0 = best candidate).

+++

### 3.5 Catalog paths

Each catalog is stored in its own HATS collection.
A collection contains the catalog itself, plus ancillary data that some libraries (such as LSDB) use to make queries faster: a margin cache and, for some catalogs, an index table.
The collection for a catalog is at `individual_catalogs/hats/<catalog>/` and the catalog inside it is named `euclid_dr1_individual_catalogs_<catalog>-hats`.

The collections can be browsed [here]() [FIXME: add the irsa link].

In the rest of this tutorial we use the MER catalog for the Deep Survey, `mer_catalogue_deep_survey`.

```{code-cell} ipython3
# AWS S3 paths.
s3_bucket = "nasa-irsa-euclid-dr1"
# [FIXME] Confirm the S3 prefix once the catalogs are in the bucket.
collection_prefix = "FIXME/individual_catalogs/hats/mer_catalogue_deep_survey"
catalog_name = "euclid_dr1_individual_catalogs_mer_catalogue_deep_survey-hats"

collection_path = f"{s3_bucket}/{collection_prefix}"
dataset_path = f"{collection_path}/{catalog_name}/dataset"
schema_path = f"{dataset_path}/_common_metadata"

# S3 pointer. Use `anonymous=True` to access without credentials.
s3 = pyarrow.fs.S3FileSystem(anonymous=True)
```

## 4. Column Names and Added Columns

+++

First we'll load the Parquet schema (column information) of the catalog.
The Parquet schema is accessible from a few locations, all of which include the column names and types.
Here, we load it from the `_common_metadata` file because it also includes the column units and descriptions.

```{code-cell} ipython3
# Load the Parquet schema.
schema = pyarrow.parquet.read_schema(schema_path, filesystem=s3)

print(f"{len(schema)} columns in the {catalog_name} catalog")
```

The original tables' column names are mostly in all caps.
In the HATS catalogs and the catalogs available through IRSA's TAP service, all column names have been lower-cased.
In addition, the characters `-`, `.`, `(`, and `)` have been replaced with an underscore.
An example that includes both transformations is: `E(B-V)` -> `e_b_v_`.
Unlike the Q1 Merged Objects catalog, the table name is not prepended to the column names.

In the MER catalogs, three columns were renamed to match other IRSA services:

- `ra` : 'RIGHT_ASCENSION'
- `dec` : 'DECLINATION'
- `tileid` : 'TILE_INDEX', the ID of the Euclid tile the object was detected in.

Two columns were added to every catalog.
They contain the object's coordinates, calculated from its `object_id` following the [MER catalog description](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dm10/merdpd/dpcards/mer_finalcatalog.html#detailed-description-of-the-data-product).
They were used to partition the catalogs spatially, including the PHZ and SPE catalogs that don't have RA and Dec columns of their own.

- `ra_from_object_id` : RA coordinate (degrees) calculated from `object_id`.
- `dec_from_object_id` : Dec coordinate (degrees) calculated from `object_id`.

**HEALPix columns:**

These HEALPix indexes correspond to the `ra_from_object_id` and `dec_from_object_id` coordinates.
They are useful for spatial queries.

- `_healpix_9` : HEALPix order 9 pixel index.
  Order 9 pixels have a resolution (square root of area) of ~400 arcseconds or ~0.1 degrees.
- `_healpix_19` : HEALPix order 19 pixel index.
  Order 19 pixels have a resolution of ~0.4 arcseconds.
- `_healpix_29` : HEALPix order 29 pixel index.
  Order 29 pixels have a resolution of ~4e-4 arcseconds.

The HEALPix columns and the added coordinate columns appear first:

```{code-cell} ipython3
schema.names[:5]
```

**HATS columns:**

These are the HATS partitioning columns.
They appear in the Parquet file names but are not included inside the files.
However, PyArrow automatically makes them available as regular columns when the dataset is loaded as demonstrated in this tutorial.

- `Norder` : (hats column) HEALPix order at which the data is partitioned.
- `Npix` : (hats column) HEALPix pixel index at order Norder.
- `Dir` : (hats column) Integer equal to 10_000 * floor[Npix / 10_000].

The HATS columns appear at the end:

```{code-cell} ipython3
schema.names[-3:]
```

### 4.1 Find columns of interest

```{code-cell} ipython3
# Access the data type using the `field` method.
schema.field("flux_y_2fwhm_aper")
```

```{code-cell} ipython3
# The column metadata includes unit and description.
# Parquet metadata is always stored as bytestrings, which are denoted by a leading 'b'.
schema.field("flux_y_2fwhm_aper").metadata
```

Euclid DR1 offers many flux measurements, both from Euclid detections and from external ground-based surveys.
Flux column names start with "flux_", and columns associated with external surveys include "ext" in the name.

```{code-cell} ipython3
# Find all flux columns.
flux_columns = [name for name in schema.names if name.startswith("flux_")]
print(f"{len(flux_columns)} flux columns. First four are:")
flux_columns[:4]
```

```{code-cell} ipython3
external_flux_columns = [name for name in flux_columns if "ext" in name]
print(f"{len(external_flux_columns)} flux columns from external surveys. First four are:")
external_flux_columns[:4]
```

## 5. Basic Query with PyArrow

To demonstrate a basic query, we'll search for bright, point-like objects detected in the VIS mosaic.
PyArrow dataset filters are described at [Filtering by Expressions](https://arrow.apache.org/docs/python/compute.html#filtering-by-expressions), and the list of available functions is at [Compute Functions](https://arrow.apache.org/docs/python/api/compute.html).

```{code-cell} ipython3
dataset = pyarrow.dataset.dataset(dataset_path, partitioning="hive", filesystem=s3, schema=schema)

columns = ["object_id", "ra", "dec", "flux_vis_psf", "point_like_prob"]
bright_point_sources = (
    (pc.field("vis_det") == 1)  # detected in the VIS mosaic
    & (pc.field("spurious_flag") == 0)  # not flagged as spurious
    & (pc.field("point_like_prob") > 0.99)  # point-like
    & (pc.field("flux_vis_psf") > 1_000)  # brighter than 1 mJy in VIS
)

bright_objects = dataset.to_table(columns=columns, filter=bright_point_sources).to_pandas()
bright_objects
```

## 6. Cone Search with LSDB

[LSDB](https://docs.lsdb.io/) is a Python library built for HATS catalogs.
It is especially good for spatial searches and cross matches.
We can pass it the collection path and it will find the catalog inside it and use the ancillary data (such as the margin cache) when needed.
We prefix the path with `s3://` to indicate that it is an S3 path.

LSDB loads data lazily.
Here we open the catalog with only the columns we need and set up a cone search with a 2 arcminute radius around the planetary nebula NGC 6543, which is in the Euclid Deep Field North.

```{code-cell} ipython3
mer_deep = lsdb.open_catalog(
    f"s3://{collection_path}",
    columns=["object_id", "ra", "dec", "flux_vis_psf", "point_like_prob"],
)

ngc6543_ra, ngc6543_dec = 269.6394799, 66.6329765  # degrees
ngc6543_cone = mer_deep.cone_search(ra=ngc6543_ra, dec=ngc6543_dec, radius_arcsec=120)
ngc6543_cone
```

The cone search is not run until we call `compute()`, which loads the result into a NestedFrame (a pandas DataFrame with support for nested columns).
LSDB always loads the columns that the catalog was partitioned on, so `ra_from_object_id` and `dec_from_object_id` are included even though we didn't ask for them.
The cone search uses these columns.
LSDB always loads the columns that the catalog was partitioned on, so `ra_from_object_id` and `dec_from_object_id` are included even though we didn't ask for them.
The cone search uses these columns.

```{code-cell} ipython3
ngc6543_objects = ngc6543_cone.compute()
ngc6543_objects
```

## About this notebook

**Updated:** 2026-11-09

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html)