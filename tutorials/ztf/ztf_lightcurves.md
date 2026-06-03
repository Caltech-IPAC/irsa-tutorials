---
authors:
- name: Jaladh Singhal
- name: Troy Raen
- name: "Brigitta Sip\u0151cz"
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.3
kernelspec:
  display_name: irsa-tutorials
  language: python
  name: python3
---

(ztf-lightcurves)=
# Access ZTF DR24 Light Curves from HATS Catalog

+++

## Learning Goals

By the end of this tutorial, you will learn how to:

- Open ZTF DR24 HATS catalogs for light curves and the Objects Table using `lsdb`.
- Retrieve light curves for specific objects by ZTF object IDs using an index search.
- Retrieve light curves for objects in a sky region using a cone search.
- Cross-reference the Objects Table to enrich cone search results with per-object variability statistics.
- Plot ZTF light curves (filtered by variability statistics).

+++

## Introduction

The ZTF DR24 enhanced data products at IRSA include two [HATS](https://irsa.ipac.caltech.edu/docs/parquet_catalogs/#hats) (Hierarchical Adaptive Tiling Scheme) catalogs hosted on AWS S3:

- **Lightcurves catalog**: one row per ZTF object, with a nested column storing the full photometry time series — timestamps, magnitudes, uncertainties, and quality flags.
- **Objects Table**: one row per ZTF object, with collapsed light curve metrics such as magnitude RMS, chi-squared variability statistic, number of good observations, and mean magnitude.

These HATS catalogs offer a scalable, cloud-native alternative to the [ZTF light curve service](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html), enabling efficient access especially when the service is overloaded.
The [lsdb](https://docs.lsdb.io/en/latest/index.html) Python library provides a convenient interface for working with HATS catalogs, including spatial queries and object-ID-based lookups.

This tutorial covers two common entry points for accessing ZTF light curves:

1. **Object IDs**: you have specific ZTF object IDs — from a previous query, a catalog crossmatch, or a published object list — and want their light curves directly.
2. **RA/Dec**: you have sky coordinates and want all ZTF objects within a given radius.

Both approaches are demonstrated below. An optional section then shows how to join the position search results with the Objects Table to select and plot the most variable objects using robust variability statistics.

For more context on ZTF DR24 data products, refer to the [ZTF DR24 release notes](https://irsa.ipac.caltech.edu/data/ZTF/docs/releases/ztf_release_notes_latest) and [explanatory supplement](https://irsa.ipac.caltech.edu/data/ZTF/docs/ztf_explanatory_supplement.pdf) at IRSA.

+++

## Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install s3fs "lsdb>=0.6.6,<0.8" pyarrow pandas astropy matplotlib
```

```{code-cell} ipython3
import s3fs
import lsdb
import pyarrow.parquet as pq
from astropy.coordinates import SkyCoord
import pandas as pd
from astropy import units as u
import os
import matplotlib.pyplot as plt
from dask.distributed import Client
```

```{code-cell} ipython3
pd.set_option("display.max_colwidth", None)
pd.set_option("display.min_rows", 18)
```

## 1. Locate ZTF DR24 HATS Catalogs in the Cloud

From IRSA's [cloud data access page](https://irsa.ipac.caltech.edu/cloud_access/), we identify the S3 bucket and path prefixes for the ZTF DR24 HATS catalogs:

```{code-cell} ipython3
ztf_bucket = "ipac-irsa-ztf"
ztf_lc_hats_prefix = "ztf/enhanced/dr24/lc/hats" # Lightcurves catalog
ztf_objects_hats_prefix = "ztf/enhanced/dr24/objects/hats" # Objects Table
```

[s3fs](https://s3fs.readthedocs.io/en/latest/) provides a filesystem-like Python interface for AWS S3 buckets.
First, we create an S3 client:

```{code-cell} ipython3
s3 = s3fs.S3FileSystem(anon=True)
```

Let's list the contents of the ZTF DR24 Lightcurves HATS **Collection**:

```{code-cell} ipython3
s3.ls(f"{ztf_bucket}/{ztf_lc_hats_prefix}")
```

In this collection, you can see collection properties, catalog, index table, and margin cache in order.
You can explore more directories to see how this HATS Collection follows the directory structure described in IRSA's documentation on [HATS partitioning and HATS Collections](https://irsa.ipac.caltech.edu/docs/parquet_catalogs/#hats).

As per the documentation, the Parquet file containing the schema for this catalog is stored in `dataset/_common_metadata`.
Let's save its path for later use (using the catalog name identified from the listing above):

```{code-cell} ipython3
ztf_lc_schema_path = "ztf_dr24_lc-hats/dataset/_common_metadata"  # ztf_dr24_lc-hats is the catalog name identified above
```

Similarly, let's list the ZTF DR24 Objects Table HATS Collection:

```{code-cell} ipython3
s3.ls(f"{ztf_bucket}/{ztf_objects_hats_prefix}")
```

```{code-cell} ipython3
ztf_objects_schema_path = "ztf_dr24_objects-hats/dataset/_common_metadata"  # ztf_dr24_objects-hats is the catalog name identified above
```

## 2. Explore the Catalog Schema

Before querying the catalogs, let's inspect what columns are available in each.
We read schemas from the `_common_metadata` files, which also contain column metadata such as units and descriptions:

```{code-cell} ipython3
def pq_schema_to_df(schema):
    """Convert a PyArrow schema to a Pandas DataFrame."""
    return pd.DataFrame(
        [
            (
                field.name,
                str(field.type),
                field.metadata.get(b"unit", b"").decode(),
                field.metadata.get(b"description", b"").decode()
            )
            for field in schema
        ],
        columns=["name", "type", "unit", "description"]
    )
```

```{code-cell} ipython3
ztf_lc_schema = pq.read_schema(
    f"s3://{ztf_bucket}/{ztf_lc_hats_prefix}/{ztf_lc_schema_path}",
    filesystem=s3
)
ztf_lc_schema_df = pq_schema_to_df(ztf_lc_schema)
ztf_lc_schema_df
```

Notice the `lightcurve` column — this is a **[nested](https://nested-pandas.readthedocs.io/) column** that stores the full photometric time series for each ZTF object.
Each element of `lightcurve` is itself a table with columns including `hmjd`, `mag`,`magerr`, `clrcoeff` and `catflags`.
We save the list of columns interesting to us for later use when opening the catalog with `lsdb`:

```{code-cell} ipython3
ztf_lc_columns = ["objectid", "objra", "objdec", "filterid", "nepochs", "lightcurve"]
```

## 3. Get Light Curves by Object ID

If you have specific ZTF object IDs, you can retrieve their light curves directly using an index search — no spatial filter needed.
This is the fastest approach for targeted lookups.

### 3.1 Open the Lightcurves Catalog

We open the ZTF DR24 Lightcurves HATS catalog. No data is read yet — lsdb opens catalogs [lazily](https://docs.lsdb.io/en/latest/tutorials/lazy_operations.html):

```{code-cell} ipython3
ztf_lc_catalog = lsdb.open_catalog(
    f"s3://{ztf_bucket}/{ztf_lc_hats_prefix}",
    columns=ztf_lc_columns
)
ztf_lc_catalog
```

### 3.2 Identify the Index Column

The ZTF DR24 Lightcurves HATS catalog ships with an ancillary index table that enables fast lookups by object ID.
Let's identify which column is indexed:

```{code-cell} ipython3
ztf_lc_idx_column = list(ztf_lc_catalog.hc_collection.all_indexes.keys())[0]
print(f"Index column: {ztf_lc_idx_column}")
```

### 3.3 Perform an Index Search

We use the same object IDs from the [ZTF light curve API docs](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html) multi-object example — you can compare results from this tutorial directly with that service.
In your workflow, these IDs might come from a previous query, a catalog crossmatch, or a published object catalog:

```{code-cell} ipython3
object_ids = [686103400034440, 686103400106565]
```

```{code-cell} ipython3
ztf_lcs_by_id = ztf_lc_catalog.id_search(values={ztf_lc_idx_column: object_ids})
ztf_lcs_by_id
```

### 3.4 Compute and Inspect the Results

Now we execute the query by calling `compute()`, which reads the data into memory as a Pandas DataFrame. We use a Dask client to parallelize across partitions, manage 
memory, and monitor progress in the Dask dashboard.

```{code-cell} ipython3
def get_nworkers(catalog):
    return min(os.cpu_count(), catalog.npartitions + 1)

with Client(n_workers=get_nworkers(ztf_lcs_by_id),
            threads_per_worker=1,
            memory_limit=None  # each partition can be several GB; avoid per-worker cap
            ) as client:
    print(f"You can monitor progress in the Dask dashboard at {client.dashboard_link}")
    ztf_lcs_by_id_df = ztf_lcs_by_id.compute()

ztf_lcs_by_id_df
```

```{code-cell} ipython3
print(f"Found {len(ztf_lcs_by_id_df)} light curves for {len(object_ids)} objects.")
```

Each row is one ZTF object. The `lightcurve` column contains a nested DataFrame per object.
Let's inspect the light curve of the first object:

```{code-cell} ipython3
ztf_lcs_by_id_df['lightcurve'].iloc[0]
```

### 3.5 Plot Light Curves
When plotting the light curves, it's important to note that we apply `catflags == 0` filter to keep only clean epochs (as described in the [explanatory supplement](https://irsa.ipac.caltech.edu/data/ZTF/docs/ztf_explanatory_supplement.pdf) section 13.6).

```{code-cell} ipython3
fig, axs = plt.subplots(len(ztf_lcs_by_id_df), 1,
                        figsize=(10, 4 * len(ztf_lcs_by_id_df)),
                        constrained_layout=True)

if len(ztf_lcs_by_id_df) == 1:
    axs = [axs]

for ax, (_, row) in zip(axs, ztf_lcs_by_id_df.iterrows()):
    lc = row['lightcurve'].query("catflags == 0")  # to keep only clean epochs
    title = f"ZTF Object {row['objectid']}  (RA={row['objra']:.4f}°, Dec={row['objdec']:.4f}°)"
    pts = ax.plot(lc['hmjd'], lc['mag'], '.', markersize=4, zorder=3)
    ax.errorbar(
        lc['hmjd'], lc['mag'], yerr=lc['magerr'],
        fmt='none', ecolor=pts[0].get_color(), elinewidth=0.8, alpha=0.3, zorder=2
    )
    ax.set_ylabel("Magnitude")
    ax.set_xlabel("HMJD")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=10)

fig.suptitle("ZTF DR24 Light Curves — Object ID Search Results", fontsize=13, y=1.02)
plt.show()
```

## 4. Get Light Curves by Sky Position

If you have sky coordinates and want all ZTF objects within a given area, use a cone search.

```{important} ZTF objects are defined per (filter, field, quadrant)
ZTF objects (i.e., unique object IDs) are defined _per_ (filter, field, quadrant).
This means that observations of a single _astrophysical_ object are usually spread out amongst several different _ZTF_ objects.

At minimum, a given astrophysical object will be represented by up to 3 ZTF objects, one per filter (g, r, and i).
The per-filter observations may themselves be separated into additional ZTF objects if the astrophysical object lies near the boundary of a ZTF field and/or quadrant.

ZTF's pixel scale is 1"/pixel (see [ZTF Technical Specifications](https://www.ptf.caltech.edu/page/ztf_technical)), so combining all ZTF objects within a 1" cone search may be reasonable for a given astrophysical object.
```

### 4.1 Define a Spatial Filter

We use the same sky position as the [ZTF light curve API docs](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html) positional example but you can specify any coordinates and search radius you want:

```{code-cell} ipython3
target = SkyCoord(ra=298.0025, dec=29.87147, unit="deg")
search_radius = 1 * u.arcsec # to keep the runtime minimal for this tutorial
```

Using lsdb, we define a cone [search object](https://docs.lsdb.io/en/latest/tutorials/region_selection.html#4.-The-Search-object) for this region:

```{code-cell} ipython3
spatial_filter = lsdb.ConeSearch(
    ra=target.ra.deg,
    dec=target.dec.deg,
    radius_arcsec=search_radius.to(u.arcsec).value
)
```

### 4.2 Define Row Filters

In addition to the spatial filter, we can pre-filter rows using Parquet column statistics.
Here we keep only objects with more than 25 epochs, focusing on well-sampled light curves:

```{code-cell} ipython3
row_filters = [
    ["nepochs", ">", 25],
    # additional filters can be added here if desired
    ]
```

### 4.3 Open the Filtered Lightcurves Catalog

We open the catalog with both filters applied. lsdb evaluates this lazily — no data is read yet:

```{code-cell} ipython3
ztf_lc_cone = lsdb.open_catalog(
    f"s3://{ztf_bucket}/{ztf_lc_hats_prefix}",
    search_filter=spatial_filter,
    columns=ztf_lc_columns,
    filters=row_filters
)
ztf_lc_cone
```

Notice that only the partitions overlapping the cone are included, avoiding reads of the full catalog.

### 4.4 Compute and Inspect the Results

Now we execute the query by calling `compute()`, which reads the data into memory as a Pandas DataFrame. We use a Dask client to parallelize across partitions, manage 
memory, and monitor progress in the Dask dashboard.

```{code-cell} ipython3
with Client(n_workers=get_nworkers(ztf_lc_cone),
            threads_per_worker=1,
            memory_limit=None  # each partition can be several GB; avoid per-worker cap
            ) as client:
    print(f"You can monitor progress in the Dask dashboard at {client.dashboard_link}")
    ztf_lc_cone_df = ztf_lc_cone.compute()
```

```{code-cell} ipython3
print(f"Found {len(ztf_lc_cone_df)} ZTF light curves for the search criteria.")
```

```{code-cell} ipython3
ztf_lc_cone_df.head(5)
```

Each row corresponds to one ZTF object. The `lightcurve` column contains a nested DataFrame per object:

```{code-cell} ipython3
ztf_lc_cone_df['lightcurve'].iloc[0]
```

## 5. [Optional] Look Up Additional Info from the Objects Table

```{note}
This section is optional — skip it if you don't need additional information beyond the raw light curves from section 4.
```

### 5.1 Explore the Objects Table Schema

The Objects Table contains summary statistics for each ZTF object.
Let's inspect its schema to identify columns of interest:

```{code-cell} ipython3
ztf_objects_schema = pq.read_schema(
    f"s3://{ztf_bucket}/{ztf_objects_hats_prefix}/{ztf_objects_schema_path}",
    filesystem=s3
)
pq_schema_to_df(ztf_objects_schema)
```

```{important} objectid == oid
ZTF's object ID column is named `objectid` in Lightcurves and `oid` in Objects Table.
Despite this difference, the two columns are the same and can be used to join the catalogs. 
```

We'll select a subset of columns useful for characterizing and annotating variable objects:

```{code-cell} ipython3
ztf_objects_columns = ['oid', 'ra', 'dec', 'filtercode', 'ngoodobsrel', 'chisq', 'magrms', 'meanmag', 'medianabsdev']
```

### 5.2 Open the Objects Table

We reuse the same `spatial_filter` from section 4 to retrieve Objects Table entries for the same sky region. This is important for ensuring we only retrieve rows relevant to the light curves we got from the position search.

```{code-cell} ipython3
ztf_objects_cone = lsdb.open_catalog(
    f"s3://{ztf_bucket}/{ztf_objects_hats_prefix}",
    search_filter=spatial_filter,
    columns=ztf_objects_columns
)
ztf_objects_cone
```

### 5.3 Compute and Inspect

```{code-cell} ipython3
with Client(n_workers=get_nworkers(ztf_objects_cone),
            threads_per_worker=1,
            memory_limit=None) as client:
    ztf_objects_cone_df = ztf_objects_cone.compute()
```

```{code-cell} ipython3
ztf_objects_cone_df
```

### 5.4 Merge Objects Table Info into Light Curves

We merge the Objects Table with the position search light curves on the shared object ID via an inner join:

```{code-cell} ipython3
combined_df = ztf_lc_cone_df.merge(
    ztf_objects_cone_df,
    left_on='objectid',
    right_on='oid',
    how='inner'
)
combined_df
```

## 6. Plot Most Variable Light Curves from the Position Search

Using the `chisq` column, we rudimentarily select the top 2 most variable objects from the position search results combined with the Objects Table.

```{code-cell} ipython3
# most_variable = ztf_lc_cone_df.iloc[:2]  # uncomment if you skipped section 5, and comment the line below
most_variable = combined_df.nlargest(2, 'chisq')
most_variable
```

Then we plot their light curves annotated with summary statistics:

```{code-cell} ipython3
fig, axs = plt.subplots(len(most_variable), 1,
                        figsize=(10, 4 * len(most_variable)),
                        constrained_layout=True)

if len(most_variable) == 1:
    axs = [axs]

for ax, (_, row) in zip(axs, most_variable.iterrows()):
    lc = row['lightcurve'].query("catflags == 0")  # to keep only clean epochs
    title = (f"ZTF Object {row['objectid']}  ({row['filtercode']} filter)\n"
             f"χ²={row['chisq']:.2f},  RMS mag={row['magrms']:.4f},  "
             f"mean mag={row['meanmag']:.3f},  N good obs={int(row['ngoodobsrel'])}")
    pts = ax.plot(lc['hmjd'], lc['mag'], '.', markersize=4, zorder=3)
    ax.errorbar(
        lc['hmjd'], lc['mag'], yerr=lc['magerr'],
        fmt='none', ecolor=pts[0].get_color(), elinewidth=0.8, alpha=0.3, zorder=2
    )
    ax.set_ylabel("Magnitude")
    ax.set_xlabel("HMJD")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=10)

fig.suptitle("Most Variable ZTF DR24 Objects from Position Search (annotated with Objects Table data)", fontsize=13, y=1.02)
plt.show()
```

## About this notebook

**Updated:** 2026-06-02

**Contact:** the [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or to report problems.

**AI Acknowledgement:** This tutorial was developed with the assistance of AI tools.
