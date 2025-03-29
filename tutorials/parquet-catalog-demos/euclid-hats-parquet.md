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

# Euclid Q1: MER Catalogs in HATS Parquet

+++

## Learning Goals
By the end of this tutorial, you will:

- Access basic metadata to understand the format and schema of this unified HATS Parquet dataset.
- Visualize the HATS partitioning of this dataset.
- Query this dataset for likely stars and create a color-magnitude diagram. (Recreate the figure from
  [Introduction to Euclid Q1 MER catalog](https://caltech-ipac.github.io/irsa-tutorials/tutorials/euclid_access/2_Euclid_intro_MER_catalog.html),
  this time with *all* likely stars.)

+++

## Introduction

+++

This notebook demonstrates accesses to a version of the
[Euclid Q1](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_q1.html) MER Catalogs
that is in Apache Parquet format, partitioned according to the
Hierarchical Adaptive Tiling Scheme (HATS), and stored in an AWS S3 bucket.

The catalog version accessed here is a single dataset which comprises all three MER Catalogs
-- MER, MER Morphology, and MER Cutouts -- which have been joined by Object ID.
Their schemas (pre-join) can be seen at
[Euclid Final Catalog description](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmq1/merdpd/dpcards/mer_finalcatalog.html).
Minor modifications were made to the parquet schema to accommodate the join (de-duplicating column names)
and for the HATS standard. These differences are shown below.

Parquet is a file format that enables flexible and efficient data access by, among other things,
supporting the application of both column and row filters when reading the data (very similar to a SQL query)
so that only the desired data is loaded into memory.

HATS is a spatial partitioning scheme based on HEALPix that aims to
produce partitions (files) of roughly equal size.
This makes the files more efficient to work with,
especially for large-scale analyses and/or parallel processing.
It does this by adapting the HEALPix order at which data is partitioned in a given catalog based
on the on-sky density of the rows it contains.
In other words, data from dense regions of sky will be partitioned at a higher order
(i.e., higher resolution; smaller pixel size) than data in sparse regions.
HATS-aware python packages are being developed to take full advantage of the partitioning.
In this notebook, we will use the [hats](https://hats.readthedocs.io/) library to visualize the
catalog and access the schema, and [lsdb](https://docs.lsdb.io/) to do a query for all likely stars.

+++

## Installs and imports

```{code-cell}
# # Uncomment the next line to install dependencies if needed.
# !pip install 'hats>=0.5' 'lsdb>=0.5' matplotlib 'numpy>=2.0' 'pyerfa>=2.0.1.3' s3fs
```

```{code-cell}
import os

import dask.distributed
import hats
import lsdb
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
```

```{tip}
If you run into an error that starts with,
"A module that was compiled using NumPy 1.x cannot be run in NumPy 2.1.3 as it may crash.",
make sure you have restarted the kernel since doing `pip install`. Then re-run the cell.
```

+++

## 1. Setup

```{code-cell}
# AWS S3 path where this dataset is stored.
s3_bucket = "irsa-fornax-testdata"
s3_key = "EUCLID/q1/mer_catalogue/hats"
euclid_s3_path = f"s3://{s3_bucket}/{s3_key}"

# Temporary try/except to handle credentials in different environments before public release.
try:
    # If running from within IPAC's network (maybe VPN'd in with "tunnel-all"),
    # your IP address acts as your credentials and this should just work.
    hats.read_hats(euclid_s3_path)
except FileNotFoundError:
    # If running from Fornax, credentials are provided automatically under the hood, but
    # hats ignores them in the call above and raises a FileNotFoundError.
    # Construct a UPath which will pick up the credentials.
    from upath import UPath

    euclid_s3_path = UPath(f"s3://{s3_bucket}/{s3_key}")
```

## 2. Visualize the on-sky density of Q1 Objects and HATS partitions

+++

Euclid Q1 covers four non-contiguous fields: Euclid Deep Field North (22.9 sq deg), Euclid Deep Field Fornax (12.1 sq deg), Euclid Deep Field South (28.1 sq deg), and LDN1641.
We can visualize the Object density in the four fields using `hats`.

```{code-cell}
# Load the dataset.
euclid_hats = hats.read_hats(euclid_s3_path)

# Visualize the on-sky distribution of objects in the Q1 MER Catalog.
hats.inspection.plot_density(euclid_hats)
```

We can see how the on-sky density maps to the HATS partitions by calling `plot_pixels`.

```{code-cell}
# Visualize the HEALPix orders of the dataset partitions.
hats.inspection.plot_pixels(euclid_hats)
```

## 3. CMD of stars in Euclid Q1

+++

In this section, we query the Euclid Q1 MER catalogs for likely stars and create a color-magnitude diagram (CMD), following
[Introduction to Euclid Q1 MER catalog](https://caltech-ipac.github.io/irsa-tutorials/tutorials/euclid_access/2_Euclid_intro_MER_catalog.html).
Here, we use `lsdb` to query the parquet files that are sitting in an S3 bucket (the intro notebook uses `pyvo` to query the TAP service).
`lsdb` enables efficient, large-scale queries on HATS catalogs, so let's look at *all* likely stars in Euclid Q1 instead of limiting to 10,000.

`lsdb` uses Dask for parallelization. So first, set up the workers.

```{code-cell}
client = dask.distributed.Client(
    n_workers=os.cpu_count(), threads_per_worker=2, memory_limit="auto"
)
```

The data will be lazy-loaded. This means that commands like `query` are not executed until the data is actually required.

```{code-cell}
# Load the dataset.
columns = [
    "TILEID",
    "FLUX_VIS_PSF",
    "FLUX_Y_TEMPLFIT",
    "FLUX_J_TEMPLFIT",
    "FLUX_H_TEMPLFIT",
    "POINT_LIKE_FLAG",
]
euclid_lsdb = lsdb.read_hats(euclid_s3_path, columns=columns)

# Set up the query for likely stars.
star_cuts = "FLUX_VIS_PSF > 0 & FLUX_Y_TEMPLFIT > 0 & FLUX_J_TEMPLFIT > 0 & FLUX_H_TEMPLFIT > 0 & POINT_LIKE_FLAG == 1"
euclid_stars = euclid_lsdb.query(star_cuts)
```

```{code-cell}
# Peek at the data.
euclid_stars.head(10)
```

We peeked at the data but we haven't loaded all of it yet.
What we really need in order to create a CMD is the magnitudes, so let's calculate those now.
Appending `.compute()` to the commands will trigger Dask to actually load this data into memory.
It is not strictly necessary, but will allow us to look at the data repeatedly without having to re-load it each time.

```{code-cell}
# Calculate magnitudes. Appending `.compute()` triggers Dask to load this data now.
mag_y = (-2.5 * np.log10(euclid_stars["FLUX_Y_TEMPLFIT"]) + 23.9).compute()
mag_h = (-2.5 * np.log10(euclid_stars["FLUX_H_TEMPLFIT"]) + 23.9).compute()

print(f"Loaded magnitudes of {len(mag_y):,} likely stars in Euclid Q1.")
```

Create the CMD

```{code-cell}
hb = plt.hexbin(mag_y - mag_h, mag_y, norm=matplotlib.colors.LogNorm(vmin=1, vmax=50_000))
plt.colorbar(hb)
plt.xlabel("Y-H")
plt.ylabel("Y")
plt.xlim(-10, 10)
plt.ylim(10, 35)
plt.title("Stars in Euclid Q1 MER Catalog")
plt.show()
```

```{code-cell}
# Close the Dask client.
client.close()
```

## 4. Schema

+++

The three catalogs MER, MER Morphology, and MER Cutouts have been joined together in this parquet version.

IRSA's
[Cloud Access](https://caltech-ipac.github.io/irsa-tutorials/tutorials/cloud_access/cloud-access-intro.html#navigate-a-catalog-and-perform-a-basic-query)
notebook shows how to work with parquet schemas.
`hats` will return the same pyarrow schema object shown in that notebook, so let's use it.

```{code-cell}
# Fetch the pyarrow schema from hats.
euclid_hats = hats.read_hats(euclid_s3_path)
schema = euclid_hats.schema

print(f"{len(schema)} columns in the combined Euclid Q1 MER Catalogs")
```

Two columns have been added to the top of the schema.
'_healpix_29' is the pixel index at HEALPix order 29.
It is used by `hats` and is generally useful for spatial queries.
'TILEID' is the Euclid MER tile index, added for convenience.

```{code-cell}
schema.names[:5]
```

Three columns have been added to the bottom: 'Norder', 'Dir', and 'Npix'. These are the HATS partitioning columns.

```{code-cell}
schema.names[-5:]
```

In the original schemas, a handful of column names appear in both the main MER catalog and one of its secondaries (MER Cutouts).
To ensure that column names are unique in the parquet, the name of the secondary catalog was appended
to the affected column names, separated by "-".

```{code-cell}
# Find the column names that were affected in the secondary catalogs.
mer_secondaries = ["MORPH", "CUTOUTS"]
[name for name in euclid_hats.schema.names if name.split("-")[-1] in mer_secondaries]
```

For each of these, there is another column without the catalog name appended, which belongs to the main catalog.

```{code-cell}
print("-- MER column --")
print(schema.field("RIGHT_ASCENSION"))
print(schema.field("RIGHT_ASCENSION").metadata)

print("-- MER Cutouts column --")
print(schema.field("RIGHT_ASCENSION-CUTOUTS"))
print(schema.field("RIGHT_ASCENSION-CUTOUTS").metadata)
```

## About this notebook

**Authors:** Troy Raen (Developer; Caltech/IPAC-IRSA) and the IRSA Data Science Team.

**Updated:** 2025-03-29

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.
