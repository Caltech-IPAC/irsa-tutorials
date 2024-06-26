---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: parquet
  language: python
  name: python3
---

# AllWISE Source Catalog Examples

+++

## Learning Goals

- Learn how to load data from the AllWISE parquet catalog that is partitioned by HEALPix pixel index at order 5.

- Learn efficient methods for performing common types of searches. This includes:

  - query/load using pandas, applying basic filters
  - query/load using pyarrow, applying advanced filters that combine and/or compare columns
  - perform nearest neighbor searches using pyarrow and astropy

+++

## Introduction

This notebook demonstrates access to the [HEALPix](https://ui.adsabs.harvard.edu/abs/2005ApJ...622..759G/abstract)-partitioned (order 5), [Apache Parquet](https://parquet.apache.org/) version of the [AllWISE Source Catalog](https://wise2.ipc.caltech.edu/docs/release/allwise/expsup/sec1_3.html#src_cat).
The catalog is available through the [AWS Open Data](https://aws.amazon.com/opendata) program, as part of the [NASA Open-Source Science Initiative](https://science.nasa.gov/open-science-overview).

Parquet is convenient for large astronomical catalogs in part because the storage format supports efficient database-style queries on the files themselves, without having to load the catalog into a database (or into memory) first.
The AllWISE catalog is fairly large at 340 GB.
The examples below demonstrate different methods that can be used to query the catalog, filtering the data and loading only the results.
Each method accesses the parquet files a bit differently and is useful for different types of queries.

+++

## Imports

```{code-cell} ipython3
!pip install -r requirements.txt
```

```{code-cell} ipython3
import sys

import hpgeom as hp
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib import colors
from matplotlib import pyplot as plt
from pyarrow.fs import S3FileSystem
```

## Setup catalog paths and query filters

+++

This AllWISE catalog is stored in an [AWS S3](https://aws.amazon.com/s3/) bucket.
To connect to an S3 bucket we just need to point the reader at S3 instead of the local filesystem, and pass in AWS credentials.
(Here, a "reader" is a python library that reads parquet files.)
We'll use [pyarrow.fs.S3FileSystem](https://arrow.apache.org/docs/python/generated/pyarrow.fs.S3FileSystem.html) for this because it is recognized by every reader in examples below, and we're already using pyarrow.
[s3fs](https://s3fs.readthedocs.io/en/latest/index.html) is another common option.
The call to `S3FileSystem` will look for AWS credentials in environment variables and/or the file ~/.aws/credentials.
Credentials can also be passed as keyword arguments.

```{code-cell} ipython3
bucket = "nasa-irsa-wise"
folder = "wise/allwise/catalogs/p3as_psd/healpix_k5"
parquet_root = f"{bucket}/{folder}/wise-allwise.parquet"

fs = S3FileSystem(region="us-west-2")  # the bucket is in region us-west-2
```

These limits will be used to query the catalog using specific filters created in examples below.
The Schema Access section (below) shows how to access column names and other schema information.

```{code-cell} ipython3
w1mpro_min = 10.0
ra_min, ra_max = 15, 25  # deg
dec_min, dec_max = 62, 72  # deg
polygon_corners = [(ra_min, dec_min), (ra_min, dec_max), (ra_max, dec_max), (ra_max, dec_min)]
radius = 5 * u.arcmin.to(u.deg)
```

The catalog is partitioned by HEALPix pixel index at order 5.
Queries can be most efficient when a filter on the partition column is included, since the reader can skip those partitions entirely.
Thus, for queries that include ra/dec constraints, we can usually speed up load times significantly by including a constraint on the HEALPix order 5 pixel index.
The required pixel indexes can be calculated using [hpgeom](https://hpgeom.readthedocs.io/en/latest/index.html), as demonstrated here.

```{code-cell} ipython3
order = 5  # the catalog is partitioned by HEALPix pixel index at order 5
nside = hp.order_to_nside(order)  # number of order 5 pixels along one side of a base pixel
```

```{code-cell} ipython3
# polygon search: get the set of pixel indexes that overlap the ra/dec polygon
polygon_pixels = hp.query_polygon(
    nside=nside,
    a=[corner[0] for corner in polygon_corners],  # ra values
    b=[corner[1] for corner in polygon_corners],  # dec values
    nest=True,  # catalog uses nested ordering scheme for pixel index
    inclusive=True,  # return all pixels that overlap with the polygon, and maybe a few more
)

print(f"polygon_pixels contains {len(polygon_pixels)} of a possible {hp.nside_to_npixel(nside)} pixels")
```

```{code-cell} ipython3
# cone search: get the set of pixel indexes that overlap a 5' circle around the ra/dec min
cone_pixels = hp.query_circle(
    nside=nside,
    a=ra_min,
    b=dec_min,
    radius=radius,
    nest=True,  # catalog uses nested ordering scheme for pixel index
    inclusive=True,  # return all pixels that overlap with the disk, and maybe a few more
)

# this can reduce the number of partitions the reader needs to look at from 12288 down to 2
print(f"cone_pixels contains {len(cone_pixels)} of a possible {hp.nside_to_npixel(nside)} pixels")
```

## Example 1:  Pandas with basic filters (magnitude limit and ra/dec polygon)

+++

Load using [pandas.read_parquet](https://pandas.pydata.org/docs/reference/api/pandas.read_parquet.html).
Filter for magnitudes above our w1mpro limit and a sky-area limited to the ra/dec polygon.

Pandas actually uses either pyarrow or fastparquet to interact with parquet files.
We'll choose pyarrow (the default).
For filter options, see the `filters` arg description in [ParquetDataset](https://arrow.apache.org/docs/python/generated/pyarrow.parquet.ParquetDataset.html#pyarrow.parquet.ParquetDataset).

```{code-cell} ipython3
%%time
# expect this to take 40-90 seconds
pandas_df = pd.read_parquet(
    parquet_root,
    filesystem=fs,
    # columns to be returned. similar to a SQL SELECT clause.
    columns=["designation", "ra", "dec", "w1mpro", "healpix_k5"],
    # row filters. similar to a SQL WHERE clause.
    # tuple conditions are joined by AND (for OR, use a list of lists)
    # supported operators: ==, !=, <, >, <=, >=, in, not in
    filters=[
        ("w1mpro", ">", w1mpro_min),
        ("ra", ">", ra_min),
        ("ra", "<", ra_max),
        ("dec", ">", dec_min),
        ("dec", "<", dec_max),
        # include filter on partition column for most efficient loading
        ("healpix_k5", "in", polygon_pixels),
    ],
)

pandas_df.describe()
```

## Example 2:  Pyarrow with advanced filters (color-color cuts for AGN)

+++

Load using [pyarrow.dataset.parquet_dataset](https://arrow.apache.org/docs/python/generated/pyarrow.dataset.parquet_dataset.html) and convert to pandas.
Useful for:

1. advanced filters that combine, compare, and/or create new columns. (this example)
2. speed. This method is more efficient at "discovering" the dataset than pandas. It also provides a persistent dataset object that can be reused for future queries, where pandas must re-discover the dataset every time. (this example and Example 3)

This example filters the catalog for AGN by making cuts in color-color space using the selection limits from [Mateos et al. (2012)](https://arxiv.org/pdf/1208.2530.pdf).
This is a more complicated filter than in Example 1 (it requires both constructing new columns and comparing values between columns) but this load is generally faster, demonstrating the efficiency of this method.

For basic info about the `columns` and `filter` arguments, see [Scanner](https://arrow.apache.org/docs/python/generated/pyarrow.dataset.Scanner.html#pyarrow.dataset.Scanner.from_dataset).
The construction of columns/filters is more involved than before because they must be passed as [Expressions](https://arrow.apache.org/docs/python/generated/pyarrow.dataset.Expression.html#pyarrow.dataset.Expression), and all operations must be done using pyarrow.compute [functions](https://arrow.apache.org/docs/python/api/compute.html) and `field`s.
This is demonstrated below.
Note that the catalog uses a file naming scheme called "hive", which the reader uses to identify partitions.
In other examples this is recognized automatically, but here we must pass it explicitly.

```{code-cell} ipython3
# define new columns for colors W1-W2 and W3-W4
w1w2 = pc.subtract(pc.field("w1mpro"), pc.field("w2mpro"))
w3w4 = pc.subtract(pc.field("w3mpro"), pc.field("w4mpro"))

# define the AGN locus, as in Mateos et al., 2012
locus = pc.multiply(pc.scalar(0.5), w3w4)
```

```{code-cell} ipython3
%%time
# expect this to take 20-60 seconds.
# notice this is generally faster than example 1 using pandas even though
# this filter is much more complicated, highlighting the efficiency of this method.

# load catalog metadata as a pyarrow dataset
pyarrow_ds = ds.parquet_dataset(f"{parquet_root}/_metadata", filesystem=fs, partitioning="hive")

# query for AGN using selection limits from Mateos et al., 2012
pyarrow_df = pyarrow_ds.to_table(
    # column filter. similar to a SQL SELECT clause.
    columns={
        "w1w2": w1w2,
        "w3w4": w3w4,
        "cntr": pc.field("cntr"),
        "ra": pc.field("ra"),
        "dec": pc.field("dec"),
        "healpix_k5": pc.field("healpix_k5"),
    },
    # row filter. similar to a SQL WHERE clause.
    filter=(
        # color-color cuts
        (w1w2 < pc.add(locus, pc.scalar(0.979)))
        & (w1w2 > pc.subtract(locus, pc.scalar(0.405)))
        & (w3w4 > pc.scalar(1.76))
        # to do an all-sky search, comment out the rest of the filter. expect it to take 30-60 min
        # same ra/dec polygon as before
        & (pc.field("ra") > ra_min)
        & (pc.field("ra") < ra_max)
        & (pc.field("dec") > dec_min)
        & (pc.field("dec") < dec_max)
        # same partition-column filter as before
        & (pc.field("healpix_k5").isin(polygon_pixels))
    ),
).to_pandas()
```

```{code-cell} ipython3
len(pyarrow_df.index)
```

```{code-cell} ipython3
colorbar_norm = colors.LogNorm(vmin=1, vmax=10)  # for an all-sky search, use vmax=100_000
pyarrow_df.plot.hexbin("w3w4", "w1w2", norm=colorbar_norm)
```

## Example 3:  Nearest-neighbor search (using pyarrow and astropy)

+++

Nearest-neighbor searches and cone searches generally use the on-sky separation distance to determine the matches.
It would be cumbersome to construct the new column and filter on it using the methods shown above because the separation distance is a fairly complicated function of ra and dec.
However, we can get pretty fast results by filtering down to the HEALPix pixels that cover the region, loading all the data in those partitions, and then using astropy to compute the separations and find the matches.

Here, we'll search for the 3 nearest neighbors to each of the 4 corners of our ra/dec polygon.
We'll load the data with pyarrow because because it makes this query significantly faster than with pandas (see explanation in Example 2).
(Note that astropy can also read parquet, but it can only read a single file at a time and so is less convenient.)
We'll use [astropy.coordinates.SkyCoord.match_to_catalog_sky](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord.match_to_catalog_sky) to do the actual nearest neighbor search.

```{code-cell} ipython3
# construct dictionary of pixels covering a cone around each polygon corner
# we did this once before but now we want all 4 corners
corner_cone_pixels = {
    (ra, dec): hp.query_circle(nside=nside, a=ra, b=dec, radius=radius, nest=True, inclusive=True)
    for (ra, dec) in polygon_corners
}
corner_cone_pixels
```

Find the 3 nearest neighbors of each corner:

```{code-cell} ipython3
%%time
# expect this to take 30-60 seconds
idcol = "cntr"
neighbor_ids = []  # store the allwise source id (cntr) of each neighbor
corner_neighbors = {}  # record neighbor info for each corner (for reference. not actually used)

# use same pyarrow dataset as before, but get it again to include the load time in this cell
pyarrow_ds = ds.parquet_dataset(f"{parquet_root}/_metadata", filesystem=fs, partitioning="hive")

for (ra, dec), cone_pixels in corner_cone_pixels.items():

    # load data from all pixels/partitions/files covering this corner
    src_tbl = pyarrow_ds.to_table(
        columns=[idcol, "ra", "dec"], filter=(pc.field("healpix_k5").isin(cone_pixels))
    )

    # get list of 3 nearest neighbors to this corner
    corner = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    allwise_sources = SkyCoord(ra=src_tbl["ra"] * u.degree, dec=src_tbl["dec"] * u.degree)
    neighbors = [corner.match_to_catalog_sky(allwise_sources, nthneighbor=i) for i in range(1, 4)]

    # get the allwise source ids. record/report some info
    corner_neighbors[(ra, dec)] = []
    for n, (idx, sep, _) in enumerate(neighbors):
        srcid = src_tbl[idcol][idx.item()].as_py()
        neighbor_ids.append(srcid)
        corner_neighbors[(ra, dec)].append((srcid, sep))
        print(f"neighbor {n+1}, corner ({ra}, {dec}): {srcid} at on-sky dist {sep.to('arcsec')}")

neighbor_ids
```

Load all the data (all columns) for all nearest neighbors:

```{code-cell} ipython3
%%time
# expect this to take 15-25 seconds
cone_pixels = [pixel for cone_pixels in corner_cone_pixels.values() for pixel in cone_pixels]
neighbors_df = pyarrow_ds.to_table(
    filter=((pc.field(idcol).isin(neighbor_ids)) & (pc.field("healpix_k5").isin(cone_pixels)))
).to_pandas()

neighbors_df
```

## Schema Access

+++

The schema can be viewed [online](http://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_1a.html) and also accessed from the parquet catalog itself.
The example below loads it using the `_common_metadata` file and [pyarrow.dataset.parquet_dataset](https://arrow.apache.org/docs/python/generated/pyarrow.dataset.parquet_dataset.html).
Note that the schema can also be accessed from the metadata in the footer of each Parquet file and from the `_metadata` file, but the method used here is generally faster and easier.
In addition, this `_common_metadata` file has extra information (units and descriptions) stored in the custom metadata of each column.

```{code-cell} ipython3
schema = ds.parquet_dataset(f"{parquet_root}/_common_metadata", filesystem=fs).schema
```

```{code-cell} ipython3
# access individual columns by name or index
fld = schema.field(1)  # equivalently: fld = schema.field("ra")

# basic column information
print(fld.name)
print(fld.type)
print(fld.nullable)

# units and descriptions are in the `metadata` attribute, always as bytestrings
print(fld.metadata)
print(fld.metadata[b"units"].decode())  # use decode to get a regular string
print(fld.metadata[b"description"].decode())
```

```{code-cell} ipython3
# full list of column names
schema.names
```

## About this notebook

This notebook was developed by Troy Raen (raen@ipac.caltech.edu) in conjunction with David Shupe, Jessica Krick and the IRSA Science Platform team at IPAC.
