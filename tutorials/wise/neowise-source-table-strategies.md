---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: science_demo
  language: python
  name: python3
---

# Strategies to Efficiently Work with NEOWISE Single-exposure Source Table in Parquet

+++

An executed version of this notebook can be seen on
[IRSA's website](https://irsa.ipac.caltech.edu/docs/notebooks/neowise-source-table-strategies.html).

+++


This notebook discusses strategies for working with the Apache Parquet version of the
[NEOWISE](https://irsa.ipac.caltech.edu/Missions/wise.html) Single-exposure Source Table
and provides the basic code needed for each approach.
This is a very large catalog -- 11 years and 42 terabytes in total with 145 columns and 200 billion rows.
Most of the work shown in this notebook is how to efficiently deal with so much data.

Learning Goals:

- Identify use cases that will benefit most from using this Parquet format.
- Understand how this dataset is organized and three different ways to efficiently
  slice it in order to obtain a subset small enough to load into memory and process.
- Feel prepared to apply these methods to a new use case and determine efficient filtering and slicing options.

+++

## 1. Introduction

+++

The NEOWISE Single-exposure Source Table comprises 11 years of data.
Each year on its own would be considered "large" compared to astronomy catalogs produced
contemporaneously, so working with the full dataset requires extra consideration.
In this Parquet version, each year is stored as an independent Parquet dataset.
This tutorial shows how to combine them and work with all years as one.
The data are partitioned by HEALPix ([Górski et al., 2005](https://doi.org/10.1086/427976)) order k=5.
HEALPix is a tessellation of the sky, so partitioning the dataset this way makes it especially
efficient for spatial queries.
In addition, the access methods shown below are expected to perform well when parallelized.

The terms "partition", "filter", and "slice" all relate to subsets of data.
In this notebook we'll use them with definitions that overlap but are not identical, as follows.
A "partition" includes all data observed in a single HEALPix pixel.
This data is grouped together in the Parquet files.
A "filter" is a set of criteria defined by the user and applied when reading data so that only the desired rows are loaded.
We'll use it exclusively to refer to a PyArrow `filter`.
The criteria can include partitions and/or any other column in the dataset.
A "slice" is a generic subset of data.
There are a few ways to obtain a slice; one is by applying a filter.

This notebook is expected to require about 2 CPUs and 50G RAM and to complete in about 10 minutes.
These estimates are based on testing in science platform environments.
Your numbers will vary based on many factors including compute power, bandwidth, and physical distance from the data.
The required RAM and runtime can be reduced by limiting the number of NEOWISE years loaded.

+++

### 1.1 When to use the Parquet version

IRSA provides large catalogs in file formats (e.g., Parquet) primarily to support use cases that require bulk access.
The Parquet version of the NEOWISE Source Table, coupled with the methods demonstrated in this tutorial,
are expected to be the most efficient option for large-ish use cases like classifying, clustering, or
building light curves for many objects.
In general, the efficiency will increase with the number of rows required from each slice because the slice
only needs to be searched once regardless of the number of rows that are actually loaded.
In addition, this access route tends to perform well when parallelized.
Note that these use cases (including this notebook) are often too large for a laptop and may perform
poorly and/or crash if attempted.

For small-ish use cases like searching for a handful of objects, other access routes like
PyVO and TAP queries will be faster.

Consider using this tutorial if either of the following is true:

- Your use case is large enough that you are considering parallelizing your code to speed it up.
- Your sample size is large enough that loading the data using a different method is likely to
  take hours, days, or longer.

+++

### 1.2 Recommended approach

The basic process is:

1. Load the catalog metadata as a PyArrow dataset.
2. Decide how to slice the dataset (e.g., by year, partition, and/or file) depending on your use case.
3. Iterate and/or parallelize over the slices. For each slice:
    1. Use the PyArrow dataset to load data of interest, applying row filters during the read.
    2. Process the data as you like (e.g., cluster, classify, etc.).
    3. Write your results to disk and/or return them for further processing.
4. (Optional) Concatenate your results and continue processing.

This notebook covers steps 1 through 3.1 and indicates where to insert your own code to proceed with 3.2.
Here we iterate over slices, but the same code can be parallelized using any multi-processing framework.
A fully-worked example is shown in the light curve notebook linked below.

+++

### 1.3 See also

- [](#cloud-access-intro)
- [](#allwise-source-catalog-parquet)
- [](#neowise-lightcurves-parquet)

+++

## 2. Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install hpgeom pandas pyarrow
```

```{code-cell} ipython3
import re  # parse strings
import sys  # check size of loaded data

import hpgeom  # HEALPix math
import pandas as pd  # store and manipulate table data
import pyarrow.compute  # construct dataset filters
import pyarrow.dataset  # load and query the NEOWISE dataset
import pyarrow.fs  # interact with the S3 bucket storing the NEOWISE catalog
```

## 3. Setup

### 3.1 Define variables and helper functions

+++

Choose which NEOWISE years to include.
Expect the notebook to require about 4G RAM and 1 minute of runtime per year.

```{code-cell} ipython3
# All NEOWISE years => about 40G RAM and 10 minutes runtime
YEARS = [f"year{yr}" for yr in range(1, 12)] + ["addendum"]

# To reduce the needed RAM or runtime, uncomment the next line and choose your own years.
# Years 1 and 8 are needed for the median_file and biggest_file (defined below).
# YEARS = [1, 8]
```

Column and partition variables:

```{code-cell} ipython3
# subset of columns to load
flux_columns = ["w1flux", "w1sigflux", "w2flux", "w2sigflux"]
COLUMN_SUBSET = ["cntr", "source_id", "ra", "dec"] + flux_columns

# partitioning info. do not change these values.
K = 5  # healpix order of the catalog partitioning
KCOLUMN = "healpix_k5"  # partitioning column name
KFIELD = pyarrow.compute.field(KCOLUMN)  # pyarrow compute field, to be used in filters
```

Paths:

```{code-cell} ipython3
# We're going to look at several different files, so make a function to return the path.
def neowise_path(year, file="_metadata"):
    """Return the path to a file. Default is "_metadata" file of the given year's dataset.

    Parameters
    ----------
    year : int
        NEOWISE year for which the path is being generated.
    file : str
        The name of the file to the returned path.

    Returns
    -------
    str
        The path to the file.
    """
    # This information can be found at https://irsa.ipac.caltech.edu/cloud_access/.
    bucket = "nasa-irsa-wise"
    base_prefix = "wise/neowiser/catalogs/p1bs_psd/healpix_k5"
    root_dir = f"{bucket}/{base_prefix}/{year}/neowiser-healpix_k5-{year}.parquet"
    return f"{root_dir}/{file}"
```

Some representative partitions and files (see dataset stats in the Appendix for how we determine these values):

```{code-cell} ipython3
# pixel index of the median partition and the biggest partition by number of rows
median_part = 11_831
biggest_part = 8_277

# path to the median file and the biggest file by file size on disk (see Appendix)
median_file = neowise_path("year8", "healpix_k0=1/healpix_k5=1986/part0.snappy.parquet")
biggest_file = neowise_path("year1", "healpix_k0=2/healpix_k5=2551/part0.snappy.parquet")
```

Convenience function for displaying a table size:

```{code-cell} ipython3
# We'll use this function throughout the notebook to see how big different tables are.
def print_table_size(table, pixel_index=None):
    """Prints the shape (rows x columns) and size (GiB) of the given table.

    Parameters
    ----------
    table : pyarrow.Table
        The table for which to print the size.
    pixel_index : int or str or None
        The pixel index corresponding to the partition this table was loaded from.
    """
    if pixel_index is not None:
        print(f"pixel index: {pixel_index}")
    print(f"table shape: {table.num_rows:,} rows x {table.num_columns} columns")
    print(f"table size: {sys.getsizeof(table) / 1024**3:.2f} GiB")
```

### 3.2 Load NEOWISE metadata as a pyarrow dataset

+++

The metadata contains column names, schema, and row-group statistics for every file in the dataset.
Later, we will use this pyarrow dataset object to slice and query the catalog in several different ways.

```{code-cell} ipython3
# This catalog is so big that even the metadata is big.
# Expect this cell to take about 30 seconds per year.
fs = pyarrow.fs.S3FileSystem(region="us-west-2", anonymous=True)

# list of datasets, one per year
year_datasets = [
    pyarrow.dataset.parquet_dataset(neowise_path(yr), filesystem=fs, partitioning="hive")
    for yr in YEARS
]

# unified dataset, all years
neowise_ds = pyarrow.dataset.dataset(year_datasets)
```

`neowise_ds` is a [UnionDataset](https://arrow.apache.org/docs/python/generated/pyarrow.dataset.UnionDataset.html).
All methods demonstrated for pyarrow datasets in the AllWISE demo notebook can be used with
`neowise_ds` and will be applied to all years as if they were a single dataset.
In addition, a separate [Dataset](https://arrow.apache.org/docs/python/generated/pyarrow.dataset.Dataset.html)
for each year is stored in the list attribute `neowise_ds.children` (== `year_datasets`, loaded above),
and the same methods can be applied to them individually.

+++

## 4. Example: Slice by partition

+++

This example shows how to load data from each partition separately.
The actual "slicing" is done by applying a filter to the pyarrow dataset `neowise_ds`.
Constructing filters was discussed in the AllWISE notebook linked above.

```{code-cell} ipython3
# number of order K pixels covering the full sky
npixels = hpgeom.nside_to_npixel(hpgeom.order_to_nside(order=K))

# iterate over all partitions
for pix in range(npixels):

    # slice and load to get all rows in this partition, subset of columns
    pixel_tbl = neowise_ds.to_table(filter=(KFIELD == pix), columns=COLUMN_SUBSET)

    # insert your code here to continue processing

    # we'll just print the table size to get a sense of how much data has been loaded
    print_table_size(table=pixel_tbl, pixel_index=pix)

    # when done, you may want to delete pixel_tbl to free the memory
    del pixel_tbl
    # we'll stop after one partition
    break
```

`pixel_tbl` is a (pyarrow) [Table](https://arrow.apache.org/docs/python/generated/pyarrow.Table.html)
containing all NEOWISE sources with an ra/dec falling within HEALPix order 5 pixel `pix`.
Use `pixel_tbl.to_pandas()` to convert the table to a pandas dataframe.

+++

How big are the partitions? (see Appendix for details)

```{code-cell} ipython3
# median partition
median_part_tbl = neowise_ds.to_table(
    filter=(KFIELD == median_part), columns=COLUMN_SUBSET
)
print_table_size(table=median_part_tbl, pixel_index=median_part)
```

Often only a few columns are needed for processing, so most partitions will fit comfortably in memory.
(The recommended maximum for an in-memory table/dataframe is typically 1GB, but there is
no strict upper limit -- performance will depend on the compute resources available.)

However, beware that the largest partitions are quite large:

```{code-cell} ipython3
# biggest partition
# this is very large, so we'll restrict the number of columns to one
biggest_part_tbl = neowise_ds.to_table(
    filter=(KFIELD == biggest_part), columns=COLUMN_SUBSET[:1]
)
print_table_size(table=biggest_part_tbl, pixel_index=biggest_part)

# Additional filters can be included to reduce the number of rows if desired.
# Another option is to load individual files.
```

```{code-cell} ipython3
# cleanup
del median_part_tbl
del biggest_part_tbl
```

## 5. Example: Slice by file

+++

If you don't need data for all years at the same time, you may want to load individual files.
Generally, there is 1 file per partition per year, but a few partitions are as large as 6+ files per year.
Most of the files are about 0.3GB (compressed on disk) but about 1% are > 1GB.
Thus it should be reasonable to load at least a subset of columns for every row in a file.

+++

The actual "slicing" here is done by using `neowise_ds` to access a
dataset [Fragment](https://arrow.apache.org/docs/python/generated/pyarrow.dataset.Fragment.html)
(`frag` in the code below) which represents a single file.

```{code-cell} ipython3
# slice by file and iterate
for frag in neowise_ds.get_fragments():
    # load the slice to get every row in the file, subset of columns
    file_tbl = frag.to_table(columns=COLUMN_SUBSET)

    # insert your code here to continue processing the file as desired

    # if you need to see which file this is, parse the path
    print(f"file path: {frag.path}")
    # let's see how much data this loaded
    print_table_size(table=file_tbl)

    # again, we'll stop after one
    del file_tbl
    break
```

This can be combined with the previous example to iterate over the files in a single partition
(left as an exercise for the reader).

+++

How big are the files?

```{code-cell} ipython3
# median file
median_file_frag = [
    frag for frag in neowise_ds.get_fragments() if frag.path == median_file
][0]
median_file_tbl = median_file_frag.to_table(columns=COLUMN_SUBSET)
print_table_size(table=median_file_tbl)
```

```{code-cell} ipython3
# biggest file
biggest_file_frag = [
    frag for frag in neowise_ds.get_fragments() if frag.path == biggest_file
][0]
biggest_file_tbl = biggest_file_frag.to_table(columns=COLUMN_SUBSET)
print_table_size(table=biggest_file_tbl)
```

```{code-cell} ipython3
# cleanup
del median_file_tbl
del biggest_file_tbl
```

## 6. Example: Slice by year

+++

If you want to handle the years independently you can work with the per-year datasets.
We actually created these "slices" in the Setup section with `year_datasets`, and that
same list is now accessible in `neowise_ds.children` used below.
Any of the techniques shown in this notebook or those listed under "See also" can also
be applied to the per-year datasets.

```{code-cell} ipython3
# slice by year and iterate. zip with YEARS so that we know which slice this is.
for year, year_ds in zip(YEARS, neowise_ds.children):
    # insert your code here to process year_ds as desired.
    # filter and load, iterate over partitions or files, etc.

    # we'll just look at some basic metadata.
    num_rows = sum(frag.metadata.num_rows for frag in year_ds.get_fragments())
    num_files = len(year_ds.files)
    print(f"NEOWISE {year} dataset: {num_rows:,} rows in {num_files:,} files")
```

## Appendix

+++

### A.1 Considerations when extending to specific use cases

Because the catalog is so large, you will need to carefully consider your specific problem and
determine how to slice and filter the data most efficiently.
There is no one right answer; it will depend on the use case.

#### A.1.1 Filtering

Filter out as much data as possible as early as possible. Ideas to consider are:

1. With the Parquet file format, you can apply filters during the read to avoid loading
   rows that you don't need.
      - Pandas (not demonstrated here) supports basic filters.
      - PyArrow (demonstrated here) also supports complex filters which allow you to compare
        values between columns and/or construct new columns on the fly (e.g., subtracting
        magnitude columns to construct a new color column, as done in the AllWISE notebook).
2. Queries (i.e., loading data by applying filters) will be *much* more efficient when they
   include a filter on the partitioning column ('healpix_k5'; demonstrated above).
      - Notice both that this is essentially equivalent to slicing by partition and that
        you can filter for more than one partition at a time.
      - This is highly recommended even if your use case doesn't explicitly care about it.
        Exceptions include situations where you're working with individual files and when
        it's impractical or counterproductive for the science.
3. You should also include filters specific to your use case if possible.
4. Exceptions: Sometimes it's not easy to write a dataset filter for the query.
  A cone search is a common example.
  In principal it could be written as a PyArrow dataset filter, but in practice the correct
  formula is much too complicated.
  In this case, it's easier to write dataset filters for broad RA and Dec limits and then
  do the actual cone search using `astropy`.
  This approach is still quite performant (see the NEOWISE light curves notebook).

#### A.1.2 Slicing

Slice the dataset in some way(s), then iterate and/or parallelize over the slices. Ideas to consider are:

1. Choose your slices so that you can:
      - Run your processing code on one slice independently. For example, if your code must
        see all the data for some target object (RA and Dec) at the same time, you may
        slice the dataset by partition, but don't slice it by year.
      - Load all data in the slice into memory at once (after applying your filters during the read).
        This notebook shows how to determine how big a slice of data is in order to guide this decision.
2. By default, slice by partition. If this is too much data, you may also want to slice
   by year and/or file.
3. If you have enough memory to load more than one slice simultaneously, parallelize over
   the slices to speed up your code.

### A.2 Inspect dataset stats

+++

When deciding how to slice and filter the dataset, it can be useful to understand
dataset statistics like partition and file sizes.

```{code-cell} ipython3
def pixel_index_from_path(path, k_column=KCOLUMN):
    """Parse the path and return the partition pixel index.

    Parameters
    ----------
    path : str
        The path to parse.
    k_column : str (optional)
        Name of the partitioning column.

    Returns
    -------
    int
        The partition pixel index parsed from the path.
    """
    pattern = rf"({k_column}=)([0-9]+)"  # matches strings like "healpix_k5=1124"
    return int(re.search(pattern, path).group(2))  # pixel index, e.g., 1124
```

```{code-cell} ipython3
# load some file statistics to a dataframe
file_stats = pd.DataFrame(
    columns=["path", KCOLUMN, "numrows"],
    data=[
        (frag.path, pixel_index_from_path(frag.path), frag.metadata.num_rows)
        for frag in neowise_ds.get_fragments()
    ],
)
```

```{code-cell} ipython3
file_stats.sample(5)
```

```{code-cell} ipython3
file_stats.describe()
```

#### A.2.1 Dataset statistics per file

```{code-cell} ipython3
# visualize distribution of file sizes (number of rows)
ax = file_stats.numrows.hist(log=True)
ax.set_xlabel("Number of rows")
ax.set_ylabel("Number of files")
```

```{code-cell} ipython3
# largest file
file_stats.loc[file_stats.numrows == file_stats.numrows.max()].head(1)
```

```{code-cell} ipython3
# median file
file_stats.sort_values("numrows").iloc[len(file_stats.index) // 2]
```

#### A.2.2 Dataset statistics per partition

```{code-cell} ipython3
# get stats per partition
k_groups = file_stats[[KCOLUMN, "numrows"]].groupby(KCOLUMN)
per_part = k_groups.sum()
per_part["numfiles"] = k_groups.count()
```

```{code-cell} ipython3
per_part.sample(5)
```

```{code-cell} ipython3
per_part.describe()
```

```{code-cell} ipython3
# visualize number of rows per partition
per_part.numrows.plot(
    logy=True, xlabel=f"{KCOLUMN} pixel index", ylabel="Number of rows per partition"
)
```

```{code-cell} ipython3
# largest partition
per_part.loc[per_part.numrows == per_part.numrows.max()]
```

```{code-cell} ipython3
# median partition
per_part.sort_values("numrows").iloc[len(per_part.index) // 2]
```

***

## About this notebook

**Author:** Troy Raen (IRSA Developer) and the IPAC Science Platform team

**Updated:** 2025-03-07

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
