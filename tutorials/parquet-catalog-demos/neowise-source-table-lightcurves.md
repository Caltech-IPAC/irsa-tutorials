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

# Make Light Curves from NEOWISE Single-exposure Source Table

+++

Learning Goals:

- Search the NEOWISE Single-exposure Source Table (Parquet version) for the light curves of a
  set of targets with RA/Dec coordinates.
  - Write a pyarrow dataset filter and use it to load the NEOWISE detections near each target (rough cut).
  - Match targets to detections using an astropy cone search (precise cut).
  - Parallelize this.
- Plot the light curves.

+++

## 1. Introduction

This notebook loads light curves from the
[NEOWISE](https://irsa.ipac.caltech.edu/Missions/wise.html) Single-exposure Source Table
for a sample of about 2000 cataclysmic variables from [Downes et al. (2001)](https://doi.org/10.1086/320802).
The NEOWISE Single-exposure Source Table is a very large catalog -- 10 years and 40 terabytes in total
with 145 columns and 200 billion rows.
When searching this catalog, it is important to consider the requirements of your use case and
the format of this dataset.
This notebook applies the techniques developed in
[Strategies to Efficiently Work with NEOWISE Single-exposure Source Table in Parquet](https://irsa.ipac.caltech.edu/docs/notebooks/neowise-source-table-strategies.html).
This is a fully-worked example that demonstrates the important steps, but note that this is a
relatively small use case for the Parquet version of the dataset.

The specific strategy we employ is:

- Choose a cone search radius that determines which NEOWISE source detections to associate
  with each target.
- Load the sample of CV targets.
- Calculate the indexes of all HEALPix order k=5 pixels within the radius of each target.
  These are the dataset partitions that need to be searched.
- Parallelize over the partitions using `multiprocessing.Pool`.
  For each pixel:
  - Construct a dataset filter for NEOWISE sources in the vicinity of the targets in the partition.
  - Load data, applying the filter. In our case, the number of rows loaded will be fairly small.
  - Do a cone search to match sources with targets in the partition.
  - Return the results.
- Concatenate the cone search results, groupby target ID, and sort by time to construct the light curves.

The efficiency of this method will increase with the number of rows needed from each partition.
For example, a cone search radius of 1 arcsec will require about 10 CPUs, 65G RAM, and
50 minutes to load the data from all 10 NEOWISE years.
Increasing the radius to 10 arcsec will return about 2.5x more rows using roughly the same resources.
Increasing the target sample size can result in similar efficiency gains.
To try out this notebook with fewer resources, use a subset of NEOWISE years.
Using one year is expected to require about 5 CPUs, 20G RAM, and 10 minutes.
These estimates are based on testing in science platform environments.
Your numbers will vary based on many factors including compute power, bandwidth, and physical distance from the data.

+++

## 2. Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install astropy astroquery hpgeom matplotlib pandas pyarrow pyvo
```

```{code-cell} ipython3
import multiprocessing  # parallelization

import astroquery.vizier  # fetch the sample of CV targets
import hpgeom  # HEALPix math
import numpy as np  # math
import pandas as pd  # manipulate tabular data
import pyarrow.compute  # construct dataset filters
import pyarrow.dataset  # load and query the NEOWISE dataset
import pyarrow.fs  # interact with the S3 bucket storing the NEOWISE catalog
import pyvo  # TAP service for the Vizier query
from astropy import units as u  # manipulate astropy quantities
from astropy.coordinates import SkyCoord  # manipulate sky coordinates
from matplotlib import pyplot as plt  # plot light curves

# copy-on-write will become the default in pandas 3.0 and is generally more performant
pd.options.mode.copy_on_write = True
```

## 3. Setup

+++

### 3.1 Define variables

First, choose which NEOWISE years to include.
Real use cases are likely to require all ten years but it can be helpful to start with
fewer while exploring to make things run faster.

```{code-cell} ipython3
YEARS = list(range(1, 11))  # all years => about 11 CPU, 65G RAM, and 50 minutes runtime

# To try out a smaller version of the notebook,
# uncomment the next line and choose your own subset of years.
# YEARS = [10]  # one year => about 5 CPU, 20G RAM, and 10 minutes runtime
```

```{code-cell} ipython3
# sets of columns that we'll need
FLUX_COLUMNS = ["w1flux", "w2flux"]
LIGHTCURVE_COLUMNS = ["mjd"] + FLUX_COLUMNS
COLUMN_SUBSET = ["cntr", "ra", "dec"] + LIGHTCURVE_COLUMNS

# cone-search radius defining which NEOWISE sources are associated with each target object
MATCH_RADIUS = 1 * u.arcsec
```

### 3.2 Load NEOWISE metadata

+++

The metadata contains column names, schema, and row-group statistics for every file in the dataset.
We'll load it as a pyarrow dataset.

```{code-cell} ipython3
# This catalog is so big that even the metadata is big.
# Expect this cell to take about 30 seconds per year.

# This information can be found at https://irsa.ipac.caltech.edu/cloud_access/.
bucket = "nasa-irsa-wise"
base_prefix = "wise/neowiser/catalogs/p1bs_psd/healpix_k5"
metadata_path = (
    lambda yr: f"{bucket}/{base_prefix}/year{yr}/neowiser-healpix_k5-year{yr}.parquet/_metadata"
)
fs = pyarrow.fs.S3FileSystem(region="us-west-2", anonymous=True)

# list of datasets, one per year
year_datasets = [
    pyarrow.dataset.parquet_dataset(metadata_path(yr), filesystem=fs, partitioning="hive")
    for yr in YEARS
]

# unified dataset, all years
neowise_ds = pyarrow.dataset.dataset(year_datasets)
```

## 4. Define functions to filter and load data

+++

These functions will be used in the next section.
Defining them here in the notebook is useful for demonstration and should work seamlessly on Linux, which includes most science platforms.
Mac and Windows users should see the note at the end of the notebook.
However, note that this use case is likely too large for a laptop and may perform poorly and/or crash if attempted.

```{code-cell} ipython3
# If you have your own list of target objects, replace this function to load your sample.
def load_targets_Downes2001(radius=1 * u.arcsec):
    """Load a sample of targets and return a pandas DataFrame.

    References:
    - Downes et al., 2001 ([2001PASP..113..764D](https://ui.adsabs.harvard.edu/abs/2001PASP..113..764D/abstract)).
    - https://cdsarc.cds.unistra.fr/ftp/V/123A/ReadMe

    Parameters
    ----------
    radius : astropy.Quantity (optional)
        Radius for the cone search around each target. This is used to determine which partition(s)
        need to be searched for a given target. Use the same radius here as in the rest of the notebook.

    Returns
    -------
    pandas.DataFrame
        The loaded targets with the following columns:
            - uid: Unique identifier of the target.
            - GCVS: Name in the General Catalogue of Variable Stars if it exists, else the constellation name.
            - RAJ2000: Right Ascension of the target in J2000 coordinates.
            - DEJ2000: Declination of the target in J2000 coordinates.
            - healpix_k5: HEALPix pixel index at order k=5.
    """
    astroquery.vizier.Vizier.ROW_LIMIT = -1
    # https://cdsarc.cds.unistra.fr/vizier/notebook.gml?source=V/123A
    # https://cdsarc.cds.unistra.fr/ftp/V/123A/ReadMe
    CATALOGUE = "V/123A"
    voresource = pyvo.registry.search(ivoid=f"ivo://CDS.VizieR/{CATALOGUE}")[0]
    tap_service = voresource.get_service("tap")

    # Query Vizier and load targets to a dataframe.
    cv_columns = ["uid", "GCVS", "RAJ2000", "DEJ2000"]
    cvs_records = tap_service.run_sync(
        f'SELECT {",".join(cv_columns)} from "{CATALOGUE}/cv"'
    )
    cvs_df = cvs_records.to_table().to_pandas()

    # Add a new column containing a list of all order k HEALPix pixels that overlap with
    # the CV's position plus search radius.
    cvs_df["healpix_k5"] = [
        hpgeom.query_circle(
            a=cv.RAJ2000,
            b=cv.DEJ2000,
            radius=radius.to_value("deg"),
            nside=hpgeom.order_to_nside(order=5),
            nest=True,
            inclusive=True,
        )
        for cv in cvs_df.itertuples()
    ]
    # Explode the lists of pixels so the dataframe has one row per target per pixel.
    # Targets near a pixel boundary will now have more than one row.
    # Later, we'll search each pixel separately for NEOWISE detections and then
    # concatenate the matches for each target to produce complete light curves.
    cvs_df = cvs_df.explode("healpix_k5", ignore_index=True)

    return cvs_df
```

```{code-cell} ipython3
# This is the main function.
def load_lightcurves_one_partition(targets_group):
    """Load lightcurves from a single partition.

    Parameters
    ----------
    targets_group : tuple
        Tuple of pixel index and sub-DataFrame (result of DataFrame groupby operation).

    Returns
    -------
    pd.DataFrame
        The lightcurves for targets found in this partition.
    """
    # These global variables will be set when the worker is initialized.
    global _neowise_ds
    global _columns
    global _radius

    # Get row filters that will limit the amount of data loaded from this partition.
    # It is important for these filters to be efficient for the specific use case.
    filters = _construct_dataset_filters(targets_group=targets_group, radius=_radius)

    # Load this slice of the dataset to a pyarrow Table.
    pixel_tbl = _neowise_ds.to_table(columns=_columns, filter=filters)

    # Associate NEOWISE detections with targets to get the light curves.
    lightcurves_df = _cone_search(
        targets_group=targets_group, pixel_tbl=pixel_tbl, radius=_radius
    )

    return lightcurves_df
```

```{code-cell} ipython3
# The filters returned by this function need to be efficient for the specific use case.
def _construct_dataset_filters(*, targets_group, radius, scale_factor=100):
    """Construct dataset filters for a box search around all targets in the partition.

    Parameters
    ----------
    targets_group : tuple
        Tuple of pixel index and sub-DataFrame (result of DataFrame groupby operation).
    radius : astropy.Quantity
        The radius used for constructing the RA and Dec filters.
    scale_factor : int (optional)
        Factor by which the radius will be multiplied to ensure that the box encloses
        all relevant detections.

    Returns
    -------
    filters : pyarrow.compute.Expression
        The constructed filters based on the given inputs.
    """
    pixel, locations_df = targets_group

    # Start with a filter for the partition. This is the most important one because
    # it tells the Parquet reader to just skip all the other partitions.
    filters = pyarrow.compute.field("healpix_k5") == pixel

    # Add box search filters. For our CV sample, one box encompassing all targets in
    # the partition should be sufficient. Make a different choice if you use a different
    # sample and find that this loads more data than you want to handle at once.
    buffer_dist = scale_factor * radius.to_value("deg")
    for coord, target_coord in zip(["ra", "dec"], ["RAJ2000", "DEJ2000"]):
        coord_fld = pyarrow.compute.field(coord)

        # Add a filter for coordinate lower limit.
        coord_min = locations_df[target_coord].min()
        filters = filters & (coord_fld > coord_min - buffer_dist)

        # Add a filter for coordinate upper limit.
        coord_max = locations_df[target_coord].max()
        filters = filters & (coord_fld < coord_max + buffer_dist)

    # Add your own additional requirements here, like magnitude limits or quality cuts.
    # See the AllWISE notebook for more filter examples and links to pyarrow documentation.
    # We'll add a filter for sources not affected by contamination or confusion.
    filters = filters & pyarrow.compute.equal(pyarrow.compute.field("cc_flags"), "0000")

    return filters
```

```{code-cell} ipython3
def _cone_search(*, targets_group, pixel_tbl, radius):
    """Perform a cone search to select NEOWISE detections belonging to each target object.

    Parameters
    ----------
    targets_group : tuple
        Tuple of pixel index and sub-DataFrame (result of DataFrame groupby operation).
    pixel_tbl : pyarrow.Table
        Table of NEOWISE data for a single pixel
    radius : astropy.Quantity
        Cone search radius.

    Returns
    -------
    match_df : pd.DataFrame
        A dataframe with all matched sources.
    """
    _, targets_df = targets_group

    # Cone search using astropy to select NEOWISE detections belonging to each object.
    pixel_skycoords = SkyCoord(ra=pixel_tbl["ra"] * u.deg, dec=pixel_tbl["dec"] * u.deg)
    targets_skycoords = SkyCoord(targets_df["RAJ2000"], targets_df["DEJ2000"], unit=u.deg)
    targets_ilocs, pixel_ilocs, _, _ = pixel_skycoords.search_around_sky(
        targets_skycoords, radius
    )

    # Create a dataframe with all matched source detections.
    match_df = pixel_tbl.take(pixel_ilocs).to_pandas()

    # Add the target IDs by joining with targets_df.
    match_df["targets_ilocs"] = targets_ilocs
    match_df = match_df.set_index("targets_ilocs").join(targets_df.reset_index().uid)

    return match_df
```

```{code-cell} ipython3
# This function will be called once for each worker in the pool.
def init_worker(neowise_ds, columns, radius):
    """Set global variables '_neowise_ds', '_columns', and '_radius'.

    These variables will be the same for every call to 'load_lightcurves_one_partition'
    and will be set once for each worker. It is important to pass 'neowise_ds' this
    way because of its size and the way it will be used. (For the other two, it makes
    little difference whether we use this method or pass them directly as function
    arguments to 'load_lightcurves_one_partition'.)

    Parameters
    ----------
    neowise_ds : pyarrow.dataset.Dataset
        NEOWISE metadata loaded as a PyArrow dataset.
    columns : list
        Columns to include in the output DataFrame of light curves.
    radius : astropy.Quantity
        Cone search radius.
    """
    global _neowise_ds
    _neowise_ds = neowise_ds
    global _columns
    _columns = columns
    global _radius
    _radius = radius

```

## 5. Load light curves

+++

Load the target objects' coordinates and other info.

```{code-cell} ipython3
targets_df = load_targets_Downes2001(radius=MATCH_RADIUS)
targets_df.head()
```

Search the NEOWISE Source Table for all targets (positional matches) and load the light curves.
Partitions are searched in parallel.
For targets located near partition boundaries, relevant partitions will be searched
independently for the given target and the results will be concatenated.
If searching all NEOWISE years, this may take 45 minutes or more.

```{code-cell} ipython3
# Group targets by partition. 'load_lightcurves_one_partition' will be called once per group.
targets_groups = targets_df.groupby("healpix_k5")
# Arguments for 'init_worker'.
init_args = (neowise_ds, COLUMN_SUBSET, MATCH_RADIUS)

# Start a multiprocessing pool and load the target light curves in parallel.
# About 1900 unique pixels in targets_df, 8 workers, 48 chunksize => ~5 chunks per worker.
nworkers = 8
chunksize = 48
with multiprocessing.Pool(nworkers, initializer=init_worker, initargs=init_args) as pool:
    lightcurves = []
    for lightcurves_df in pool.imap_unordered(
        load_lightcurves_one_partition, targets_groups, chunksize=chunksize
    ):
        lightcurves.append(lightcurves_df)
neowise_lightcurves_df = pd.concat(lightcurves).sort_values("mjd").reset_index(drop=True)
```

```{code-cell} ipython3
neowise_lightcurves_df.head()
```

## 6. Plot NEOWISE light curves

The light curves will have large gaps due to the observing cadence, so we'll plot each
"epoch" separately to see them better.

```{code-cell} ipython3
# get the light curves of the target with the most data
target_uid = neowise_lightcurves_df.groupby("uid").mjd.count().sort_values().index[-1]
target_df = neowise_lightcurves_df.loc[neowise_lightcurves_df.uid == target_uid]

# list of indexes that separate epochs (arbitrarily at delta mjd > 30)
epoch_idxs = target_df.loc[target_df.mjd.diff() > 30].index.to_list()
epoch_idxs = epoch_idxs + [target_df.index[-1]]  # add the final index

# make the figure
ncols = 4
nrows = int(np.ceil(len(epoch_idxs) / ncols))
fig, axs = plt.subplots(nrows, ncols, sharey=True, figsize=(3 * ncols, 2.5 * nrows))
axs = axs.flatten()
idx0 = target_df.index[0]
for i, (idx1, ax) in enumerate(zip(epoch_idxs, axs)):
    epoch_df = target_df.loc[idx0 : idx1 - 1, LIGHTCURVE_COLUMNS].set_index("mjd")
    for col in FLUX_COLUMNS:
        ax.plot(epoch_df[col], ".", markersize=3, label=col)
    ax.set_title(f"epoch {i}")
    ax.xaxis.set_ticks(  # space by 10
        range(int((ax.get_xlim()[0] + 10) / 10) * 10, int(ax.get_xlim()[1]), 10)
    )
    idx0 = idx1
axs[0].legend()
fig.supxlabel("MJD")
fig.supylabel("RAW FLUX")
fig.suptitle(f"NEOWISE light curves for target CV {target_uid}")
fig.tight_layout()
plt.show(block=False)
```

-----

[*] Note to Mac and Windows users:

You will need to copy the functions and imports from this notebook into a separate '.py' file and then import them in order to use the multiprocessing pool for parallelization.
In addition, you may need to load `neowise_ds` separately for each child process (i.e., worker) by adding that code to the `init_worker` function instead of passing it in as an argument.
This has to do with differences in what does / does not get copied into the child processes on different platforms.

About this notebook:

- Author: Troy Raen (Applications Developer, IRSA) and the IPAC Science Platform team
- Contact: [https://irsa.ipac.caltech.edu/docs/help_desk.html](https://irsa.ipac.caltech.edu/docs/help_desk.html)
- Updated: 2024-08-08
