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

# Euclid Q1 HATS query - parallelization methods

This notebook demonstrates a couple of ways that a query for a quality Euclid Q1 redshift sample could be parallelized.

+++

## Setup

```{code-cell}
# # Uncomment the next line to install dependencies if needed.
%pip install 'lsdb>=0.5.2' 'numpy>=2.0' pyarrow s3fs
```

```{code-cell}
import os
import sys
import lsdb
import dask
import dask.distributed
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset
from upath import UPath
```

```{code-cell}
s3_bucket = "irsa-fornax-testdata"
euclid_prefix = "EUCLID/q1/catalogues"

euclid_hats_collection_uri = UPath(f"s3://{s3_bucket}/{euclid_prefix}")  # for lsdb
euclid_parquet_metadata_path = f"{s3_bucket}/{euclid_prefix}/hats/dataset/_metadata"  # for pyarrow

max_magnitude = 24.5
min_flux = 10 ** ((max_magnitude - 23.9) / -2.5)
```

```{code-cell}
# Columns we actually want to load.
OBJECT_ID = "OBJECT_ID"
PHZ_Z = "PHZ_PHZ_MEDIAN"
columns = [OBJECT_ID, PHZ_Z]
```

## Load with pyarrow

```{code-cell}
%%time
# Construct filter for quality PHZ redshifts.
phz_filter = (
    (pc.field("MER_VIS_DET") == 1)  # No NIR-only objects.
    & (pc.field("MER_FLUX_DETECTION_TOTAL") > min_flux)  # I < 24.5
    & (pc.divide(pc.field("MER_FLUX_DETECTION_TOTAL"), pc.field("MER_FLUXERR_DETECTION_TOTAL")) > 5)  # I band S/N > 5
    & ~pc.field("PHZ_PHZ_CLASSIFICATION").isin([1, 3, 5, 7])  # Exclude objects classified as star.
    & (pc.field("MER_SPURIOUS_FLAG") == 0)  # MER quality
)

# Load.
dataset = pyarrow.dataset.parquet_dataset(euclid_parquet_metadata_path, partitioning="hive", filesystem=pyarrow.fs.S3FileSystem())
pa_df = dataset.to_table(columns=columns, filter=phz_filter).to_pandas()
# 1 - 1.5 min

print(f"Pyarrow loaded {sys.getsizeof(pa_df) / 1024**3:.3f}G")
```

```{code-cell}
pa_df = pa_df.sort_values(OBJECT_ID, ignore_index=True)
pa_df
```

## Load with pyarrow + dask

The above query probably doesn't need to be parallelized, but we will want to parallelize similar queries in the future.
One option is to have dask workers execute the pyarrow filter and load.
This query has no spatial component and needs to look at all the files, so below is a basic implementation that just distributes the files.
(I wonder if this could be made better by batching the files, to reduce overhead?)

```{code-cell}
%%time
@dask.delayed
def load_fragment(frag):
    table = frag.to_table(filter=phz_filter, columns=columns)
    return table.to_pandas()


client = dask.distributed.Client(n_workers=os.cpu_count(), threads_per_worker=2, memory_limit=None)
delayed_dfs = [load_fragment(frag) for frag in dataset.get_fragments()]
padask_df = pd.concat(dask.compute(*delayed_dfs), ignore_index=True)
client.close()
# 40 sec - 1 min (4, 8, and 16 workers)

print(f"Pyarrow + dask loaded {sys.getsizeof(padask_df) / 1024**3:.3f}G")
# Ignoring the warnings about the s3 connection for now.
```

```{code-cell}
padask_df = padask_df.sort_values(OBJECT_ID, ignore_index=True)
padask_df
```

```{code-cell}
# Check for equality.
pa_df.equals(padask_df)
```

## Load with lsdb

```{code-cell}
%%time
# Construct the query equivalent of phz_filter.
query = (
    "MER_VIS_DET == 1"
    f" & MER_FLUX_DETECTION_TOTAL > {min_flux}"
    " & MER_FLUX_DETECTION_TOTAL / MER_FLUXERR_DETECTION_TOTAL > 5"
    " & PHZ_PHZ_CLASSIFICATION not in [1,3,5,7]"
    " & MER_SPURIOUS_FLAG == 0"
)

# We don't want to load these columns, but we have to in order to use them in the filter.
extra_columns = ["MER_VIS_DET", "MER_FLUX_DETECTION_TOTAL", "PHZ_PHZ_CLASSIFICATION", "MER_SPURIOUS_FLAG", "MER_FLUXERR_DETECTION_TOTAL"]

# Load.
client = dask.distributed.Client(n_workers=os.cpu_count(), threads_per_worker=2, memory_limit=None)
lsdb_catalog = lsdb.read_hats(euclid_hats_collection_uri, columns=columns + extra_columns)
lsdb_df = lsdb_catalog.query(query).compute()
client.close()
# 4.5 - 5.5 min (16 workers)
# 8 - 9 min (8 workers)
# 21 min (4 workers)

print(f"Pyarrow loaded {sys.getsizeof(lsdb_df) / 1024**3:.3f}G")
```

```{code-cell}
lsdb_df = lsdb_df.sort_values(OBJECT_ID, ignore_index=True)
lsdb_df
```

```{code-cell}
# Check for equality.
lsdb_df[PHZ_Z].astype("float32").equals(pa_df[PHZ_Z])
```
