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
import top

RUN_ID = "euclid-hats-query-methods"
top.tag(run_id=RUN_ID, time="start-run")
```

```{code-cell}
s3_bucket = "nasa-irsa-euclid-q1"
euclid_prefix = "contributed/q1/merged_objects/hats"

euclid_hats_collection_uri = f"s3://{s3_bucket}/{euclid_prefix}"  # for lsdb
euclid_parquet_metadata_path = f"{s3_bucket}/{euclid_prefix}/euclid_q1_merged_objects-hats/dataset/_metadata"  # for pyarrow

max_magnitude = 24.5
min_flux = 10 ** ((max_magnitude - 23.9) / -2.5)
```

```{code-cell}
# Columns we actually want to load.
OBJECT_ID = "object_id"
PHZ_Z = "phz_phz_median"
columns = [OBJECT_ID, PHZ_Z]
```

## Load with pyarrow

```{code-cell}
%%time
top.tag(run_id=RUN_ID, time="pyarrow")

# Construct filter for quality PHZ redshifts.
phz_filter = (
    (pc.field("mer_vis_det") == 1)  # No NIR-only objects.
    & (pc.field("mer_flux_detection_total") > min_flux)  # I < 24.5
    & (pc.divide(pc.field("mer_flux_detection_total"), pc.field("mer_fluxerr_detection_total")) > 5)  # I band S/N > 5
    & ~pc.field("phz_phz_classification").isin([1, 3, 5, 7])  # Exclude objects classified as star.
    & (pc.field("mer_spurious_flag") == 0)  # MER quality
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
top.tag(run_id=RUN_ID, time="pyarrow + dask")


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
top.tag(run_id=RUN_ID, time="lsdb")

# Construct the query equivalent of phz_filter.
query = (
    "mer_vis_det == 1"
    f" & mer_flux_detection_total > {min_flux}"
    " & mer_flux_detection_total / mer_fluxerr_detection_total > 5"
    " & phz_phz_classification not in [1,3,5,7]"
    " & mer_spurious_flag == 0"
)

# We don't want to load these columns, but we have to in order to use them in the filter.
extra_columns = ["mer_vis_det", "mer_flux_detection_total", "phz_phz_classification", "mer_spurious_flag", "mer_fluxerr_detection_total"]

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

top.tag(run_id=RUN_ID, time="end-run")


tl = top.load_top_output(run_id=RUN_ID, named_pids_only=False)
fig = tl.plot_overview()
fig.savefig(tl.base_dir / "top.png")
```
