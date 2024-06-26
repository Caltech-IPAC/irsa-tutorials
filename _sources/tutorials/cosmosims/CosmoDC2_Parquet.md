---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Working with CosmoDC2 Parquet files

This notebook shows how to work with the bulk-downloadable CosmoDC2 files in Parquet format. Goals include:

* how to read in a single file
* how to plot galaxy properties
* how to read in selected columns
* how to read from all of the files

+++

### Imports

The pyarrow package is the underlying engine for reading these files

```{code-cell} ipython3
import os

from astropy.table import Table
import dask.dataframe as dd
import numpy as np
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import pandas as pd
import pyarrow.parquet as pq
```

```{code-cell} ipython3
%matplotlib inline
```

### Set data location

+++

Modify the `inpath` to your location for the Parquet files.

```{code-cell} ipython3
inpath = '/stage/irsa-data-download10/parquet-work/CosmoDC2'
```

Select a single file for the first demonstrations. The files are numbered by HEALPix number with NSIDE=32. The example here is near the center of the 440 sq deg covered by the simulations.

```{code-cell} ipython3
infile = 'cosmoDC2_healpix09560.parquet'
```

```{code-cell} ipython3
joined = os.path.join(inpath, infile)
```

### Get all column names and types

+++

Read a table of zero length to get the column names.

```{code-cell} ipython3
schema = Table.read(joined, schema_only=True)
```

The column names are in case-sensitive alphabetical order here.

```{code-cell} ipython3
schema
```

Alternatively, use the pyarrow package directly.

```{code-cell} ipython3
schema = pq.read_schema(os.path.join(inpath, infile))
```

```{code-cell} ipython3
len(schema)
```

Capitalized magnitudes are absolute magnitudes. Descriptions of the columns are not included in the Parquet files, but are in a separate file (to be specified).

```{code-cell} ipython3
print(schema)
```

### Read only the redshifts as a single Astropy table

Reading a single column is faster with Parquet files.

```{code-cell} ipython3
tab = Table.read(joined, include_names=['redshift'])
```

```{code-cell} ipython3
tab
```

### Read in all data from a single file

```{code-cell} ipython3
d=pd.read_parquet(joined)
```

```{code-cell} ipython3
len(d)
```

### Plot number counts vs halo mass

```{code-cell} ipython3
plt.figure()
h,xbins = np.histogram(np.log10(d['halo_mass']),bins=40)
xbins_avg = (xbins[1:]+xbins[:-1])/2.0
plt.semilogy(xbins_avg, h)
plt.ylabel(r'Galaxy Count')
plt.xlabel(r'log10( M$_{\rm{halo}}$ / M$_\odot)$')
plt.show()
```

### Plot g-r colors vs redshift

```{code-cell} ipython3
plt.figure()
gal_clr = d['mag_g_sdss']-d['mag_r_sdss']
plt.hist2d(d['redshift'], gal_clr, bins=100, cmap='PuBu', norm=clr.LogNorm())
plt.colorbar(label='population density')
plt.ylabel('Observed g-r')
plt.xlabel('redshift')
plt.title('Galaxy Colors in Clusters')
plt.tight_layout()
```

### Plot r-i colors vs redshift

```{code-cell} ipython3
plt.figure()
gal_clr = d['mag_r_sdss']-d['mag_i_sdss']
plt.hist2d(d['redshift'], gal_clr, bins=100, cmap='PuBu',norm=clr.LogNorm())
plt.colorbar(label='population density')
plt.ylabel('r-i')
plt.xlabel('redshift')
plt.title('Galaxy Colors in Clusters')
plt.tight_layout()
plt.show()
```

### Read one Parquet file using Pandas, and plot a histogram of redshifts

```{code-cell} ipython3
df = pd.read_parquet(joined, columns=['redshift'])
```

```{code-cell} ipython3
num_bins = 20
n, bins, patches = plt.hist(df.redshift, num_bins,
                            facecolor='blue', alpha=0.5)
plt.xlabel('Redshift')
plt.ylabel('Number')
plt.title('Redshift Histogram CosmoDC2 Mock Catalog V1 Abridged')
```

### Count the number of records in all the files

+++

This section assumes you have downloaded all the files, or at least multiple files.

```{code-cell} ipython3
ddf = dd.read_parquet(os.path.join(inpath, '*parquet'),
                      engine='pyarrow',
                      columns=['stellar_mass'])
```

```{code-cell} ipython3
%%time
len(ddf)
```

### Read multiple files using Dask, and plot the galaxy main sequence in a redshift interval

+++

This section assumes you have downloaded all the files, or at least multiple files.
Read in stellar mass and absolute magnitudes in a redshift interval from 0.5 to 0.54.

```{code-cell} ipython3
ddf = dd.read_parquet(os.path.join(inpath, '*parquet'),
                      engine='pyarrow',
                      columns=['stellar_mass',
                               'Mag_true_g_lsst_z0',
                               'Mag_true_r_lsst_z0',
                               'Mag_true_i_lsst_z0',
                               'Mag_true_z_lsst_z0',
                               'Mag_true_y_lsst_z0',
                               'redshift'],
                      filters=[('redshift', '>', 0.5), 
                               ('redshift', '<', 0.54)])
```

The next cell does the data read and can take several minutes.

```{code-cell} ipython3
%%time
df2d = ddf.compute()
```

```{code-cell} ipython3
df2d.head()
```

```{code-cell} ipython3
len(df2d)
```

Since this results in almost 4 million galaxies, we will construct 2D histograms 
rather than scatter plots.

```{code-cell} ipython3

plt.hist2d(df2d.Mag_true_r_lsst_z0, 
           (df2d.Mag_true_g_lsst_z0 - df2d.Mag_true_r_lsst_z0),
           bins=200, cmap='plasma', cmax = 500)

# Plot a colorbar with label.
cb = plt.colorbar()
cb.set_label('Number')

# Add title and labels to plot.
plt.title('True magnitudes and colors in 0.5 < z < 0.54')
plt.xlabel('LSST Mag r')
plt.ylabel('LSST rest-frame g-r color')

# Show the plot.
plt.show()
```

```{code-cell} ipython3

```
