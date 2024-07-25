---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---



# CosmoDC2 Mock v1 catalogs with IRSA TAP

This tutorial demonstrates how to access the CosmoDC2 Mock V1 catalogs. More information about these catalogs can be found here: https://irsa.ipac.caltech.edu/Missions/cosmodc2.html

These catalogs can be accessed through IRSA's Virtual Ovservatory Table Access Protocol (TAP) service. See https://www.ivoa.net/documents/TAP/ for details on the protocol. This service can be accessed through Python using the PyVO library.

```{code-cell} ipython3
import pyvo as vo
```

```{code-cell} ipython3
service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
```

## List the available DC2 tables

```{code-cell} ipython3
tables = service.tables
for tablename in tables.keys():
    if not "tap_schema" in tablename:
        if "dc2" in tablename:
            tables[tablename].describe()
```

## Choose the DC2 catalog you want to work with.

IRSA currently offers 3 versions of the DC2 catalog.

* ``cosmodc2mockv1_new`` has been optimized to make searches with constraints on stellar mass and redshift fast.

* ``cosmodc2mockv1`` has been optimized to make searches with spatial constraints fast.

* ``cosmodc2mockv1_heavy`` is the same as ``cosmodc2mockv1_new``, except that it does not contain galaxies with stellar masses <= 10^7 solar masses.

If you are new to the DC2 catalog, we recommend that you start with ``cosmodc2mockv1_heavy``

```{code-cell} ipython3
# Choose the abridged table to start with.
# Queries should be faster on smaller tables.

tablename = 'cosmodc2mockv1_heavy'
```

## How many rows are in the chosen table?

With TAP, you can query catalogs with constraints specified in IVOA Astronomical Data Query Language (ADQL; https://www.ivoa.net/documents/latest/ADQL.html), which is based on SQL.

```{code-cell} ipython3
# For example, this snippet of ADQL counts the number of elements in
# the redshift column of the table we chose.
adql = f"SELECT count(redshift) FROM {tablename}"
adql
```

```{code-cell} ipython3
# In order to use TAP with this ADQL string using pyvo, you can do the following:
result = service.search(adql)
result
```

The above query shows that there are 597,488,849 redshifts in this table.

+++

## What is the default maximum number of rows returned by the service?

This service will return a maximum of 2 billion rows by default.

```{code-cell} ipython3
service.maxrec
```

This default maximum can be changed, and there is no hard upper limit to what it can be changed to.

```{code-cell} ipython3
print(service.hardlimit)
```

## List the columns in the chosen table

This table contains 301 columns.

```{code-cell} ipython3
columns = tables[tablename].columns
print(len(columns))
```

Let's learn a bit more about them.

```{code-cell} ipython3
for col in columns:
    print(f'{f"{col.name}":30s}  {col.description}')
```

## Create a histogram of redshifts

Let's figure out what redshift range these galaxies cover. Since we found out above that it's a large catalog, we can start with a spatial search over a small area of 0.1 deg. The ADQL that is needed for the spatial constraint is:

```{code-cell} ipython3
adql = f"SELECT redshift FROM {tablename} WHERE CONTAINS(POINT('ICRS', RAMean, DecMean), CIRCLE('ICRS',54.218205903,-37.497959343,.1))=1"
adql
```

Now we can use the previously-defined service to execute the query with the spatial contraint.

```{code-cell} ipython3
cone_results = service.search(adql)
```

```{code-cell} ipython3
# Plot a histogram
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

num_bins = 20
# the histogram of the data
n, bins, patches = plt.hist(cone_results['redshift'], num_bins,
                            facecolor='blue', alpha = 0.5)
plt.xlabel('Redshift')
plt.ylabel('Number')
plt.title('Redshift Histogram CosmoDC2 Mock Catalog V1 abridged')
```

We can easily see form this plot that the simulated galaxies go out to z = 3.

+++

## Visualize galaxy colors at z ~ 0.5

Now let's visualize the galaxy main sequence at z = 2.0. First, we'll do a narrow redshift cut with no spatial constraint.

Let's do it as an asynchronous search since this might take awhile.

```{code-cell} ipython3
service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
adql = f"SELECT Mag_true_r_sdss_z0, Mag_true_g_sdss_z0, redshift FROM {tablename} WHERE redshift > 0.5 and redshift < 0.54"
results = service.run_async(adql)
```

```{code-cell} ipython3
len(results['mag_true_r_sdss_z0'])
```

```{code-cell} ipython3
# Since this results in almost 4 million galaxies,
# we will construct a 2D histogram rather than a scatter plot.
plt.hist2d(results['mag_true_r_sdss_z0'], results['mag_true_g_sdss_z0']-results['mag_true_r_sdss_z0'],
           bins=200, cmap='plasma', cmax=500)

# Plot a colorbar with label.
cb = plt.colorbar()
cb.set_label('Number')

# Add title and labels to plot.
plt.xlabel('SDSS Mag r')
plt.ylabel('SDSS rest-frame g-r color')

# Show the plot.
plt.show()
```

## About this notebook

 * Author: Vandana Desai (IRSA Science Lead)
 * Updated: 2024-07-24
 * Contact: https://irsa.ipac.caltech.edu/docs/help_desk.html
