---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
execution:
  timeout: 2600
---

# Querying the CosmoDC2 Mock v1 Catalogs

This tutorial demonstrates how to access and query the **CosmoDC2 Mock v1** catalogs using IRSA’s Table Access Protocol (TAP) service. Background information on the catalogs is available on the [IRSA CosmoDC2 page](https://irsa.ipac.caltech.edu/Missions/cosmodc2.html).

The catalogs are served through IRSA’s Virtual Observatory–standard **TAP** interface (see the [IVOA TAP specification](https://www.ivoa.net/documents/TAP/)), which you can access programmatically in Python via the **PyVO** library. TAP queries are written in the **Astronomical Data Query Language (ADQL)** — a SQL-like language designed for astronomical catalogs (see the [ADQL specification](https://www.ivoa.net/documents/latest/ADQL.html)).

If you are new to PyVO’s query modes, the documentation provides a helpful comparison between **synchronous** and **asynchronous** execution:  [PyVO: Synchronous vs. Asynchronous Queries](https://pyvo.readthedocs.io/en/latest/dal/index.html#synchronous-vs-asynchronous-query)


## Tips for Working with CosmoDC2 via TAP

- **Use indexed columns for fast queries.**
  CosmoDC2 is indexed on the following fields:
  `ra`, `dec`, `redshift`, `mag*_lsst`, `halo_mass`, `stellar_mass`
  Queries involving these columns generally return much faster.

- **Ensure your positional queries fall within the survey footprint.**
  CosmoDC2 roughly covers:
  - **RA:** 0° → 60°
  - **Dec:** –45° → 0°
  (Coverage is not perfectly rectangular, so some edges may be sparse.)

- **Avoid overloading the TAP service.**
  Prefer **asynchronous** queries (`submit_job` / `run`) so they can be monitored or cancelled.
  Using `run_sync()` makes it easy to fire off queries that cannot be interrupted and may continue running on the server long after your client session stops.

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy matplotlib pyvo
```

```{code-cell} ipython3
import pyvo as vo
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import time
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

# testing

```{code-cell} ipython3
# Test the TAP service health
test_adql = f"SELECT TOP 1 redshift FROM {tablename}"
job = service.submit_job(test_adql)
job.run()
results = job.fetch_result()
print(len(results))
```

```{code-cell} ipython3
job.phase
```

```{code-cell} ipython3
import time

def wait_for_job(job, timeout=60, poll_interval=5.0):
    """
    Wait until a TAP job reaches COMPLETED or ERROR.
    If it is still EXECUTING after timeout seconds, cancel it.

    Parameters
    ----------
    job : pyvo.dal.AsyncTAPJob
        The TAP job object from service.submit_job().
    timeout : int or float
        Max time (seconds) to wait before auto-cancelling.
    poll_interval : float
        Seconds between job.phase checks.

    Returns
    -------
    bool
        True if COMPLETED; False if ERROR, ABORTED, or timeout-cancelled.
    """
    start = time.time()

    while True:
        phase = job.phase
        print("Phase:", phase)

        # Normal success
        if phase == "COMPLETED":
            return True

        # Hard failure
        if phase in ("ERROR", "ABORTED", "UNKNOWN"):
            return False

        # Timeout reached — cancel it
        if time.time() - start > timeout:
            print(f"Timeout reached ({timeout} s). Aborting TAP job...")
            job.abort()
            return False

        time.sleep(poll_interval)
```

```{code-cell} ipython3
import time

# ------- Test A: (baseline speed) -------

test_adql_A = f"SELECT TOP 1 redshift FROM {tablename}"
t0 = time.time()
jobA = service.submit_job(test_adql_A)
t_submit = time.time()
print("submit job took", t_submit - t0, "s")

jobA.run()
if wait_for_job(jobA, timeout=60):
    results = jobA.fetch_result()
    print("Rows:", len(results))
else:
    print("Job did NOT finish — cancelled or failed.")

tA = time.time() - t_submit

print("No-WHERE TOP1 query time:", tA, "seconds")
```

```{code-cell} ipython3
test_adql_A = f"SELECT TOP 1 redshift FROM {tablename}"
results = service.run_async(test_adql_A)
```

```{code-cell} ipython3
results
```

```{code-cell} ipython3
# ------- Test B: RA/Dec CONTAINS query (index test) -------
test_adql_B = f"""
SELECT TOP 1 redshift
FROM {tablename}
WHERE CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', 54.2, -37.4, 0.05)
) = 1
"""

jobB = service.submit_job(test_adql_B)

t0 = time.time()
jobB.run()

if wait_for_job(jobB, timeout=150):
    results = jobB.fetch_result()
    print("Rows:", len(results))
else:
    print("Job did NOT finish — cancelled or failed.")

tB = time.time() - t0

print("Cone-search TOP1 query time:", tB, "seconds")
```

```{code-cell} ipython3
results
```

In order to use TAP with this ADQL string using pyvo, you can do the following:

```{code-cell} ipython3
# Uncomment the next line to run the query. Beware that it can take awhile.
# service.run_async(adql)
```

The above query shows that there are 597,488,849 redshifts in this table.
Running ``count`` on an entire table is an expensive operation, therefore we ran it asynchronously to avoid any potential timeout issues.
To learn more about synchronous versus asynchronous PyVO queries please read the [relevant PyVO documentation](https://pyvo.readthedocs.io/en/latest/dal/index.html#synchronous-vs-asynchronous-query).

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

## Get a list of galaxies within a small area

Since we know that cosmoDC2 is a large catalog, we can start with a spatial search over a small square area. The ADQL that is needed for the spatial constraint is shown below.  We then show how to make a redshift histogram of the sample generated.

```{code-cell} ipython3
# Setup the query
adql = f"""
SELECT TOP 100 redshift
FROM {tablename}
WHERE CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', 54.2, -37.4, 0.05)
) = 1
"""
job = service.submit_job(adql)
job.run()
if wait_for_job(job, timeout=100):
    spatial_results = job.fetch_result()
    print("Rows:", len(results))
else:
    print("Job did NOT finish — cancelled or failed.")
```

```{code-cell} ipython3
spatial_results
```

```{code-cell} ipython3
if spatial_results:
    ]# Plot a histogram
    num_bins = 20
    # the histogram of the data
    n, bins, patches = plt.hist(spatial_results['redshift'], num_bins,
                                facecolor='blue', alpha = 0.5)
    plt.xlabel('Redshift')
    plt.ylabel('Number')
    plt.title('Redshift Histogram CosmoDC2 Mock Catalog V1 abridged')
```

We can see form this plot that the simulated galaxies go out to z = 3.

+++

## Visualize galaxy colors: redshift search

First, we'll do a narrow redshift cut with no spatial constraint.  Then, from that redshift sample we will visualize the galaxy main sequence at z = 2.0.

```{code-cell} ipython3
# Setup the query
adql = f"""
SELECT TOP 50000
    Mag_r_LSST,
    Mag_g_LSST,
    redshift
FROM {tablename}
WHERE redshift > 1.95 and redshift < 2.05

"""
# Run the query
job = service.submit_job(adql)
job.run()

#if the job does not finish in a reasonable amount of time, cancel it
if wait_for_job(job, timeout=1000):
    redshift_results = job.fetch_result()
    print("Rows:", len(results))
else:
    print("Job did NOT finish — cancelled or failed.")
```

```{code-cell} ipython3
redshift_results
```

```{code-cell} ipython3
if redshift_results:
    # Construct a 2D histogram of the galaxy colors
    plt.hist2d(redshift_results['mag_r_lsst'], redshift_results['mag_g_lsst']-redshift_results['mag_r_lsst'],
               bins=100, cmap='plasma', cmax=500)

    # Plot a colorbar with label.
    cb = plt.colorbar()
    cb.set_label('Number')

    # Add title and labels to plot.
    plt.xlabel('LSST Mag r')
    plt.ylabel('LSST rest-frame g-r color')

    # Show the plot.
    plt.show()
```

+++ {"jp-MarkdownHeadingCollapsed": true}

## Suggestions for further queries:
TAP queries are extremely powerful and provide flexible ways to explore large catalogs like CosmoDC2, including spatial searches, photometric selections, cross-matching, and more. However, many valid ADQL queries can take minutes or longer to complete due to the size of the catalog, so we avoid running those directly in this tutorial. Instead, the examples here have so far focused on fast, lightweight queries that illustrate the key concepts without long wait times. If you are interested in exploring further, here are some additional query ideas that are scientifically useful but may take longer to run depending on server conditions.

### How many redshifts are in the chosen table?
`adql = f"SELECT count(redshift) FROM {tablename}"  #answer: 597,488,849 redshifts`

### Retrieve only a subset of columns (recommended for speed)
This use of "TOP 5000" just limits the number of rows returned. Remove it if you want all rows

`adql = f"SELECT TOP 5000 ra, dec, redshift, stellar_mass FROM {tablename}"`

### Cone search around a specific position
This search is slower than the spatial search above because it uses "contains" which does not take advantage of position indexing.

`adql = f""" SELECT TOP 50000 redshift FROM {tablename} WHERE CONTAINS(POINT('ICRS', RAMean, DecMean), CIRCLE('ICRS',54.2, -37.5,.1))=1`

+++

***

## About this notebook

**Author:** IRSA Data Science Team, including Vandana Desai, Jessica Krick, Troy Raen, Brigitta Sipőcz, Andreas Faisst, Jaladh Singhal

**Updated:** December 2025

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
