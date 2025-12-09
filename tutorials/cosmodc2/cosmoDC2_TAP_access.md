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

# Querying CosmoDC2 Mock v1 catalogs

This tutorial demonstrates how to access the CosmoDC2 Mock V1 catalogs. More information about these catalogs can be found here: https://irsa.ipac.caltech.edu/Missions/cosmodc2.html

These catalogs can be accessed through IRSA's Virtual Ovservatory Table Access Protocol (TAP) service. See https://www.ivoa.net/documents/TAP/ for details on the protocol. This service can be accessed through Python using the PyVO library.

Tips:
- This catalog is spatially indexed so searching on position will be the fastest searches.
- If searching on position, make sure to search in the area that cosmoDC3 covers.
  - CosmoDC2 covers an area of roughly: RA ≈ 0° to 60°, Dec ≈ –45° to 0°, but not exactly
- run queries in a way that can be cancelled.  It is very easy to overload the TAP server
  with multiple queries and not be able to cancel them if you just do a run_sync.

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy matplotlib pyvo
```

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
            try:
                job.abort()
            except Exception as e:
                print("Abort raised exception:", e)
            return False

        time.sleep(poll_interval)
```

```{code-cell} ipython3
test_small = "SELECT TOP 1 * FROM tap_schema.tables"

t0 = time.time()
job = service.submit_job(test_small)
print("Submit time:", time.time() - t0)
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
print("Rows returned:", len(resA))
```

```{code-cell} ipython3
# ------- Test B: RA/Dec CONTAINS query (index test) -------
test_adql_B = f"""
SELECT TOP 1 redshift
FROM {tablename}
WHERE CONTAINS(
    POINT('ICRS', RAMean, DecMean),
    CIRCLE('ICRS', 30.0, -36.0, 0.2)
) = 1
"""

jobB = service.submit_job(test_adql_B)

t0 = time.time()
jobB.run()

if wait_for_job(jobB, timeout=20):
    results = jobB.fetch_result()
    print("Rows:", len(results))
else:
    print("Job did NOT finish — cancelled or failed.")

tB = time.time() - t0

print("Cone-search TOP1 query time:", tB, "seconds")
print("Rows returned:", len(resB))
```

```{code-cell} ipython3
print(jobB.phase)
```

## How many rows are in the chosen table?

With TAP, you can query catalogs with constraints specified in IVOA Astronomical Data Query Language (ADQL; https://www.ivoa.net/documents/latest/ADQL.html), which is based on SQL.

```{code-cell} ipython3
# For example, this snippet of ADQL counts the number of elements in
# the redshift column of the table we chose.
adql = f"SELECT count(redshift) FROM {tablename}"
adql
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

## Spatial Search

Since we found out above that cosmoDC2 a large catalog, we can start with a spatial search over a small area of 0.1 deg. The ADQL that is needed for the spatial constraint is shown below.  We then show how to make a redshift histogram of the sample generated.

```{code-cell} ipython3
adql = f"SELECT redshift FROM {tablename} WHERE CONTAINS(POINT('ICRS', RAMean, DecMean), CIRCLE('ICRS',54.218205903,-37.497959343,.1))=1"
adql
```

Now we can use the previously-defined service to execute the query with the spatial contraint.

```{code-cell} ipython3
cone_results = service.run_sync(adql)
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

## Visualize galaxy colors at z ~ 2.0

Now let's visualize the galaxy main sequence at z = 2.0. First, we'll do a narrow redshift cut with no spatial constraint.

Let's do it as an asynchronous search since this might take awhile, too.

```{code-cell} ipython3
adql = f"SELECT TOP 5 redshift FROM {tablename}"
results = service.run_sync(adql)
```

```{code-cell} ipython3
service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
result = service.run_sync("SELECT TOP 1 * FROM tap_schema.tables")
```

```{code-cell} ipython3
result = service.run_sync("SELECT TOP 5 redshift FROM cosmodc2mockv1_heavy")
```

```{code-cell} ipython3
service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
jobs = service.get_jobs()
jobs
```

```{code-cell} ipython3
service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")

#setup the query
adql = f"""
SELECT TOP 50000
    Mag_true_r_sdss_z0,
    Mag_true_g_sdss_z0,
    redshift
FROM {tablename}
WHERE CONTAINS(
    POINT('ICRS', RAMean, DecMean),
    CIRCLE('ICRS', 0.5, -43.0, 0.2)
    ) = 1

"""
#run the query
results = service.run_sync(adql)
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

***

## About this notebook

**Author:** Vandana Desai (IRSA Science Lead)

**Updated:** 2024-07-24

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
