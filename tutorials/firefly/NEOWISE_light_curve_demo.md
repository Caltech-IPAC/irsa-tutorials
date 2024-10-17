---
jupytext:
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

# Using Firefly visualization tools to understand the light curves of Solar System objects

## Learning Goals

By the end of this tutorial, you will:

- Construct a TAP query to download the necessary data and visualize it via the web browser with an instantiated Firefly environment.
- Plot light curves from NEOWISE data using the Firefly Python API.
- Format cells containing tables, charts and images viewed in the client.
- Overlay a catalog of data onto a HiPS image.

+++

## Introduction

This tutorial demonstrates how to plot light curves from NEOWISE data while also flaunting the many useful capabilities of the Firefly Python API. Using the 'Known Solar System Object Possible Association List' catalog from the NEOWISE-R database, we can easily compose a light curve of the faint asteroid `558 Carmen` and show its observed positions in the sky solely through the `firefly_client` package. Minor planet light curves are important in determining their size, spectral class, rotation period and many other properties.

Firefly is an astronomical data access and visualization written by Caltech/IPAC-IRSA. The visualization provides user with an integrated experience with brushing and linking capabilities among images, catalogs, and plots. Firefly is used in IRSA GUIs to query and visualize data from missions such as WISE, Spitzer, SOFIA, ZTF, PTF, etc. and a large number of highly-used contributed data products from a diverse set of astrophysics projects.

The `firefly_client` package provides a lightweight client class that includes a Python interface to Firefly’s Javascript API.

For documentation on the firefly client visit https://caltech-ipac.github.io/firefly_client/.

+++

## Imports

- *firefly_client FireflyClient* - Python API to Firefly for displaying tables, images and charts
- *astropy.utils.data* for downloading the catalog data via TAP query
- *urllib.parse* for converting regular query string to url-safe string

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install firefly_client astropy
```

```{code-cell} ipython3
from firefly_client import FireflyClient
import astropy.utils.data
import urllib.parse
```

## Step 1

Instantiate the client and view it in a different tab in the web browser.

In this example, we use the IRSA Viewer - a public firefly server. The firefly server can also be run locally, e.g. via a Firefly Docker image obtained from https://hub.docker.com/r/ipac/firefly/tags/.

```{code-cell} ipython3
url = 'https://irsa.ipac.caltech.edu/irsaviewer'
# url='http://127.0.0.1:8080/firefly'  # if you have firefly server running locally (preferably through docker)

fc = FireflyClient.make_client(url)
```

You can re-initizialize the viewer to return to a clean state with [`reinit_viewer`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.reinit_viewer).

```{code-cell} ipython3
# fc.reinit_viewer();  # The semi-colon suppresses the output of the method when ran
```

## Step 2

TAP search the 'Known Solar System Object Possible Association List' catalog from the NEOWISE-R database. The specific target we are looking for is minor planet `558 Carmen`. We can query this target using a TAP search through IRSA; the `table_url` is broken down as follows:

- We want to search the data through IRSA, which supports TAP querying, and we want it streamed directly to us via a synchronous search: <br>"https://<!---->irsa.ipac.caltech.edu/TAP/sync"<br><br>
- Next, we want to structure query to only retrieve (558) Carmen data from the NEOWISE-R 'Known Solar System Object Possible Association List' catalog. The table name of the catalog can be found using [IRSAViewer](https://irsa.ipac.caltech.edu/irsaviewer/?__action=layout.showDropDown&view=MultiTableSearchCmd) and clicking the **VO TAP Search** tab and changing the 'Project' to **neowiser**. We query all columns of data and we search the target by its object id, which is its name, and use the 'like' condition to only write (558) with a wildcard as shown in the cell below.

Construction of the query can be found in the [`IRSA TAP documentation page`](https://irsa.ipac.caltech.edu/docs/program_interface/TAP.html).

```{code-cell} ipython3
BASE_URL = "https://irsa.ipac.caltech.edu/TAP/sync"
QUERY = """
SELECT * 
FROM neowiser_p1ba_mch AS n 
WHERE n.objid LIKE '(558)%'
"""
table_url = f"{BASE_URL}?QUERY={urllib.parse.quote_plus(QUERY)}"
table_url
```

Now, we can request the necessary data from the catalog and display the data as a table in the Firefly client, using the [`show_table`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_table) method.

Alternatively, we can download data from the catalog using [`astropy.utils.data.download_file`](https://docs.astropy.org/en/stable/api/astropy.utils.data.download_file.html) and upload it to the Firefly client shown in the cell below the first method.

```{code-cell} ipython3
fc.show_table(table_url, tbl_id='tableneo', title='558 Carmen NeoWise Catalog', page_size=50)
```

```{code-cell} ipython3
# tablename = astropy.utils.data.download_file(table_url, timeout=120, cache=True)
# file = fc.upload_file(tablename)
# fc.show_table(file, tbl_id='tableneo', title='558 Carmen Catalog', page_size=50)
```

Note that along with the table, firefly also displays the coverage and chart associated with the table. It overlays colored squares for each row of the table onto a HiPS image, because the table contains recognizable celestial coordinates. It also creates a scatter plot of ra and dec from the table.

+++

## Step 3

After retrieving the data and displaying it in the client, we can now create a light curve by plotting the Modified Julian Date ('mjd') in the abscissa and the magnitude from band W1 ('w1mpro') in the ordinate. We also flip the ordinate to accurately display magnitude.

```{code-cell} ipython3
fc.show_xyplot(tbl_id='tableneo', xCol='mjd', yCol='w1mpro', yOptions='flip')
```

## Step 4

Finally, we can overlay the catalog of data in the table onto a HiPS image of our choice, using the method [`show_hips`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_hips). However, this method requires target coordinates for the object you want to analyze.

```{code-cell} ipython3
target='229.851396;-9.720647;EQ_J2000'
viewer_id = 'hipsDiv'
hips_url = 'http://alasky.u-strasbg.fr/AllWISE/RGB-W4-W2-W1'

fc.show_hips(viewer_id=viewer_id, plot_id='aHipsID1-1', hips_root_url = hips_url,
                          Title='HiPS-WISE', WorldPt=target)
```

## Summary

Firefly allows you to visualize data for specific targets. In conjuction with Astropy, one can manipulate a catalog of observations to display a light curve in an instantiated Firefly enviroment on their web browser.

1. We import all necessary modules to create a Firefly client and to download the catalog of data for our target.

2. We start the client in our web browser to appropiately display our tables, plots and images.

3. We use the TAP schema to display the data for our target &mdash; [`558 Carmen`](https://irsa.ipac.caltech.edu/irsaviewer/?__action=table.search&request=%7B%22startIdx%22%3A0%2C%22SearchMethod%22%3A%22AllSky%22%2C%22RequestedDataSet%22%3A%22NEOWISE%20Reactivation%20Database%22%2C%22id%22%3A%22GatorQuery%22%2C%22tbl_id%22%3A%22tbl_id-cf48-45%22%2C%22META_INFO%22%3A%7B%22title%22%3A%22WISE-neowiser_p1ba_mch%20(AllSky)%22%2C%22tbl_id%22%3A%22tbl_id-cf48-45%22%2C%22tbl_pref_key%22%3A%22WISE-neowiser_p1ba_mch%22%7D%2C%22catalogProject%22%3A%22WISE%22%2C%22catalog%22%3A%22neowiser_p1ba_mch%22%2C%22constraints%22%3A%22objid%20like%20%27%25(558)%20Carmen%25%27%22%2C%22pageSize%22%3A100%7D&options=%7B%22backgroundable%22%3Atrue%2C%22pageSize%22%3A100%7D) &mdash; via a table and visualize such data through charts.

4. We finally overlay the catalog onto a HiPS image to dynamically view where our target has been observed in space.

+++

***

## About This Notebook

+++

**Author:** Eric Bratton II (IRSA Scientist) in conjunction with the IRSA Science Team<br>
**Updated On:** 2024-10-17<br>
**Contact:** irsasupport@ipac.caltech.edu or  https://irsa.ipac.caltech.edu/docs/help_desk.html
