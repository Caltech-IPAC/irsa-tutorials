---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Introduction to Euclid Q1 MER catalog

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- Understand the basic characteristics of Euclid Q1 MER catalogs.
- What columns are available in the MER catalog.
- How to query with ADQL in the MER catalog.
- How to make a simple color-magnitude diagram with the data.

+++

## Introduction

+++

Euclid is a European Space Agency (ESA) space mission with NASA participation, to study the geometry and nature of the dark Universe.
The Quick Data Release 1 (Q1) are the first data release from the Euclid mission after the Early Release Observations (ERO).
On March 19, 2025 the data will be available on the [ESA archive](https://easidr.esac.esa.int/sas/) and on the [IRSA archive](https://irsa.ipac.caltech.edu).

These Q1 notebooks focus on how to access, download, and process Euclid Q1 data from the IRSA archive.
If you have any issues accessing data from the archives, please contact the helpdesk directly: [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) and [ESA Euclid Helpdesk](https://support.cosmos.esa.int/euclid).

Each entry in the MER catalog is a single source containing all its photometry from the MER Mosaics (VIS, Y, J, H and any accompanying external ground observations) along with other basic measurements, like size and shape.

This notebook provides an introduction to the MER catalog released as part of Euclid Q1.
Other Euclid notebooks show how to use other data products released as part of Euclid Q1.

+++

## Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed
# !pip install numpy matplotlib pyvo
```

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

import pyvo as vo
```

## 1. Download MER catalog from IRSA directly to this notebook

```{code-cell} ipython3
service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
```

```{code-cell} ipython3
tables = service.tables
for tablename in tables.keys():
    if "tap_schema" not in tablename and "euclid_q1" in tablename:
            tables[tablename].describe()
```

### Choose the Euclid MER table

```{code-cell} ipython3
table_mer = 'euclid_q1_mer_catalogue'
```

### Learn some information about the table:
- How many columns are there?
- List the column names

```{code-cell} ipython3
columns = tables[table_mer].columns
print(len(columns))
```

```{code-cell} ipython3
for col in columns:
    print(f'{f"{col.name}":30s}  {col.unit}  {col.description}')
```

### Define the following ADQL query to find the first 10k stars in the MER catalog

Since we are just using the MER catalog alone, it does not have a column for classification. We can use the point_like_flag = 1 or point_like_prob>0.99 for stars.

Set all the fluxes to be greater than 0 so the object is detected in all four Euclid MER mosaic images

```{code-cell} ipython3
adql_stars = ("SELECT TOP 10000 mer.ra, mer.dec, mer.flux_vis_psf, mer.fluxerr_vis_psf, mer.flux_y_templfit,mer.fluxerr_y_templfit, "
    "mer.flux_j_templfit, mer.fluxerr_j_templfit, mer.flux_h_templfit, mer.fluxerr_h_templfit, mer.point_like_prob, mer.extended_prob "
    f"FROM {table_mer} AS mer "
    "WHERE  mer.flux_vis_psf > 0 "
    "AND mer.flux_y_templfit > 0 "
    "AND mer.flux_j_templfit > 0 "
    "AND mer.flux_h_templfit > 0 "
    "AND mer.point_like_flag = 1 ")

# Run the query

result_stars = service.search(adql_stars)
```

```{code-cell} ipython3
df_s_irsa = result_stars.to_table().to_pandas()   # Convert to Pandas DataFrame

# Display first few rows
df_s_irsa.head()
```

## 2. Make a color-magnitude diagram using the catalogs pulled from IRSA

- Convert from flux in uJy to magnitudes using the zero point correction
- Convert the error bars to magnitudes as well
- Plot the color-magnitude diagram

```{code-cell} ipython3
mag_y_s_irsa=-2.5*np.log10(df_s_irsa["flux_y_templfit"]) + 23.9 # Y
mag_h_s_irsa=-2.5*np.log10(df_s_irsa["flux_h_templfit"]) + 23.9 # H

x_s_irsa = mag_y_s_irsa - mag_h_s_irsa # Y - H
y_s_irsa = mag_y_s_irsa

xerr_s_irsa= 2.5 / np.log(10) * np.sqrt((df_s_irsa["fluxerr_y_templfit"] / df_s_irsa["flux_y_templfit"])**2
                                 + (df_s_irsa["fluxerr_h_templfit"] / df_s_irsa["flux_h_templfit"])**2)
yerr_s_irsa= 2.5 / np.log(10) * (df_s_irsa["fluxerr_y_templfit"] / df_s_irsa["flux_y_templfit"])

plt.errorbar(x_s_irsa, y_s_irsa, xerr=xerr_s_irsa, yerr=yerr_s_irsa, fmt='o', markersize=1.5, ecolor='lightgrey', elinewidth=0.5, capsize=2)

plt.xlabel('Y-H')
plt.ylabel('Y')
plt.xlim(-10,10)
plt.ylim(10,35)
plt.title('10k Stars in MER catalog -- IRSA')
```

## About this Notebook

**Author**: Tiffany Meshkat (IPAC Scientist)

**Updated**: 2025-03-19

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
