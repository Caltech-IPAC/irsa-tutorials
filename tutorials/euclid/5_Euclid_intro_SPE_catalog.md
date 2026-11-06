---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  name: python3
  display_name: python3
  language: python
authors:
  - name: Tiffany Meshkat
  - name: Anahita Alavi
  - name: Anastasia Laity
  - name: Andreas Faisst
  - name: Brigitta Sipőcz
  - name: Dan Masters
  - name: Harry Teplitz
  - name: Jaladh Singhal
  - name: Shoubaneh Hemmati
  - name: Vandana Desai
  - name: Troy Raen
  - name: Jessica Krick
---

# Euclid DR1: SPE Catalogs

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- Understand the basic characteristics of Euclid DR1 SPE catalogs.
- Understand what SPE catalogs are available and how to view the columns in those catalogs.
- How to query with ADQL in the SPE lines catalog to find strong H-alpha detections.
- How to make a plot the detected line features over the 1D spectra.

+++

## Introduction

+++
Euclid launched in July 2023 as a European Space Agency (ESA) mission with contributions from NASA.
The primary science goals of Euclid are to better understand the composition and evolution of the dark Universe.
The Euclid mission provides space-based imaging and spectroscopy, complemented by ground-based imaging, to achieve these goals.
These data will be archived by multiple global repositories, including IRSA, where they will support a wide range of astrophysical research.

Euclid Data Release 1 (DR1) will include imaging, spectroscopy, and catalogs covering approximately 1,900 square degrees of the Euclid Wide Survey, along with observations of the Euclid Deep Fields. Further information about Euclid DR1 can be found at: https://euclid.caltech.edu/page/euclid-dr1

Among the DR1 data products are multiple spectroscopic catalogs produced by the Spectroscopic Analysis (SPE) Processing Function, including spectroscopic redshifts, emission-line measurements, and other derived spectral properties.

This notebook provides an introduction to the Euclid DR1 SPE catalogs.
If you have questions about this notebook, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## Imports

```{important}
We rely on ``astroquery`` features that have been recently added, so please make sure you have version v0.4.10 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed
# %pip install matplotlib astropy 'astroquery>=0.4.11'
```

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import QTable
from astropy import units as u
from astropy.utils.data import conf, download_file
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch, quantity_support

from astroquery.ipac.irsa import Irsa

# The Euclid spectrum files are large and the time it takes to read
# them can exceed astropy's default timeout limit. Increase it.
conf.remote_timeout = 120

#suppress warnings about deprecated units and cache
import warnings
warnings.filterwarnings(
    "ignore",
    message="The unit 'Angstrom' has been deprecated",
    category=u.UnitsWarning,
)

warnings.filterwarnings(
    "ignore",
    message="The unit 'erg' has been deprecated",
    category=u.UnitsWarning,
)
warnings.filterwarnings(
    "ignore",
    message="XDG_CACHE_HOME is set",
)
```

```{tip}
The IRSA SIA collections can be listed using using the ``list_collections`` method, we can filter on the ones containing "euclid" in the collection name:

    Irsa.list_collections(filter='euclid')
```



## 1. Download SPE catalog from IRSA directly to this notebook

Search for all Euclid catalogs and metadata tables in IRSA:

```{note}
In the example below, we use ``include_metadata_tables=True``.
This parameter was introduced in astroquery 0.4.11. For earlier versions, please omit this parameter to receive the combined list of available catalogs and metadata tables.
```

```{code-cell} ipython3
Irsa.list_catalogs(filter='euclid', include_metadata_tables=True)
```

```{code-cell} ipython3
table_mer = 'euclid_dr1_mer_catalogue_deep_survey'
table_galaxy_candidates = 'euclid_dr1_spectro_z_spe_galaxy_candidates_deep'
table_quality = 'euclid_dr1_spectro_z_spe_quality_deep'
table_lines = 'euclid_dr1_spectro_line_features_spe_line_features_deep'
table_photo_z = 'euclid_dr1_phz_photo_z_deep_survey'
```

### Learn some information about the table:
- How many columns are there?
- List the column names

```{code-cell} ipython3
columns_info = Irsa.list_columns(catalog=table_mer)
print(len(columns_info))
```

```{code-cell} ipython3
Irsa.list_columns(catalog=table_1dspectra, full=True)
```

```{code-cell} ipython3
# Full list of columns and their description
columns_info
```

## 2. Find objects within a specified sky region defined by RA and Dec ranges
Euclid DR1 covers a large sky area of approximately 1,900 square degrees. To reduce query execution time and avoid potential timeout issues, we recommend restricting searches to smaller sky regions.

As an example, let's select objects within a 0.1 square-degree region of the Euclid Deep Field North (EDF-N).

```{code-cell} ipython3
query_mer = f"""
SELECT DISTINCT
    object_id,
    ra,
    dec
FROM {table_mer}
WHERE
     ra BETWEEN 269.34391 AND 270.12192
AND dec BETWEEN 65.85958 AND 66.17581
"""

mer = Irsa.query_tap(query_mer).to_table()

object_ids = mer["object_id", "ra", "dec"].as_array().tolist()
object_ids = [row[0] for row in object_ids]
```


## 3. Find some objects with spectra with conditions  

We specify the following conditions on our search, based on recommendations in the Euclid DR1 spectroscopic paper:
- Signal to noise ratio column (_gf = gaussian fit) should be greater than 10
- We want to detect H-alpha
- We select objects with well-measured H-alpha flux based on the spe_line_qual_gf quality flags
- Choose spectroscopic redshift (spe_z) between 0.85 and 1.85 and spe_z_prob greater than 0.999
- We select objects whose spectroscopic redshifts (spe_z) fall within the 95% confidence intervals of their photometric redshifts (phz_95_int)
- Choose spectrum spe_npix and spe_n_dith_med greater than 15
- H-alpha line flux should be more than 2x10^16 erg s^-1 cm^-2
- We exclude point-like objects (primarily stars; mer.point_like_prob < 0.05), excessively large objects (mer.kron_radius < 50), and potentially spurious detections (mer.spurious_prob < 0.5). The size criterion helps minimize contamination from galaxies at $z < 0.85$ with incorrectly measured spectroscopic redshifts.
- Join the lines, mer, photometric redshift, and galaxy candidates tables on object_id and spe_rank

We limit our search to the first 10,000 objects to reduce query execution time and avoid potential timeout issues.
Finally we sort the data by descending spe_line_snr_gf to have the largest SNR H-alpha lines detected at the top.

```{code-cell} ipython3
batch_size = 500

line_tables = []
for i in range(0, 10000, batch_size):

    batch = object_ids[i:i + batch_size]

    ids = ", ".join(str(int(x)) for x in batch)

    query_lines = f"""
    SELECT
        lines.object_id,
        lines.spe_rank,
        lines.spe_line_name,
        lines.spe_line_qual_gf,
        lines.spe_line_central_wl_gf,
        lines.spe_line_snr_gf,
        galaxy.object_id,
        galaxy.spe_z
    FROM {table_lines} AS lines
    JOIN {table_galaxy_candidates} AS galaxy
    ON lines.object_id = galaxy.object_id
    AND lines.spe_rank = galaxy.spe_rank
    JOIN {table_quality} AS quality
    ON lines.object_id = quality.object_id
    JOIN {table_mer} AS mer
    ON lines.object_id = mer.object_id
    JOIN {table_photo_z} AS phz
    ON lines.object_id = phz.object_id
    
    WHERE
    lines.object_id IN ({ids})
    AND lines.spe_line_name = 'Halpha'
    AND lines.spe_line_qual_gf IN (212, 222, 252)
    AND lines.spe_line_snr_gf > 10
    AND lines.spe_rank = 0
    AND lines.spe_line_flux_gf > 2E-16
    AND galaxy.spe_z_prob > 0.99
    AND galaxy.spe_z_rel > 0.05
    AND galaxy.spe_rank = 0
    AND galaxy.spe_z > 0.85
    AND galaxy.spe_z < 1.85
    AND quality.spe_npix > 800
    AND quality.spe_n_dith_med > 15
    AND (mer.point_like_prob < 0.05) 
    AND (mer.spurious_prob < 0.5) 
    AND (mer.kron_radius < 50)
    AND (galaxy.spe_z > phz.phz_95_int_1)
    AND (galaxy.spe_z < phz.phz_95_int_2)
    ORDER BY lines.spe_line_snr_gf DESC
        """

    result = Irsa.query_tap(
        query_lines,
        async_job=True
    ).to_table()

    line_tables.append(result)

if line_tables:
    lines = vstack(line_tables)
else:
    lines = Table()


# --------------------------------------------------
# Add RA and Dec from the MER catalog
# --------------------------------------------------

if len(lines) > 0:
    lines = join(
        lines,
        mer["object_id", "ra", "dec"],
        keys="object_id",
        join_type="left"
    )

# Sort final table by SNR (highest first)
lines.sort("spe_line_snr_gf", reverse=True)
```

### Choose an object of interest, lets look at an object with a strong Halpha line detected with high SNR.

```{code-cell} ipython3
obj_id = 2700933273661649762

obj_row = result_table[(result_table['object_id'] == obj_id)]

obj_row
```

### Pull the spectrum of this object

```{code-cell} ipython3
# Query SSA for the 1D spectrum near this object's sky position
euclid_ssa_collection = "euclid_dr1_DpdSirCombinedSpectra"

# Use the object's MER coordinates from obj_row
coord_obj = SkyCoord(obj_row["ra"][0], obj_row["dec"][0], unit=u.deg)

#complete the query
ssa_result = Irsa.query_ssa(
    pos=coord_obj,
    radius=1.0 * u.arcsec,
    collection=euclid_ssa_collection,
)

ssa_result
```

### The following steps to read in the spectrum follows the 3_Euclid_intro_1D_spectra notebook.

```{code-cell} ipython3
# Read in the spectrum for this object from the SSA access URL
row = ssa_result[1]
row_blue = ssa_result[0]

spectrum_path = row["access_url"]
spectrum_path

spectrum_path_blue = row_blue["access_url"]
spectrum_path_blue
```

```{code-cell} ipython3
spectrum = QTable.read(spectrum_path)
spectrum_blue = QTable.read(spectrum_path_blue)
```

### Now the data are read in, plot the spectrum with the H-alpha line labeled

```{tip}
As we use astropy.visualization's ``quantity_support``, matplotlib automatically picks up the axis units from the quantities we plot.
```

```{code-cell} ipython3
quantity_support()
```

```{code-cell} ipython3
# Note that the units are missing from the lines table, we manually add Angstrom
bad_mask = (spectrum['MASK'].value % 2 == 1) | (spectrum['MASK'].value >= 64)
bad_mask_blue = (spectrum_blue['MASK'].value % 2 == 1) | (spectrum_blue['MASK'].value >= 64)

plt.plot(spectrum['WAVELENGTH'].to(u.micron), np.ma.masked_where(bad_mask, spectrum['SIGNAL']), alpha=0.6, color='red', label='Red Spectrum')
plt.plot(spectrum['WAVELENGTH'], np.ma.masked_where(~bad_mask, spectrum['SIGNAL']), color='red', label='Do not use (Red)', linestyle=':', linewidth=2.0)
plt.plot(spectrum['WAVELENGTH'], np.sqrt(spectrum['VAR']), color='grey', label='Error')
plt.plot(spectrum_blue['WAVELENGTH'].to(u.micron), np.ma.masked_where(bad_mask_blue, spectrum_blue['SIGNAL']), alpha=0.5, color='blue', label='Blue Spectrum')#, #color='black', label='Spectrum')
plt.plot(spectrum_blue['WAVELENGTH'], np.ma.masked_where(~bad_mask_blue, spectrum_blue['SIGNAL']), color='blue', label='Do not use (Blue)', linestyle=':', linewidth=2.0)
plt.plot(spectrum_blue['WAVELENGTH'], np.sqrt(spectrum_blue['VAR']), color='grey', label='Error')

plt.legend(loc='upper right')
plt.ylim(-0.05, 0.25)
plt.title(f'Object ID {obj_id}')
```

## About this Notebook

**Updated**: 2026-01-13

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.

**Runtime:** As of the date above, this notebook takes about 90 seconds to run to completion on a machine with 8GB RAM and 4 CPU.
