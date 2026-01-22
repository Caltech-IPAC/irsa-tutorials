---
short_title: "WISE"
---

# WISE Tutorial Notebooks

The Wide-field Infrared Survey Explorer ([WISE](https://irsa.ipac.caltech.edu/Missions/wise.html)) is a NASA infrared space telescope launched in December 2009 that performed a sensitive all-sky survey at 3.4, 4.6, 12, and 22 µm, cataloging hundreds of millions of stars, galaxies, and Solar System objects and enabling discoveries of cool brown dwarfs and luminous infrared galaxies. 
After exhausting its cryogen, the mission was repurposed as NEOWISE in 2013 to detect and characterize near-Earth asteroids and comets using the remaining infrared channels.


WISE and NEOWISE data are released publicly through the NASA/IPAC Infrared Science Archive (IRSA), including calibrated images, source catalogs, and single-exposure source tables that together enable multi-epoch photometry, light curves, and motion studies for a wide range of astrophysical and Solar System applications. 
Successive NEOWISE data releases (with annual updates) provide users with increasingly deep coverage and time-domain information across the infrared sky.

- [ALLWISE Images](sia_allwise_atlas.md) - Retrieves coadded AllWISE mid-infrared images via IRSA’s SIA service and make coordinate-based cutouts.

- [ALLWISE Catalog](wise-allwise-catalog-demo.md) - Demonstrates querying, filtering, and working with the AllWISE Source Catalog to explore mid-infrared source properties.

- [NEOWISE visualization](NEOWISE_light_curve_demo.md) - Visualizes and analyzes light curves of Solar System objects using NEOWISE data and Firefly.

- [NEOWISE source table](neowise-source-table-strategies.md) - Outlines efficient strategies for accessing and handling large NEOWISE single-exposure source tables

- [NEOWISE source table light curves](neowise-source-table-lightcurves.md) - Builds multi-epoch photometric light curves from NEOWISE single-exposure source table data for given target coordinates.