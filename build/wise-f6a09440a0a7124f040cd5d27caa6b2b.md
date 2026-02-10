# WISE Tutorial Notebooks

The Wide-field Infrared Survey Explorer ([WISE](https://irsa.ipac.caltech.edu/Missions/wise.html)) is a NASA infrared space telescope launched in December 2009 that performed a sensitive all-sky survey at 3.4, 4.6, 12, and 22 µm, cataloging hundreds of millions of stars, galaxies, and Solar System objects and enabling discoveries of cool brown dwarfs and luminous infrared galaxies. 
After exhausting its cryogen, the mission was repurposed as NEOWISE in 2013 to detect and characterize near-Earth asteroids and comets using the remaining infrared channels until the mission concluded in August 2024.


WISE and NEOWISE data are released publicly through the NASA/IPAC Infrared Science Archive (IRSA), including calibrated images, source catalogs, and single-exposure source tables that together enable multi-epoch photometry, light curves, and motion studies for a wide range of astrophysical and Solar System applications. 
Successive NEOWISE data releases, issued with annual updates, provided users with increasingly deep coverage and time-domain information across the infrared sky.

- [AllWISE Images](irsa-sia-examples/sia_allwise_atlas.md) - Retrieve coadded images and make coordinate-based cutouts.

- [AllWISE Catalog](parquet-catalog-demos/wise-allwise-catalog-demo.md) - Querying, filter, and work with a HEALPix-partitioned Parquet catalog using the AllWISE dataset.

- [NEOWISE Firefly](firefly/NEOWISE_light_curve_demo.md) - Visualize and analyze light curves of Solar System objects using Firefly.

- [NEOWISE Strategies](parquet-catalog-demos/neowise-source-table-strategies.md) - Use efficient strategies for accessing and handling the very large Parquet dataset.

- [NEOWISE Light Curves](parquet-catalog-demos/neowise-source-table-lightcurves.md) - Build multi-epoch photometric light curves for given coordinates.