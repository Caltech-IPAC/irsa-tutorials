----
short_title: "Euclid"
----

# Euclid Tutorial Notebooks

[Euclid](https://irsa.ipac.caltech.edu/Missions/euclid.html) launched in July 2023 with the primary science goals of better understanding the composition and evolution of the dark Universe.
It carries two instruments: the VISible instrument (VIS) and the Near-Infrared Spectrometer and Photometer (NISP).

[Quick Release 1 (Q1)](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_q1.html) was released in March 2025 and consists of approximately 35 TB of imaging, spectroscopy, and catalogs, covering four non-contiguous fields totaling 63 square degrees.
Data products include MERged mosaics of calibrated and stacked frames; combined infrared spectra (SIR); and catalogs of MER objects, photometric redshifts and classifications (PHZ), and spectroscopic redshifts and line measurements (SPE).

## Images

- [MER Mosaics](1_Euclid_intro_MER_images.md) — Retrieve both a full MER mosaic image and multi-wavelength cutouts, then subtract the background from the cutouts and extract sources.

## Spectra

- [SIR 1D Spectra](3_Euclid_intro_1D_spectra.md) — Load a galaxy spectrum and plot it. Understand the wavelength, flux, and mask values.

## Catalogs

- [MER Catalogs](2_Euclid_intro_MER_catalog.md) — Explore the columns in the MER final catalog, query for stars, and create a color-magnitude diagram.
- [PHZ Catalogs](4_Euclid_intro_PHZ_catalog.md) — Join the PHZ and MER catalogs and do a box search for galaxies with quality redshifts, load a MER mosaic cutout of the box, and plot the cutout with the catalog results overlaid.
  Then plot the SIR spectrum of the brightest galaxy and look at a MER mosaic cutout of the galaxy in Firefly.
- [SPE Catalogs](5_Euclid_intro_SPE_catalog.md) — Join the SPE and MER catalogs and query for galaxies with H-alpha line detections, then plot the SIR spectrum of a galaxy with a high SNR H-alpha line measurement.
- [Merged Objects HATS Catalog](../parquet-catalog-demos/euclid-q1-hats/1-euclid-q1-hats-intro.md) — Understand the content and format of the Euclid Q1 Merged Objects HATS Catalog, then perform a basic query.

## Special Topics

- [Cloud Access](../cloud_access/euclid-cloud-access.md) — Browse the on-cloud copy of Q1, then efficiently retrieve a MER mosaic cutout and a SIR spectrum.
- [Early Release Observations (ERO) - Star Clusters](Euclid_ERO.md) (deprecated; ERO is superseded by Q1) — Create multi-wavelength ERO image cutouts of a globular cluster, extract sources, and measure photometry.
  Then load Gaia sources at the location of the globular cluster, match with Euclid ERO catalogs, and visualize the results with Firefly.
