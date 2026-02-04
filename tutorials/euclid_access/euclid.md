# Euclid Tutorial Notebooks

[Euclid](https://irsa.ipac.caltech.edu/Missions/euclid.html) launched in July 2023 with the primary science goals of better understanding the composition and evolution of the dark Universe.
It carries two instruments: the VISible instrument (VIS) and the Near-Infrared Spectrometer and Photometer (NISP).

[Quick Release 1 (Q1)](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_q1.html) was released in March 2025 and consists of approximately 35 TB of imaging, spectroscopy, and catalogs, covering four non-contiguous fields totaling 63 square degrees.
Data products include MERged mosaics of calibrated and stacked frames; combined infrared spectra (SIR); and catalogs of MER objects, photometric redshifts and classifications (PHZ), and spectroscopic redshifts and line measurements (SPE).

Data products include:
- MER (merged) mosaic images of calibrated and stacked frames;
- SIR combined infrared spectra;
- Catalogs of MER objects, photometric redshifts and classifications (PHZ), and spectroscopic redshifts and line measurements (SPE);
- Merged Objects Catalog (created by IRSA) containing the MER, PHZ, and SPE catalogs in a single HATS Catalog.

````{grid} 1 2 2 3

```{card}
:link: 1_Euclid_intro_MER_images.md
:header: [MER Image Mosaics →](1_Euclid_intro_MER_images.md)
Retrieve both a full MER mosaic image and multi-wavelength cutouts, then subtract the background from the cutouts and extract sources.
```

```{card}
:link: 3_Euclid_intro_1D_spectra.md
:header: [SIR 1D Spectra →](3_Euclid_intro_1D_spectra.md)
Load a galaxy spectrum and plot it. Understand the wavelength, flux, and mask values.
```

```{card}
:link: 2_Euclid_intro_MER_catalog.md
:header: [MER Catalogs →](2_Euclid_intro_MER_catalog.md)

Explore the columns in the MER final catalog, query for stars, and create a color-magnitude diagram.
```

```{card}
:link: 4_Euclid_intro_PHZ_catalog.md
:header: [PHZ Catalogs →](4_Euclid_intro_PHZ_catalog.md)

Join the PHZ and MER catalogs to query galaxies with quality redshifts in a box region, create MER mosaic cutouts with catalog overlays, and plot the brightest galaxy's SIR spectrum.
```

```{card}
:link: 5_Euclid_intro_SPE_catalog.md
:header: [SPE Catalogs →](5_Euclid_intro_SPE_catalog.md)

Join the SPE and MER catalogs and query for galaxies with H-alpha line detections, then plot the SIR spectrum of a galaxy with a high SNR H-alpha line measurement.
```

```{card}
:link: ../parquet-catalog-demos/euclid-q1-hats/1-euclid-q1-hats-intro.md
:header: [Merged Objects Catalog →](../parquet-catalog-demos/euclid-q1-hats/1-euclid-q1-hats-intro.md)
Introduction: Understand the content and format of the Euclid Q1 Merged Objects HATS Catalog, then perform a basic query.
```

```{card}
:link: ../parquet-catalog-demos/euclid-q1-hats/4-euclid-q1-hats-magnitudes.md
:header: [Merged Objects Catalog →](../parquet-catalog-demos/euclid-q1-hats/4-euclid-q1-hats-magnitudes.md)
Magnitudes: Review the types of flux measurements available, load template-fit and aperture magnitudes, and plot distributions and comparisons for different object types.
```

```{card}
:link: ../cloud_access/euclid-cloud-access.md
:header: [Cloud Access →](../cloud_access/euclid-cloud-access.md)
Browse the on-cloud copy of Q1, then efficiently retrieve a MER mosaic cutout and a SIR spectrum.
```

```{card}
:link: Euclid_ERO.md
:header: [ERO Star Clusters →](Euclid_ERO.md)
Create multi-wavelength ERO image cutouts of a globular cluster, extract sources, and measure photometry. Match Gaia sources with Euclid ERO catalogs, then visualize with Firefly.
```

````
