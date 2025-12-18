----
short_title: "Quick Release 1 (Q1)"
----

(euclid-q1)=
# Euclid Quick Release 1 (Q1) Tutorial Notebooks

The tutorials in this section describe the [Euclid Q1](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_q1.html) data products and demonstrate access using various methods.

(euclid-q1-vo-series)=
## VO Access to Images, Spectra, and Catalogs

These tutorials use the Python library [astroquery.ipac.irsa](https://astroquery.readthedocs.io/en/latest/ipac/irsa/irsa.html) to access data through IRSA's Virtual Observatory (VO) services: Simple Image Access ([SIA](https://irsa.ipac.caltech.edu/docs/program_interface/api_images.html#siav2)) and Table Access Protocol ([TAP](https://irsa.ipac.caltech.edu/docs/program_interface/TAP.html)).
These methods are well-suited for most use cases, especially exploratory analyses, cone searches, and working with individual objects or small- to medium-sized subsets of data.

- [MER Mosaics](#euclid-q1-vo-mer-mosaics) — Retrieve both a full MER mosaic image and multi-wavelength cutouts, then subtract the background from the cutouts and extract sources.
- [MER Catalogs](#euclid-q1-vo-mer-catalogs) — Explore the columns in the MER final catalog, query for stars, and create a color-magnitude diagram.
- [SIR 1D Spectra](#euclid-q1-vo-sir-1d-spectra) — Load a galaxy spectrum and plot it. Understand the wavelength, flux, and mask values.
- [PHZ Catalogs](#euclid-q1-vo-phz-catalogs) — Join the main photometric redshift catalog with the MER final catalog and do a box search for galaxies with quality redshifts, load a MER mosaic cutout of the box, and plot the cutout with the catalog results overlaid.
  Then plot the SIR spectrum of the brightest galaxy and look at a MER mosaic cutout of the galaxy in Firefly.
- [SPE Catalogs](#euclid-q1-vo-spe-catalogs) — Join spectroscopic and MER catalogs and query for galaxies with H-alpha line detections, then plot the SIR spectrum of a galaxy with a high SNR H-alpha line measurement.

(euclid-q1-cloud-series)=
## Cloud Access to Images, Spectra, and Catalogs

IRSA also serves the Euclid Q1 data in FITS and Apache Parquet file formats from an Amazon Web Services (AWS) S3 cloud storage bucket.
This tutorial demonstrates access using the Python libraries [S3Fs](https://s3fs.readthedocs.io/en/latest/) and [Astropy](https://astropy.readthedocs.io/en/stable/).
S3 is well-suited for most use cases and is especially efficient for large-scale data access.

- [Cloud Access](#euclid-q1-cloud) — Browse the S3 bucket, then efficiently retrieve a MER mosaic cutout and a SIR spectrum.
- See [Merged Objects HATS Catalog](#euclid-q1-hats-series) for efficient access to the on-cloud catalogs.

(euclid-q1-hats-series)=
## Merged Objects HATS Catalog

The Euclid Q1 Merged Objects HATS Catalog comprises fourteen of the original Euclid Q1 catalogs, joined on the Euclid object ID, into a single dataset.
The catalog is in Apache Parquet file format and is stored on AWS S3.
These tutorials use the Python [Apache PyArrow](https://arrow.apache.org/docs/python/index.html) library for efficient data access.
HATS, Parquet, and PyArrow are described further on [IRSA Parquet Catalogs](https://irsa.ipac.caltech.edu/docs/parquet_catalogs/).
These methods are well-suited for most use cases and are especially good for large-scale analyses.

- [Introduction](#euclid-q1-hats-intro) — Understand the content and format of the Euclid Q1 Merged Objects HATS Catalog, then perform a basic query.
- [Redshifts](#euclid-q1-hats-redshifts) — Obtain quality samples of multiple redshift estimates, plot distributions and comparisons, and learn to choose appropriate redshifts for different use cases.
- [Classifications](#euclid-q1-hats-classifications) — blah, blah
- [Magnitudes](#euclid-q1-hats-magnitudes) — blah, blah
- [Morphology](#euclid-q1-hats-morphology) — blah, blah
- [NIR-only Objects](#euclid-q1-hats-nir-only) — blah, blah
