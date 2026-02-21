# SPHEREx Tutorial Notebooks

[SPHEREx](https://irsa.ipac.caltech.edu/Missions/spherex.html) (Spectro-Photometer for the History of the Universe, Epoch of Reionization, and Ices Explorer) is a NASA space mission designed to perform the first all-sky spectral survey in the near-infrared.
SPHEREx observes the sky from roughly 0.75–5.0 µm using a single instrument that provides low-resolution spectroscopy (R ≈ 40–150) in hundreds of spectral channels for every point on the sky.
Its science goals span cosmology, galaxy evolution, and the interstellar medium, enabling measurements of large-scale structure, the cosmic history of star formation, and the distribution of key molecules and ices in the Milky Way and nearby galaxies.

SPHEREx data releases include weekly [Quick Release spectral image products](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/) (multi-extension FITS files containing calibrated near-infrared surface brightness, variance, flags, modeled backgrounds, PSFs, and wavelength WCS) along with ancillary calibration and metadata files such as gain matrices, dark current maps, solid angle pixel maps, and detailed spectral WCS products for each detector.

- [Data Overview](spherex_intro.md) - Identify available data products and select the appropriate FITS extensions for different use cases.

- [Spectral Image Cutouts](spherex_cutouts.md) - Generate and work with spatial and spectral cutouts.

- [PSF Models](spherex_psf.md) - Understand how SPHEREx point spread function (PSF) information is organized and accessed.

- [Source Discovery Tool](spherex_sdt/sdt_irsa.md) - Discover, extract, and visualize sources from SPHEREx Spectral Images.
