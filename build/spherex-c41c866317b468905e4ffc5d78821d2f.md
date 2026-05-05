# SPHEREx Tutorial Notebooks

[SPHEREx](https://irsa.ipac.caltech.edu/Missions/spherex.html) (Spectro-Photometer for the History of the Universe, Epoch of Reionization, and Ices Explorer) is a NASA space mission designed to perform the first all-sky spectral survey in the near-infrared.
SPHEREx observes the sky from roughly 0.75–5.0 µm using a single instrument that provides low-resolution spectroscopy (R ≈ 40–150) in hundreds of spectral channels for every point on the sky.
Its science goals span cosmology, galaxy evolution, and the interstellar medium, enabling measurements of large-scale structure, the cosmic history of star formation, and the distribution of key molecules and ices in the Milky Way and nearby galaxies.

SPHEREx data releases include weekly [Quick Release spectral image products](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/) (multi-extension FITS files containing calibrated near-infrared surface brightness, variance, flags, modeled backgrounds, PSFs, and wavelength WCS) along with ancillary calibration and metadata files such as gain matrices, dark current maps, solid angle pixel maps, and detailed spectral WCS products for each detector.

```{notebook-gallery} notebook_metadata.yml
tutorials/spherex/spherex_intro.md
tutorials/spherex/spherex_cutouts.md
tutorials/spherex/spherex_psf.md
tutorials/spherex/spherex_source_discovery/spherex_source_discovery_tool_demo.md
```
