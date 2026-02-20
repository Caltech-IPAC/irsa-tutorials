# SPHEREx Source Discovery Tool

## Overview

The SPHEREx Source Discovery Tool is the Python package `spx_sdt` (included in this directory), which is used to discover and extract sources from SPHEREx Spectral Images and visualize their spectra.
The notebook [](sdt_irsa.md) demonstrates how to use it.

### Directory Contents

- `sdt_irsa.md`: Demo notebook.
- `conda-sdt_env.yml`: Conda environment definition file.
- `spx_sdt/`: Python package, including functions and configuration files to extract and visualize sources.
  - Python modules
    - `aperture_photometry.py`: Functions adapted from Zafar Rustamkulov to perform aperture photometry.
    - `bokeh_viz.py`: Functions for dynamic `bokeh` visualizations.
    - `firefly_viz.py`: Functions for `Firefly` visualizations.
    - `source_extraction.py`: Functions to work with SExtractor.
    - `sdt_utils.py`: Miscellaneous utility functions for working with IRSA tables and SPHEREx image headers.
  - SExtractor files
    - `default.conv`: Convolution mask used as detection filter.
    - `default.nnw`: Table of neural-network weights for star/galaxy separation.
    - `default_sdt.sex`: SExtractor default config file. Specifies catalog, extraction, photometry, etc. options.
    - `default_sdt.param`: SExtractor default parameters file. Specifies columns to include in output SExtractor catalog.

## Setup

### On local machine

All necessary packages and tools (SExtractor) are included in the `conda-sdt_env.yml` file.

To create a conda environment with the dependencies on your local machine use:
  ```
  conda env create --f conda-sdt_conda.yml
  ```

Then, make the environment available in the list of kernels for your JupyterLab:
  ```
  conda activate sdt_env
  python -m ipykernel install --user --name=sdt_env
  ```

To use the environment in your Jupyter notebooks, either start JupyterLab in that environment by typing

```
conda activate std_env
jupyter-lab
```

or select the environment `sdt_env` in your Jupyter Notebook using the dropdown on the upper left.

### On Fornax

Installing the conda environment on the [NASA Fornax Science Console](https://science.nasa.gov/astrophysics/programs/physics-of-the-cosmos/community/the-fornax-initiative/) needs slightly different steps. These can be reviewed in the documentation [create a new environment](https://docs.fornax.sciencecloud.nasa.gov/compute-environments/#create-new-env).

In order to install this specific conda environment on Fornax, the file name of the `yml` file specifically needs to be in the following format `conda-*.yml`. The `yml` file distributed here is already in that format (`conda-sdt_env.yml`). Once this is set, open a new terminal (click on the large "+" button right under the "File" menu tab) and type to following command _inside_ the same directory where the `yml` file is located:
  ```
  setup-conda-env --user
  ```

Note that we use the `--user` option here, which will keep the environment available for subsequent Fornax sessions.

To use the environment in your Jupyter notebooks on Fornax, directly select the environment `sdt_env` in your Jupyter Notebook using the dropdown on the upper left.
