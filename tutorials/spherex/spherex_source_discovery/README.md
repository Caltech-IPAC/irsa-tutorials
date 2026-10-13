# SPHEREx Source Discovery Tool

## Overview

The SPHEREx Source Discovery Tool is the Python package `spherex_source_discovery_tool` (included in this directory), which is used to discover and extract sources from SPHEREx Spectral Images and visualize their spectra.
The notebook [](spherex_source_discovery_tool_demo.md) demonstrates how to use it.

### Directory Contents

- `spherex_source_discovery_tool_demo.md`: Demo notebook.
- `conda-spherex_sdt.yml`: Conda environment definition file.
- `spherex_source_discovery_tool/`: Python package, including functions and configuration files to extract and visualize sources.
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

All necessary packages and tools (SExtractor) are listed in the `conda-spherex_sdt.yml` file.

To create a conda environment with the dependencies on your local machine use:
  ```
  conda env create --file conda-spherex_sdt.yml
  ```

Then, make the environment available in the list of kernels for your JupyterLab:
  ```
  conda activate spherex_sdt
  python -m ipykernel install --user --name=spherex_sdt
  ```

To use the environment in your Jupyter notebooks, either start JupyterLab in that environment by typing

```
conda activate spherex_sdt
jupyter-lab
```

or select the environment `spherex_sdt` in your Jupyter Notebook using the dropdown on the upper left.

### On Fornax

Installing the conda environment on the [NASA Fornax Science Console](https://science.nasa.gov/astrophysics/programs/physics-of-the-cosmos/community/the-fornax-initiative/) needs slightly different steps. These can be reviewed in the documentation [create a new environment](https://docs.fornax.sciencecloud.nasa.gov/compute-environments/#create-new-env).

In order to install this specific conda environment on Fornax, the file name of the `yml` file specifically needs to be in the following format `conda-*.yml`. The `yml` file distributed here is already in that format (`conda-spherex_sdt.yml`). Once this is set, open a new terminal (click on the large "+" button right under the "File" menu tab) and type to following command _inside_ the same directory where the `yml` file is located:
  ```
  setup-conda-env --user
  ```

Note that we use the `--user` option here, which will keep the environment available for subsequent Fornax sessions.

To use the environment in your Jupyter notebooks on Fornax, directly select the environment `spherex_sdt` in your Jupyter Notebook using the dropdown on the upper left.
