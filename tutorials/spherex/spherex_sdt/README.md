# __Source Discovery Tool__

## Overview
- `spx_sdt`: package including functions to run notebook
  - SExtractor files
    - `default.conv`: Convolution mask used as detection filter.
    - `default.nnw`: Table of neural-network weights for star/galaxy separation.
    - `default_sdt.sex`: SExtractor default config file. Specifies catalog, extraction, photometry, etc. options.
    - `default_sdt.param`: SExtractor default parameters file. Specifies columns to include in output SExtractor catalog.
  - Python modules
    - `aperture_photometry.py`: Functions adapted from Zafar Rustamkulov to perform aperture photometry.
    - `bokeh_viz.py`: Functions for dynamic `bokeh` visualizations.
    - `firefly_viz.py`: Functions for `Firefly` visualizations.
    - `source_extraction.py`: Functions to work with SExtractor.
    - `sdt_utils.py`: Miscellaneous utility functions for working with IRSA tables and SPHEREx image headers.
- `conda-sdt_env.yml`: Conda environment definition file.
- `sdt_irsa.ipynb`: Demo notebook.

## Setup

#### On local machine

All necessary packages and tools (SExtractor) are included in the `conda-sdt_env.yml` file.

To create a conda environment with the dependencies on your local machine use:
  ```
  conda env create --f conda-sdt_conda.yml
  ```

Then, make the environment available in the list of kernels for your Jupyter lab:
  ```
  conda activate sdt_env
  python -m ipykernel install --user --name=sdt_env
  ```

To use the environment in your Jupyter Lab notebooks, either start Jupyter Lab in that environment by typing

```
conda activate std_env
jupyter-lab
```

or select the environment `sdt_env` in your Jupyter Notebook using the dropdown on the upper left.

#### On Fornax

Installing the conda environment on the [Fornax NASA science platform](https://science.nasa.gov/astrophysics/programs/physics-of-the-cosmos/community/the-fornax-initiative/) needs slightly different steps. These can be reviewed in [this documentation](https://docs.fornax.sciencecloud.nasa.gov/compute-environments/#create-new-env).

In order to install this specific conda environment on Fornax, the file name of the `yml` file specifically needs to be in the following format `conda-*.yml`. The `yml` file distributed here is already in that format (`conda-sdt_env.yml`). Once this is set, open a new terminal (click on the large "+" button right under the "File" menu tab) and type to following command _inside_ the same directory where the `yml` file is located:
  ```
  setup-conda-env --user
  ```

Note that we use the `--user` option here, which will keep the environment available for subsequent Fornax sessions.

To use the environment in your Jupyter Lab notebooks on Fornax, directly select the environment `sdt_env` in your Jupyter Notebook using the dropdown on the upper left.
