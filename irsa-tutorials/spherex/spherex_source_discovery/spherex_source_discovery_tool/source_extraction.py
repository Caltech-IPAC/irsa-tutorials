"""Extracts sources from Level 2 SPHEREx images.
"""

# Standard library imports
from pathlib import Path
from subprocess import run, CalledProcessError
import importlib.resources as resources

__all__ = ["run_sextractor", "get_sextractor_file"]


def run_sextractor(cut_filepath: str, sxt_config: str, sxt_params: str, sxt_cat_path: str, sxt_nnw: str,
                   sxt_conv: str):
    """Runs Source Extractor on a cutout .fits file.

    Parameters
    ----------
    cut_filepath : `str`
        Path to cutout .fits file.
    sxt_config : `str`
        Path to default.sex configuration file.
    sxt_params : `str`
        Path to default.param file.
    sxt_nnw: `str`
        Path to default.nnw file.
    sxt_conv: `str`
        Path to default.conv file.
    sxt_cat_path : `str`
        Path to output SExtractor catalog.

    Returns
    -------
    None.
    """

    try:
        run(["sex", cut_filepath,
             "-c", sxt_config,
             "-PARAMETERS_NAME", sxt_params,
             "-CATALOG_NAME", sxt_cat_path,
             "-STARNNW_NAME", sxt_nnw,
             "-FILTER_NAME", sxt_conv
             ], check=True)
    except CalledProcessError as e:
        print(e.stderr)


def get_sextractor_file(rel_filename: str, package: str = "spherex_source_discovery_tool") -> Path:
    """Returns path to package file.

    Parameters
    ----------
    rel_filename : `str`
        The filename.
    package : `str`, optional
        The package to look in, by default "spherex_source_discovery_tool".

    Returns
    -------
    `~pathlib.Path`
        The path to the package file.

    Raises
    ------
    FileNotFoundError
        Raised if the file is not found.
    """
    pkg = resources.files(package)
    resource = pkg / rel_filename
    if not resource.is_file():
        raise FileNotFoundError(f"{rel_filename} not found.")
    else:
        with resources.as_file(resource) as file_path:
            return file_path
