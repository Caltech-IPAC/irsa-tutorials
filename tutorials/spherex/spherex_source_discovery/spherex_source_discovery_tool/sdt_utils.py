"""Utility functions to support working with SExtractor catalogs, IRSA
tables, and SPHEREx Level 2 images.
"""
# Standard library imports
import os
import psutil

# 3rd-party imports
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.units import Quantity

# Local imports

__all__ = ["format_extracted", "get_filename", "get_obs_id", "get_exp_id",
    "results_summary", "get_lambda_range"]

# Function definitions
def format_extracted(sxt_cat_filepath: str, sxt_cat_format: str, cut_filepath: str) -> Table:
    """Formats and converts SExtractor catalog data values.

    Parameters
    ----------
    sxt_cat_filepath : `str`
        Path to SExtractor output catalog file.
    sxt_cat_format : `str`
        Format of SExtractor output catalog file (see default.sex).
    cut_filepath : `str`
        Path to cutout .fits file.

    Returns
    -------
    sxt_tab : `astropy.table.Table`
        Table of extracted sources, including named IDs and converted fluxes.
    """
    # Read SExtractor catalog file
    sxt_tab = Table.read(sxt_cat_filepath, format=sxt_cat_format)

    # Format for SDT display table
    sxt_tab.rename_column("ALPHA_J2000", "ra")
    sxt_tab.rename_column("DELTA_J2000", "dec")
    sxt_tab.sort(keys="NUMBER")

    return sxt_tab


def get_filename(access_url: str) -> str:
    """Extracts Level 2 image filename from IRSA access URL.

    Parameters
    ----------
    access_url : `str`
        URL to SPHEREx Level 2 image.

    Returns
    -------
    l2_name : `str`
        Filename of SPHEREx Level 2 image.
    """
    l2_name = access_url.split("/")[-1]
    return l2_name


def get_obs_id(hdu: fits.ImageHDU) -> str:
    """Extracts observation ID from image header.

    Parameters
    ----------
    hdu : `astropy.io.fits.ImageHDU`
        SPHEREx Level 2 image HDU.

    Returns
    -------
    obs_id : `str`
        Observation ID of SPHEREx Level 2 image.
    """
    return hdu.header.get("OBSID")


def get_exp_id(hdu: fits.ImageHDU) -> str:
    """Extracts exposure ID from image header.

    Parameters
    ----------
    hdu : `astropy.io.fits.ImageHDU`
        SPHEREx Level 2 image HDU.

    Returns
    -------
    exp_id : `str`
        Exposure ID of SPHEREx Level 2 image.
    """
    return hdu.header.get("EXPIDN")


def results_summary(table: Table) -> None:
    """Prints a summary of the results table.

    Parameters
    ----------
    table : `astropy.table.Table`
        Table of results.
    """
    grouped_tab = table.group_by('energy_bandpassname')
    print(f"Found {len(table)} images.")
    print("By SPHEREx band:")
    for i, band in enumerate(grouped_tab.groups.keys):
        print(f"\t{band[0]}: {len(grouped_tab.groups[i]):>4} images")


def get_lambda_range(result: Table) -> tuple[Quantity, Quantity]:
    """Extracts the minimum and maximum wavelength from a results table.

    Parameters
    ----------
    result : `astropy.table.Table`
        Table of results.

    Returns
    -------
    lambda_min, lambda_max : `astropy.units.Quantity`
        Minimum and maximum wavelength [microns].
    """
    lambda_min = result["em_min"]*u.m.to(u.um)
    lambda_max = result["em_max"]*u.m.to(u.um)
    return lambda_min, lambda_max
