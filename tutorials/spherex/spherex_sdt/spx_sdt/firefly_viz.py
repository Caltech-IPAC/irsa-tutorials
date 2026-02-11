"""Functions to display images, plots, and tables in Firefly."""

# Related third-party imports
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from firefly_client import FireflyClient

# Local application/library specific imports
from spx_sdt.sdt_utils import get_filename


def preview_query(fc: FireflyClient, access_url: str, coord: SkyCoord, search_radius: Quantity):
    """Preview the results of an SDT query in Firefly.

    Parameters
    ----------
    fc : `firefly_client.FireflyClient`
        An instance of the FireflyClient to use for visualization.
    access_url : `str`
        The IRSA access URL of the FITS image to visualize.
    coord : `astropy.coordinates.SkyCoord`
        The sky coordinates of the target in degrees.
    search_radius : `astropy.units.Quantity`
        The search radius in angular units.
    """
    # Fetch filename
    filename = get_filename(access_url)

    # Show image
    fc.show_fits_image(
        file_input=access_url,
        plot_id=filename,
        Title=filename
    )

    # Add layer for query position
    position = f"J2000; point {coord.ra.degree} {coord.dec.degree} # point=cross color=red width=5"\
               f" text='Position'"
    fc.overlay_region_layer(region_data=position, region_layer_id="layer1",
                            title="Query Position", plot_id=filename)

    # Add a layer for search radii beyond typical SPHEREx pixel scale
    if search_radius > 6.2*u.arcsec:
        # NOTE: Assuming angular units in DS9 format -> " for arcsec,
        #      ' for arcmin, d for deg
        if search_radius.unit == u.arcsec:
            unit_symbol = "\""
        elif search_radius.unit == u.arcmin:
            unit_symbol = "'"
        elif search_radius.unit == u.deg:
            unit_symbol = "d"
        rad_reg = f"J2000; circle {coord.ra.degree} {coord.dec.degree} {search_radius.value}{unit_symbol}"\
                  f" # color=red width=5 text='Search Radius'"
        fc.overlay_region_layer(region_data=rad_reg, region_layer_id="layer2", title="Query Search Radius",
                                plot_id=filename)
