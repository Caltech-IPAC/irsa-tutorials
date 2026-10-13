"""Functions to render dynamic visualizations of SPHEREx images/spectra.
"""

# Standard library imports
from typing import Literal

# Related third party imports
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.units import Quantity
from astropy.wcs import WCS
from bokeh.layouts import column, layout
from bokeh.models import (
    Button, ColorBar, ColumnDataSource, CustomJS, Div, glyphs, HoverTool, Legend, LinearColorMapper,
    NumberFormatter, RangeSlider, ScientificFormatter, TapTool
)
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.palettes import RdBu11, interp_palette
from bokeh.plotting import figure, gridplot, show
import astropy.units as u
import numpy as np

# Local application/library specific imports

__all__ = ["xyz_hover", "copy_selected_button", "plot_preview", "plot_apertures",
    "plot_spectrum", "plot_subtracted_trio", "plot_overlap_trio"]


def xyz_hover(renderers: list[glyphs.Image], z_units: str) -> HoverTool:
    """Creates a hover tool to display the x, y position and z value in an
    image.

    Parameters
    ----------
    renderers : `list` of `bokeh.models.glyphs.Image`
        List of image glyph(s) to read values from.
    z_units : `str`
        Units of z-value in image.

    Returns
    -------
    hover : `bokeh.models.HoverTool`
        Hover tool to add to a bokeh figure.
    """
    hover = HoverTool(
        renderers=renderers,
        tooltips=[
            ("(x, y)", "($x{0.00}, $y{0.00})"),
            (f"Flux ({z_units})", "@image")
        ],
        formatters={"$x": "numeral", "$y": "numeral"},
    )

    return hover


def copy_selected_button(data_source: ColumnDataSource) -> Button:
    """Creates a button to copy selected source IDs from a dataset.

    Parameters
    ----------
    data_source : `bokeh.models.ColumnDataSource`
        Selectable dataset. Must include the column "NUMBER".

    Returns
    -------
    button : `bokeh.models.Button`
        Copies selected source IDs to the clipboard upon click.
    """
    button = Button(label="Copy Selected Source IDs", button_type="primary")
    callback = CustomJS(
        args=dict(s1=data_source),
        code="""
            var inds = s1.selected.indices;
            var d1 = s1.data;
            var idArr = [];
            for (var i=0; i < inds.length; i++) {
                idArr.push(d1['NUMBER'][inds[i]]);
                console.log(d1['NUMBER'][inds[i]]);
            }
            navigator.clipboard.writeText('[' + idArr.join(', ') + ']');
        """,
    )
    button.js_on_event("button_click", callback)

    return button


def plot_preview(filepath: str, coords: SkyCoord, search_radius: Quantity, description: str, clip=False):
    """Plots SPHEREx images with search region displayed.

    Parameters
    ----------
    filepath : `str`
        Path to .fits file.
    coords : `astropy.coordinates.SkyCoord`
        Search position coordinates.
    search_radius : `astropy.units.Quantity`
        Search radius in angular units.
    description : `str`
        Plot title.
    clip : `bool`, optional
        If True, quantile-clips the image for visualizations.
    """
    # Read image HDU data
    with fits.open(filepath) as hdul:
        image_data = hdul[1].data
        wcs = WCS(hdul[1].header)

    image_name = filepath.split("/")[-1]

    # Set up bokeh figure
    plot = figure(
        title=f"{image_name}",
        x_range=(0, np.shape(image_data)[1]),
        y_range=(0, np.shape(image_data)[0]),
        x_axis_label="Pixels",
        y_axis_label="Pixels",
        width=880,
        height=800,
        resizable=True,
        output_backend="webgl"
    )
    plot.title.align = "center"
    plot.title.text_font_size = "16px"
    plot.xaxis.major_tick_line_width = 3
    plot.yaxis.major_tick_line_width = 3

    # Add image
    if clip:
        image_data = np.clip(image_data, np.nanquantile(image_data, .1), np.nanquantile(image_data, .99))
    color = LinearColorMapper(palette="Viridis256", low=np.nanmin(image_data), high=np.nanmax(image_data))
    plot.image(image=[image_data], x=0, y=0, dw=np.shape(image_data)[1],
               dh=np.shape(image_data)[0], color_mapper=color)
    cb = ColorBar(color_mapper=color, location=(0, 0), title="MJy/sr")
    plot.add_layout(cb, "right")

    # Add search region
    x, y = wcs.world_to_pixel(coords)
    radius = (search_radius.to(u.arcsec))/(6.2*u.arcsec)
    region = plot.circle(x=float(x), y=float(y), radius=float(radius), line_color="red", fill_alpha=0,
                         line_width=2.5)
    legend = Legend(
        items=[(description, [region])],
        location="bottom_right",
        border_line_color="black",
        label_text_color="black",
    )
    plot.add_layout(legend, "center")


def plot_apertures(filepath: str, source_tab: Table, label: str, clip=False):
    """Plots SExtractor-derived elliptical apertures over a .fits image.

    Parameters
    ----------
    filepath : `str`
        Path to .fits file.
    source_tab : `astropy.table.Table`
        Table of extracted or selected sources to visualize.
    label : `str`
        Describes sources in source tab (e.g. "Extracted", "Selected").
    clip : `bool`, optional
        Determines whether to quantile-clip the image for visualizations.

    Returns
    -------
    None
    """
    # Read image HDU data
    with fits.open(filepath) as hdul:
        image_data = hdul[1].data

    image_name = filepath.split("/")[-1]

    # Format source table columns for plotting bokeh ellipses
    ellipse_df = source_tab.to_pandas()
    ellipse_df.sort_values(by="NUMBER", inplace=True)

    # SExtractor follows FITS conventions; center of 1st pixel is (1, 1)
    ellipse_df["X_PLOT"] = ellipse_df["X_IMAGE"] - 0.5
    ellipse_df["Y_PLOT"] = ellipse_df["Y_IMAGE"] - 0.5

    # SExtractor output is deg & measures increasing theta CCW from horizontal;
    # Bokeh input is radians & measures increasing theta CW from horizontal
    ellipse_df["THETA_PLOT"] = np.deg2rad(-ellipse_df["THETA_IMAGE"])

    ellipse_df["width"] = ellipse_df["A_IMAGE"] * 2
    ellipse_df["height"] = ellipse_df["B_IMAGE"] * 2

    # Set up data table
    ellipse_source = ColumnDataSource(data=ellipse_df)
    ra_dec_formatter = NumberFormatter(format="0.000000")
    flux_formatter = ScientificFormatter(precision=3, power_limit_high=3)
    columns = [
        TableColumn(field="NUMBER", title="SExtractor ID", formatter=NumberFormatter(font_style="bold"),
                    width=45),
        TableColumn(field="ra", title="RA", formatter=ra_dec_formatter, width=75),
        TableColumn(field="dec", title="DEC", formatter=ra_dec_formatter, width=75),
        TableColumn(field="FLUX_AUTO", title="Flux [uJy]", formatter=flux_formatter, width=70),
    ]
    discovery_table = DataTable(
        source=ellipse_source,
        columns=columns,
        sortable=True,
        selectable="checkbox",
        scroll_to_selection=True,
        multi_selectable=True
    )

    # Copy source IDs from disovery table to clipboard
    button = copy_selected_button(ellipse_source)

    # Set up bokeh figure
    plot = figure(
        title=f"{image_name}",
        x_range=(0, np.shape(image_data)[1]),
        y_range=(0, np.shape(image_data)[0]),
        x_axis_label="Pixels",
        y_axis_label="Pixels",
        width=115,  # wider to accommodate colorbar
        height=100,
        output_backend="webgl"
    )
    plot.title.align = "center"
    plot.title.text_font_size = "16px"
    plot.xaxis.major_tick_line_width = 3
    plot.yaxis.major_tick_line_width = 3

    # Add image
    if clip:
        image_data = np.clip(image_data, np.nanquantile(image_data, .1), np.nanquantile(image_data, .99))
    color = LinearColorMapper(palette="Viridis256", low=np.nanmin(image_data), high=np.nanmax(image_data))
    image = plot.image(image=[image_data], x=0, y=0, dw=np.shape(image_data)[1],
                       dh=np.shape(image_data)[0], color_mapper=color)
    cb = ColorBar(color_mapper=color, location=(0, 0), title="uJy", width=10)
    plot.add_layout(cb, "right")

    # Add SExtractor ellipses
    ellipse = plot.ellipse("X_PLOT", "Y_PLOT", width="width", height="height", angle="THETA_IMAGE",
                           source=ellipse_source,

                           # initial glyph properties
                           line_color="red", line_width=2.5, fill_alpha=0,

                           # selected glyph properties
                           selection_line_color="fuchsia", selection_line_width=2.5,
                           selection_fill_color="fuchsia", selection_fill_alpha=1.0,

                           # non-selected glyph properties
                           nonselection_line_color="red", nonselection_line_width=2.5,
                           nonselection_fill_alpha=0, nonselection_line_alpha=1.0,

                           # muted glyph properties
                           muted_fill_alpha=0, muted_line_color="red", muted_line_alpha=0.25)
    legend = Legend(
        items=[(label, [ellipse])],
        location="bottom_right",
        click_policy="mute",
        border_line_color="black",
        label_text_color="black",
    )
    plot.add_layout(legend, "center")

    # View source data on mouse hover
    image_hover = HoverTool(
        renderers=[image],
        tooltips=[
            ("(X, Y)", "($x{0.00}, $y{0.00})"),
        ],
        formatters={"$x": "numeral", "$y": "numeral"},
        visible=False,
    )
    ellipse_hover = HoverTool(
        renderers=[ellipse],
        tooltips=[
            ("SExtractor ID", "@NUMBER"),
            ("(X, Y)", "($x{0.00}, $y{0.00})"),
            ("(RA, DEC)", "(@ra{0.000000}, @dec{0.000000})"),
            ("Flux (uJy)", "@FLUX_AUTO"),
        ],
        formatters={"$x": "numeral", "$y": "numeral", "@ra": "numeral", "@dec": "numeral"},
        visible=False,
    )
    ellipse_tap = TapTool(
        renderers=[ellipse],
        mode="append"
    )
    plot.add_tools(image_hover, ellipse_hover, ellipse_tap)

    # Add instructions for interaction
    div = Div(
        text="""
            <h2>Source Discovery</h2>
            <p>
                To <b>inspect</b> sources, hover your mouse over the image or pan and zoom by activating the
                toolbar icons.</br>
                To <b>select</b> sources, click an ellipse in the image or click a checkbox in the table.</br>
                To copy the selected source IDs to your clipboard, click the button below.
            </p>
        """
    )

    # Assemble elements together & display them
    widget = layout(
        [[plot, column(div, button, discovery_table, margin=(5, 5, 5, 20))]],
        sizing_mode="scale_both",
        aspect_ratio=2,
    )
    show(widget)


def plot_spectrum(wavelength: np.ndarray, bandwidth: np.ndarray, flux: np.ndarray, variance: np.ndarray,
                  ab_mag: np.ndarray, ab_mag_err: np.ndarray, n_bad: np.ndarray,
                  type: Literal["flux", "magnitude"], source_id: str, ap_size: float):
    """Plots interactive source spectrum.

    Parameters
    ----------
    wavelength : `numpy.ndarray`
        Wavelengths [microns].
    bandwidth : `numpy.ndarray`
        Wavelength bin bandwidths [microns].
    flux : `numpy.ndarray`
        Fluxes [mJy].
    variance : `numpy.ndarray`
        Flux variances [mJy^2].
    ab_mag : `numpy.ndarray`
        AB magnitudes.
    ab_mag_err : `numpy.ndarray`
        AB magnitude errors.
    n_bad : `numpy.ndarray`
        Number of bad pixels per measurement.
    type : `str`
        Type of data to plot: "flux" or "magnitude".
    source_id : `str`
        Source ID.
    ap_size : `float`
        Aperture size [pixels].

    Notes
    -------
    - For mathematical/symbolic labels, we use unicode instead of LaTeX. bokeh
      supports LaTeX, but it expects MathJax to be in the notebook environment.
      This is true in JupyterLab, but not in IDE notebook renderers.
    - The bokeh API does not support error bars, so we use line segments with
      the same legend labels as the scatter points to fake this.
    """
    # Only plot points with finite flux/AB mag and minimal bad pixels
    mask = ((n_bad < 1)  # fewer than this many "bad" pixels
            & np.isfinite(ab_mag)
            & np.isfinite(flux))
    wavelength = wavelength[mask]
    bandwidth = bandwidth[mask]
    flux = flux[mask]
    variance = variance[mask]
    ab_mag = ab_mag[mask]
    ab_mag_err = ab_mag_err[mask]
    n_bad = n_bad[mask]

    # Determine plot type
    if type == "magnitude":
        y_label = "AB magnitude"
        y_value = ab_mag
        y_err = ab_mag_err
    elif type == "flux":
        y_label = "Flux (mJy)"
        y_value = flux
        y_err = np.sqrt(np.array(variance))

    # Set up bokeh figure
    plot = figure(
        x_range=(0.60, 5.15),
        title=f"{source_id} Spectrum",
        x_axis_label="\u03BB (\u03BCm)",
        y_axis_label=y_label,
        width=640,
        height=360,
        output_backend="webgl"
    )
    plot.title.align = "center"
    plot.title.text_font_size = "16px"
    plot.xaxis.major_tick_line_width = 2
    plot.yaxis.major_tick_line_width = 2
    plot.y_range.flipped = True if type == "magnitude" else False

    # Create column data source to share data across glyphs
    spec_data = ColumnDataSource(data=dict(
        x=wavelength,
        bdw=bandwidth,
        y=y_value,
        y_low=y_value - y_err,
        y_high=y_value + y_err,
        y_pm=[f"{y:.2f} \u00B1 {e:.2f}" for y, e in zip(y_value, y_err)],
    ))

    # Add measurement & measurement error
    plot.scatter(x='x', y='y', source=spec_data, color="darkgreen",
                 legend_label=f"r = {ap_size} pix")
    plot.segment(y0='y_low', y1='y_high', x0='x', x1='x', source=spec_data, color="darkgreen", line_width=1.5,
                 legend_label=f"r = {ap_size} pix")  # y error; shared legend label to mute glyphs together

    # Add hover tooltip
    hover = HoverTool(
        tooltips=[
            ("\u03BB (\u03BCm)", "@x{0.000000}"),
            ("Bandwidth (\u03BCm)", "@bdw{0.000000}"), 
            (y_label, "@y_pm"),
        ],
        formatters={"@x": "numeral", "@bdw": "numeral"}
    )
    plot.add_tools(hover)

    # Style legend
    plot.legend.location = "top_right"
    plot.legend.click_policy = "mute"

    # Set up wavelength range slider
    div = Div(
        text="""
            <h2>Wavelength Range</h2>
            <p>To adjust the wavelengths plotted, click the slider below.</p>
        """
    )
    range_slider = RangeSlider(
        title="\u03BB (\u03BCm)",
        start=0.60,
        end=5.15,
        step=0.1,
        value=(plot.x_range.start, plot.x_range.end),
        align="start",
        aspect_ratio=10,
        sizing_mode="scale_width",
    )
    range_slider.js_link("value", plot.x_range, "start", attr_selector=0)
    range_slider.js_link("value", plot.x_range, "end", attr_selector=1)
    range_slider.align = "center"

    # Assemble elements together & display them
    bokeh_spectrum = layout(
        children=[
            [plot],
            [column(div, range_slider, margin=(5, 5, 5, 5))],
        ],
        sizing_mode="scale_width",
    )
    show(bokeh_spectrum)


def plot_subtracted_trio(lambda_bin1_info: dict, lambda_bin2_info: dict, subtracted_info: dict,
                         cbar_tune_across: bool = False,
                         clip: bool = False):
    """Plots subtracted color cutout.

    Parameters
    ----------
    lambda_bin1_info, lambda_bin2_info, subtracted_info : `dict`
        Keys "img_data", "title", "units" should map to objects
        `numpy.ndarray`, `str`, `str`.
    cbar_tune_across : `bool`, optional
        Whether to tune the vmin/vmax of the colorbars across bins
        or by individual image. The default is `False`.
    clip : `bool`, optional
        Whether to quantile-clip the image for visualizations. The
        default is `False`.

    Returns
    -------
    None.
    """
    # Clip data for visualizations if necessary
    if clip:
        data1 = np.clip(lambda_bin1_info["img_data"],
                        np.nanquantile(lambda_bin1_info["img_data"], .1),
                        np.nanquantile(lambda_bin1_info["img_data"], .99))
        data2 = np.clip(lambda_bin2_info["img_data"],
                        np.nanquantile(lambda_bin2_info["img_data"], .1),
                        np.nanquantile(lambda_bin2_info["img_data"], .99))
        subtracted = np.clip(subtracted_info["img_data"],
                             np.nanquantile(subtracted_info["img_data"], .1),
                             np.nanquantile(subtracted_info["img_data"], .99))
    else:
        data1 = lambda_bin1_info["img_data"]
        data2 = lambda_bin2_info["img_data"]
        subtracted = subtracted_info["img_data"]

    # Tune colorbars across wavelength bins
    if cbar_tune_across:
        vmin = np.nanmin(np.stack([data1, data2]))
        vmax = np.nanmax(np.stack([data1, data2]))

    # Plot image in 1st wavelength bin
    bin1_fig = figure(
        title=lambda_bin1_info['title'],
        x_range=(0, np.shape(data1)[1]),
        y_range=(0, np.shape(data1)[0]),
        x_axis_label="Pixels",
        y_axis_label="Pixels",
        width=115,  # wider to accommodate colorbar
        height=100,
        output_backend="webgl"
    )
    bin1_fig.title.align = "center"

    color = LinearColorMapper(palette="Viridis256",
                              low=vmin if cbar_tune_across else np.nanmin(data1),
                              high=vmax if cbar_tune_across else np.nanmax(data1))
    img1 = bin1_fig.image(image=[data1], x=0, y=0, dw=np.shape(data1)[1], dh=np.shape(data1)[0],
                          color_mapper=color)
    cb = ColorBar(color_mapper=color, location=(0, 0), title=lambda_bin1_info["units"], width=10)
    bin1_fig.add_layout(cb, "right")
    bin1_fig.add_tools(xyz_hover(renderers=[img1], z_units=lambda_bin1_info["units"]))

    # Plot image in 2nd wavelength bin
    bin2_fig = figure(
        title=lambda_bin2_info['title'],
        x_range=bin1_fig.x_range,  # link panning/zoom to other plots
        y_range=bin1_fig.y_range,  # link panning/zoom to other plots
        x_axis_label="Pixels",
        y_axis_label="Pixels",
        width=115,  # wider to accommodate colorbar
        height=100,
        output_backend="webgl"
    )
    bin2_fig.title.align = "center"

    color = LinearColorMapper(palette="Viridis256",
                              low=vmin if cbar_tune_across else np.nanmin(data2),
                              high=vmax if cbar_tune_across else np.nanmax(data2))
    img2 = bin2_fig.image(image=[data2], x=0, y=0, dw=np.shape(data2)[1], dh=np.shape(data2)[0],
                          color_mapper=color)
    cb = ColorBar(color_mapper=color, location=(0, 0), title=lambda_bin2_info["units"], width=10)
    bin2_fig.add_layout(cb, "right")
    bin2_fig.add_tools(xyz_hover(renderers=[img2], z_units=lambda_bin2_info["units"]))

    # Plot subtracted color overview
    subtracted_fig = figure(
        title=subtracted_info["title"],
        x_range=bin1_fig.x_range,  # link panning/zoom to other plots
        y_range=bin1_fig.y_range,  # link panning/zoom to other plots
        x_axis_label="Pixels",
        y_axis_label="Pixels",
        width=115,  # wider to accommodate colorbar
        height=100,
        output_backend="webgl"
    )
    subtracted_fig.title.align = "center"

    bokeh_rdbu = interp_palette(RdBu11, 256)
    color = LinearColorMapper(
        palette=bokeh_rdbu[::-1], low=-np.nanmax(subtracted), high=np.nanmax(subtracted))
    img_subtracted = subtracted_fig.image(image=[subtracted], x=0, y=0, dw=np.shape(subtracted)[1],
                                          dh=np.shape(subtracted)[0], color_mapper=color)
    cb = ColorBar(color_mapper=color, location=(0, 0), title=subtracted_info["units"], width=10)
    subtracted_fig.add_layout(cb, "right")
    subtracted_fig.add_tools(xyz_hover(renderers=[img_subtracted], z_units=subtracted_info["units"]))

    # Assemble elements together & display them
    linked_plots = gridplot(
        [
            [bin1_fig, bin2_fig, subtracted_fig]
        ],
        sizing_mode="scale_both",
    )
    show(linked_plots)


def plot_overlap_trio(plot1_info: dict, plot2_info: dict, plot3_info: dict,
                      clip: bool = False):
    """Plots overlapping SPHEREx images from different observations.

    Parameters
    ----------
    plot1_info, plot2_info, plot3_info : `dict`
        Keys "img_data", "title", "units" should map to objects
        `numpy.ndarray`, `str`, `str`.
    clip : `bool`, optional
        Whether to quantile-clip the image for visualizations. The
        default is `False`.
    """
    # If requested, clip data for visualizations
    if clip:
        data1 = np.clip(plot1_info["img_data"],
                        np.nanquantile(plot1_info["img_data"], .1),
                        np.nanquantile(plot1_info["img_data"], .99))
        data2 = np.clip(plot2_info["img_data"],
                        np.nanquantile(plot2_info["img_data"], .1),
                        np.nanquantile(plot2_info["img_data"], .99))
        data3 = np.clip(plot3_info["img_data"],
                        np.nanquantile(plot3_info["img_data"], .1),
                        np.nanquantile(plot3_info["img_data"], .99))
    else:
        data1 = plot1_info["img_data"]
        data2 = plot2_info["img_data"]
        data3 = plot3_info["img_data"]

    # 1st plot
    plot1_fig = figure(
        title=plot1_info['title'],
        x_range=(0, np.shape(data1)[1]),
        y_range=(0, np.shape(data1)[0]),
        x_axis_label="Pixels",
        y_axis_label="Pixels",
        width=115,  # wider to accommodate colorbar
        height=100,
        output_backend="webgl"
    )
    plot1_fig.title.align = "center"

    color = LinearColorMapper(palette="Viridis256", low=np.nanmin(data1), high=np.nanmax(data1))
    img1 = plot1_fig.image(image=[data1], x=0, y=0, dw=np.shape(data1)[1], dh=np.shape(data1)[0],
                           color_mapper=color)
    cb = ColorBar(color_mapper=color, location=(0, 0), title=plot1_info["units"], width=10)
    plot1_fig.add_layout(cb, "right")
    plot1_fig.add_tools(xyz_hover(renderers=[img1], z_units=plot1_info["units"]))

    # 2nd plot
    plot2_fig = figure(
        title=plot2_info['title'],
        x_range=plot1_fig.x_range,  # link panning/zoom to other plots
        y_range=plot1_fig.y_range,  # link panning/zoom to other plots
        x_axis_label="Pixels",
        y_axis_label="Pixels",
        width=115,  # wider to accommodate colorbar
        height=100,
        output_backend="webgl"
    )
    plot2_fig.title.align = "center"
    color = LinearColorMapper(palette="Viridis256", low=np.nanmin(data2), high=np.nanmax(data2))
    img2 = plot2_fig.image(image=[data2], x=0, y=0, dw=np.shape(data2)[1], dh=np.shape(data2)[0],
                           color_mapper=color)
    cb = ColorBar(color_mapper=color, location=(0, 0), title=plot2_info["units"], width=10)
    plot2_fig.add_layout(cb, "right")
    plot2_fig.add_tools(xyz_hover(renderers=[img2], z_units=plot2_info["units"]))

    # 3rd plot
    plot3_fig = figure(
        title=plot3_info['title'],
        x_range=plot1_fig.x_range,  # link panning/zoom to other plots
        y_range=plot1_fig.y_range,  # link panning/zoom to other plots
        x_axis_label="Pixels",
        y_axis_label="Pixels",
        width=115,  # wider to accommodate colorbar
        height=100,
        output_backend="webgl"
    )
    plot3_fig.title.align = "center"
    color = LinearColorMapper(palette="Viridis256", low=np.nanmin(data3), high=np.nanmax(data3))
    img3 = plot3_fig.image(image=[data3], x=0, y=0, dw=np.shape(data3)[1], dh=np.shape(data3)[0],
                           color_mapper=color)
    cb = ColorBar(color_mapper=color, location=(0, 0), title=plot3_info["units"], width=10)
    plot3_fig.add_layout(cb, "right")
    plot3_fig.add_tools(xyz_hover(renderers=[img3], z_units=plot3_info["units"]))

    # Assemble elements together & display them
    linked_plots = gridplot(
        [
            [plot1_fig, plot2_fig, plot3_fig]
        ],
        sizing_mode="scale_both",
    )
    show(linked_plots)
