import logging

import numpy as np
from stdatamodels.jwst.datamodels import CubeModel, GainModel

from jwst.lib import reffile_utils
from jwst.lib.catalog_utils import replace_suffix_ext
from jwst.stpipe import Step
from jwst.tso_photometry.tso_photometry import (
    convert_data_units,
    tso_aperture_photometry,
    tso_source_centroid,
)

__all__ = ["TSOPhotometryStep"]

log = logging.getLogger(__name__)


class TSOPhotometryStep(Step):
    """Perform circular aperture photometry on imaging Time Series Observations (TSO)."""

    class_alias = "tso_photometry"

    spec = """
        save_catalog = boolean(default=False)  # Save exposure-level catalog
        radius = float(default=3.0) # Aperture radius in pixels
        radius_inner = float(default=4.0) # Background annulus inner radius in pixels
        radius_outer = float(default=5.0) # Background annulus outer radius in pixels
        centroid_source = boolean(default=True)  # Centroid source before photometry
        search_box_width = integer(default=41)  # Box width for initial source search; must be odd.
        fit_box_width = integer(default=11)  # Box width for centroid fit; must be odd.
        moving_centroid = boolean(default=False)  # Fit centroid values for each integration
    """  # noqa: E501

    reference_file_types = ["gain"]

    def process(self, input_data):
        """
        Do the TSO photometry processing.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.CubeModel`
            Filename for a FITS image, or a `~stdatamodels.jwst.datamodels.CubeModel`.

        Returns
        -------
        catalog : `~astropy.table.QTable`
            Astropy QTable (Quantity Table) containing the source photometry.
        """
        # Open the input as a CubeModel
        with CubeModel(input_data) as model:
            # Need the FITS WCS X/YREF_SCI values for setting the
            # photometry aperture location
            if model.meta.wcsinfo.siaf_xref_sci is None:
                raise ValueError("XREF_SCI is missing.")
            if model.meta.wcsinfo.siaf_yref_sci is None:
                raise ValueError("YREF_SCI is missing.")
            if model.meta.bunit_data is None:
                raise ValueError("BUNIT for data array is missing.")
            if model.meta.bunit_err is None:
                raise ValueError("BUNIT for error array is missing.")

            xcenter = model.meta.wcsinfo.siaf_xref_sci - 1  # 1-based origin
            ycenter = model.meta.wcsinfo.siaf_yref_sci - 1  # 1-based origin

            # Get the gain reference file
            gain_filename = self.get_reference_file(model, "gain")
            gain_m = GainModel(gain_filename)
            # Get the relevant 2D gain values from the model
            if reffile_utils.ref_matches_sci(model, gain_m):
                gain_2d = gain_m.data
            else:
                log.info("Extracting gain subarray to match science data")
                gain_2d = reffile_utils.get_subarray_model(model, gain_m).data

            log.debug(f"radius = {self.radius}")
            log.debug(f"radius_inner = {self.radius_inner}")
            log.debug(f"radius_outer = {self.radius_outer}")
            log.debug(f"initial xcenter = {xcenter}")
            log.debug(f"initial ycenter = {ycenter}")

            # Convert the data units as needed
            convert_data_units(model, gain_2d)

            # Centroid the source if desired
            if self.centroid_source:
                # Check the box sizes: they must be odd integers
                boxes = ["search_box_width", "fit_box_width"]
                for box in boxes:
                    box_input = getattr(self, box)
                    if box_input % 2 == 0:
                        setattr(self, box, box_input - 1)
                        log.warning(
                            "Even box widths are not supported."
                            f" Rounding the {box} down to {box_input - 1}."
                        )

                centroid_x, centroid_y, psf_width_x, psf_width_y, psf_flux = tso_source_centroid(
                    model,
                    xcenter,
                    ycenter,
                    search_box_width=self.search_box_width,
                    fit_box_width=self.fit_box_width,
                    source_radius=self.radius_inner,
                )

                if np.all(np.isnan(centroid_x)) or np.all(np.isnan(centroid_y)):
                    log.warning("Centroid fit failed. Using initial estimate for source location.")
                    xc = xcenter
                    yc = ycenter
                    log.info(f"Using planned center x,y = {xc:.2f},{yc:.2f}")
                elif not self.moving_centroid:
                    xc = np.nanmedian(centroid_x)
                    yc = np.nanmedian(centroid_y)
                    log.info(f"Using median centroid x,y = {xc:.2f},{yc:.2f}")
                else:
                    xc = centroid_x
                    yc = centroid_y
                    log.info(
                        "Using moving centroid. "
                        f"Median x,y = {np.nanmedian(xc):.2f},{np.nanmedian(yc):.2f}"
                    )

            else:
                xc = xcenter
                yc = ycenter
                centroid_x = None
                centroid_y = None
                psf_width_x = None
                psf_width_y = None
                psf_flux = None
                log.info(f"Using planned center x,y = {xc:.2f},{yc:.2f}")

            # Compute the aperture photometry
            catalog = tso_aperture_photometry(
                model,
                xc,
                yc,
                self.radius,
                self.radius_inner,
                self.radius_outer,
                centroid_x=centroid_x,
                centroid_y=centroid_y,
                psf_width_x=psf_width_x,
                psf_width_y=psf_width_y,
                psf_flux=psf_flux,
            )

            # Save the photometry in an output catalog
            if self.save_catalog:
                old_suffixes = ["calints", "crfints"]
                output_dir = self.search_attr("output_dir")
                cat_filepath = replace_suffix_ext(
                    model.meta.filename,
                    old_suffixes,
                    "phot",
                    output_ext="ecsv",
                    output_dir=output_dir,
                )
                catalog.write(cat_filepath, format="ascii.ecsv", overwrite=True)
                log.info(f"Wrote TSO photometry catalog: {cat_filepath}")

        return catalog
