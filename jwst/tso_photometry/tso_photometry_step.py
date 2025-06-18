import logging

from stdatamodels.jwst.datamodels import CubeModel, GainModel

from jwst.lib import reffile_utils
from jwst.lib.catalog_utils import replace_suffix_ext
from jwst.stpipe import Step
from jwst.tso_photometry.tso_photometry import tso_aperture_photometry

__all__ = ["TSOPhotometryStep"]

log = logging.getLogger(__name__)


class TSOPhotometryStep(Step):
    """Perform circular aperture photometry on imaging Time Series Observations (TSO)."""

    class_alias = "tso_photometry"

    spec = """
        radius = float(default=3.0) # Aperture radius in pixels
        radius_inner = float(default=4.0) # Background annulus inner radius in pixels
        radius_outer = float(default=5.0) # Background annulus outer radius in pixels
        save_catalog = boolean(default=False)  # save exposure-level catalog
    """  # noqa: E501

    reference_file_types = ["tsophot", "gain"]

    def process(self, input_data):
        """
        Do the tso photometry processing.

        Parameters
        ----------
        input_data : str or CubeModel
            Filename for a FITS image, or a `CubeModel`.

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
            log.debug(f"xcenter = {xcenter}")
            log.debug(f"ycenter = {ycenter}")

            # Compute the aperture photometry
            catalog = tso_aperture_photometry(
                model, xcenter, ycenter, self.radius, self.radius_inner, self.radius_outer, gain_2d
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
