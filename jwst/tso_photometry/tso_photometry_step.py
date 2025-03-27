#!/usr/bin/env python
from stdatamodels.jwst.datamodels import CubeModel, GainModel
import numpy as np

from ..stpipe import Step
from ..lib.catalog_utils import replace_suffix_ext
from .tso_photometry import tso_aperture_photometry

__all__ = ["TSOPhotometryStep"]


class TSOPhotometryStep(Step):
    """Perform circular aperture photometry on imaging Time Series
    Observations (TSO).

    Parameters
    -----------
    input : str or `CubeModel`
        Filename for a FITS image or association table, or a `CubeModel`.

    centroid_results : numpy.ndarray
        The centroid results from the `TSOCentroidingStep`.

    psfWidth_results : numpy.ndarray
        The PSF width results from the `TSOCentroidingStep`.

    psfFlux_results : numpy.ndarray
        The PSF flux results from the `TSOCentroidingStep`.

    Returns
    -------
    catalog : `~astropy.table.QTable`
        Astropy QTable (Quantity Table) containing the source photometry.
    """

    class_alias = "tso_photometry"

    spec = """
        radius = float(default=5.0)  # radius of circular source aperture
        radius_inner = float(default=16.0)  # inner radius of background annulus
        radius_outer = float(default=36.0)  # outer radius of background annulus
        moving = boolean(default=False)  # aperture moves with fitted centroids?
        save_catalog = boolean(default=False)  # save exposure-level catalog
    """  # noqa: E501

    reference_file_types = ['gain']

    def process(self, input_data, centroid_results, psfWidth_results, psfFlux_results):
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

            # Get the tsophot reference file
            tsophot_filename = self.get_reference_file(model, "tsophot")
            self.log.debug(f"Reference file name = {tsophot_filename}")
            if tsophot_filename == "N/A":
                self.log.warning("No TSOPHOT reference file found;")
                self.log.warning("the tso_photometry step will be skipped.")
                return None

            # Get the gain reference file
            gain_filename = self.get_reference_file(model, 'gain')
            gain_model = GainModel(gain_filename)

            # Get the x/y center from the centroid results
            # Meta values are 1-based, convert to 0-based
            start = model.meta.exposure.integration_start - 1
            end = model.meta.exposure.integration_end
            xcenter = centroid_results[start:end, 0]
            ycenter = centroid_results[start:end, 1]
            xWidth = psfWidth_results[start:end, 0]
            yWidth = psfWidth_results[start:end, 1]
            psfFlux = psfFlux_results[start:end]

            if self.moving:
                # Move the aperture along with the moving fitted centroids
                # from the relevant integrations
                x = xcenter
                y = ycenter
                # FINDME: The above code likely won't work if there are multiple exposures
            else:
                # Use the median of the fitted centroids from all integrations
                x = np.ones(end-start)*np.median(centroid_results[:, 0])
                y = np.ones(end-start)*np.median(centroid_results[:, 1])

            # Compute the aperture photometry
            catalog = tso_aperture_photometry(model, x, y, self.radius,
                                              self.radius_inner, self.radius_outer,
                                              gain_model, xcenter, ycenter,
                                              xWidth, yWidth, psfFlux, self.log)

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
                self.log.info(f"Wrote TSO photometry catalog: {cat_filepath}")

        return catalog
