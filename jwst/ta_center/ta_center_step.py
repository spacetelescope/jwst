import logging

import jwst.datamodels as dm
from jwst.assign_wcs.miri import retrieve_filter_offset
from jwst.stpipe import Step
from jwst.ta_center.ta_center import NoFinitePixelsError, center_from_ta_image

__all__ = ["TACenterStep"]

log = logging.getLogger(__name__)


class TACenterStep(Step):
    """Determine position of target source from TA verification image."""

    class_alias = "ta_center"

    spec = """
    ta_file = string(default=None)  # Target acquisition image file name
    skip = boolean(default=True)  # Skip this step by default
    """  # noqa: E501

    reference_file_types = ["specwcs", "pathloss", "filteroffset"]

    def process(self, step_input):
        """
        Process target acquisition data.

        Parameters
        ----------
        step_input : str, `~jwst.datamodels.container.ModelContainer`, \
                `~stdatamodels.jwst.datamodels.ImageModel`, \
                `~stdatamodels.jwst.datamodels.CubeModel`
            The input data model or association.

        Returns
        -------
        result : `~jwst.datamodels.container.ModelContainer`, \
                `~stdatamodels.jwst.datamodels.ImageModel`, \
                `~stdatamodels.jwst.datamodels.CubeModel`
            The output data model or association with TA centering applied.
        """
        result = self.prepare_output(step_input)
        if isinstance(result, dm.ModelContainer):
            # Extract science and TA image from container
            result, ta_model = self._ta_image_from_container(result)
            self.ta_file = ta_model

        # Ensure TA file is provided
        if str(self.ta_file).lower() == "none":
            log.error("No target acquisition file provided. Step will be SKIPPED.")
            result.meta.cal_step.ta_center = "SKIPPED"
            return result

        # Check that this is a point source
        if result.meta.target.source_type != "POINT":
            log.error("TA centering is only implemented for point sources. Step will be SKIPPED.")
            result.meta.cal_step.ta_center = "SKIPPED"
            return result

        # Check exposure type
        exp_type = result.meta.exposure.type
        if exp_type not in ["MIR_LRS-FIXEDSLIT", "MIR_LRS-SLITLESS"]:
            log.error(
                "TA centering is only implemented for MIR_LRS-FIXEDSLIT and"
                " MIR_LRS-SLITLESS modes. Step will be SKIPPED."
            )
            result.meta.cal_step.ta_center = "SKIPPED"
            return result

        # Read the TA image
        with dm.open(self.ta_file) as ta_model:
            log.info(f"Performing TA centering using file: {self.ta_file}")
            ta_image = ta_model.data
            wavelength = _get_wavelength(ta_model.meta.instrument.filter)
            if wavelength is None:
                log.error(
                    "Unknown filter; cannot determine wavelength for TA centering."
                    " Step will be SKIPPED."
                )
                result.meta.cal_step.ta_center = "SKIPPED"
                return result

        # read specwcs to get necessary reference points on detector
        reffile = self.get_reference_file(result, "specwcs")
        refmodel = dm.MiriLRSSpecwcsModel(reffile)

        # Get subarray information from TA image
        subarray_name = ta_model.meta.subarray.name
        xstart = ta_model.meta.subarray.xstart  # 1-indexed in FITS
        ystart = ta_model.meta.subarray.ystart  # 1-indexed in FITS
        log.info(f"TA subarray: {subarray_name}, origin: ({xstart}, {ystart})")

        is_slit = exp_type == "MIR_LRS-FIXEDSLIT"
        if is_slit:
            log.info("Fitting Airy disk with slit mask model for LRS FIXEDSLIT mode")
            ref_center = (refmodel.meta.x_ref, refmodel.meta.y_ref)
            pathloss_file = self.get_reference_file(result, "pathloss")
        else:
            log.info("Fitting Airy disk for slitless mode")
            ref_center = (refmodel.meta.x_ref_slitless, refmodel.meta.y_ref_slitless)
            pathloss_file = None

        try:
            x_center, y_center = center_from_ta_image(
                ta_image,
                wavelength,
                ref_center,
                subarray_origin=(xstart, ystart),
                pathloss_file=pathloss_file,
            )
        except NoFinitePixelsError as e:
            log.error(f"Error during TA centering: {e}. Step will be SKIPPED.")
            result.meta.cal_step.ta_center = "SKIPPED"
            return result

        # Apply filter offsets
        filteroffset_file = self.get_reference_file(ta_model, "filteroffset")
        with dm.FilteroffsetModel(filteroffset_file) as filteroffset:
            col_offset, row_offset = retrieve_filter_offset(
                filteroffset, ta_model.meta.instrument.filter
            )
            log.info(f"Applying filter offsets: column={col_offset}, row={row_offset}")
            x_center += col_offset
            y_center += row_offset

        # Set completion status
        result.source_xpos = x_center
        result.source_ypos = y_center
        result.meta.cal_step.ta_center = "COMPLETE"

        return result

    def _ta_image_from_container(self, container):
        """
        Extract the TA image from a container (association or ModelContainer).

        Parameters
        ----------
        container : ModelContainer
            The input container.

        Returns
        -------
        sci_model : DataModel
            The science data model.
        ta_model : DataModel
            The TA image data model.
        """
        sci_idx = container.ind_asn_type("science")
        sci_model = container[sci_idx[0]]
        ta_idx = container.ind_asn_type("target_acquisition")
        if not len(ta_idx):
            ta_model = None
        else:
            ta_model = container[ta_idx[0]]
        return sci_model, ta_model


def _get_wavelength(filter_name):  # numpydoc ignore=RT01
    """Map filter name to central wavelength in microns."""
    filter_wavelengths = {
        "F560W": 5.6,
        "F770W": 7.7,
        "F1000W": 10.0,
        "F1130W": 11.3,
        "F1280W": 12.8,
        "F1500W": 15.0,
        "F1800W": 18.0,
        "F2100W": 21.0,
        "F2550W": 25.5,
    }
    return filter_wavelengths.get(filter_name, None)
