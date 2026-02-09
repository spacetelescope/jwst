import logging

import jwst.datamodels as dm
from jwst.assign_wcs.miri import retrieve_filter_offset
from jwst.stpipe import Step
from jwst.targ_centroid.targ_centroid import (
    BadFitError,
    NoFinitePixelsError,
    WCSError,
    center_from_ta_image,
    find_dither_position,
)

__all__ = ["TargCentroidStep"]

log = logging.getLogger(__name__)


class TargCentroidStep(Step):
    """Determine position of target source from TA verification image."""

    class_alias = "targ_centroid"

    spec = """
    ta_file = string(default=None)  # Target acquisition image file name
    skip = boolean(default=True)  # Skip this step by default
    """  # noqa: E501

    reference_file_types = ["filteroffset"]

    def process(self, step_input):
        """
        Determine target position from TA verification image.

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
            container = result
            # Extract science and TA image from container
            result, ta_model = self._ta_image_from_container(result)
            self.ta_file = ta_model
        else:
            container = None

        # Ensure TA file is provided
        if str(self.ta_file).lower() == "none":
            log.error("No target acquisition file provided. Step will be SKIPPED.")
            result.meta.cal_step.targ_centroid = "SKIPPED"
            result = self._rebuild_container(container, result)
            return result

        # Check that this is a point source
        if result.meta.target.source_type != "POINT":
            log.warning(
                "TargCentroidStep is only intended for point sources. "
                f"Input data has source_type={result.meta.target.source_type}. "
                "Step may not perform as expected."
            )

        # Check exposure type
        exp_type = result.meta.exposure.type
        if exp_type not in ["MIR_LRS-FIXEDSLIT", "MIR_LRS-SLITLESS"]:
            log.error(
                "TA centering is only implemented for MIR_LRS-FIXEDSLIT and"
                " MIR_LRS-SLITLESS modes. Step will be SKIPPED."
            )
            result.meta.cal_step.targ_centroid = "SKIPPED"
            result = self._rebuild_container(container, result)
            return result

        log.info(f"Performing TA centering using file: {self.ta_file}")
        with dm.open(self.ta_file) as ta_model:
            # check TA verification image and exposure type
            ta_exptype = ta_model.meta.exposure.type
            if (not isinstance(ta_model, dm.ImageModel)) or (ta_exptype != "MIR_TACONFIRM"):
                log.error(
                    "TA verification image must be an ImageModel with exposure type MIR_TACONFIRM. "
                    f"Input TA file is of type {type(ta_model)} with exposure type {ta_exptype}. "
                    "Step will be SKIPPED."
                )
                result.meta.cal_step.targ_centroid = "SKIPPED"
                result = self._rebuild_container(container, result)
                return result

            # Find expected source position by assigning WCS to TA model
            # and translating dither offsets to detector coordinates
            try:
                ref_center, ta_dither_offset = find_dither_position(ta_model)
            except WCSError as e:
                log.error(f"Error when assigning WCS to {self.ta_file}: {e}. Step will be SKIPPED.")
                result.meta.cal_step.targ_centroid = "SKIPPED"
                result = self._rebuild_container(container, result)
                return result
            log.debug(f"Reference center in detector subarray (0-indexed): {ref_center}")

            # Put attributes needed later into local variables so we can exit context manager
            ta_image = ta_model.data
            ta_filter = ta_model.meta.instrument.filter
            filteroffset_file = self.get_reference_file(ta_model, "filteroffset")
            subarray_name = ta_model.meta.subarray.name
            xstart = ta_model.meta.subarray.xstart  # 1-indexed in FITS
            ystart = ta_model.meta.subarray.ystart  # 1-indexed in FITS
            log.debug(f"TA subarray: {subarray_name}, origin: ({xstart}, {ystart})")

        # Compute centroid from TA image
        try:
            (x_center, y_center), (x_center_sub, y_center_sub) = center_from_ta_image(
                ta_image,
                ref_center,
                subarray_origin=(xstart, ystart),
            )
        except (NoFinitePixelsError, BadFitError) as e:
            log.error(f"Error during TA centering: {e}. Step will be SKIPPED.")
            result.meta.cal_step.targ_centroid = "SKIPPED"
            result = self._rebuild_container(container, result)
            return result

        # Store TA centering results in output model
        result.ta_xpos = x_center_sub
        result.ta_ypos = y_center_sub

        log.info(
            "Fitted source position on TA verification image: "
            f"({result.ta_xpos:.2f}, {result.ta_ypos:.2f})"
        )

        # undo dither offset of TA image.
        # At this stage we should be near the reference point from the specwcs,
        # (326.13, 300.7) for slit data
        # (38.5, 829.0) for slitless data
        x_center -= ta_dither_offset[0]
        y_center -= ta_dither_offset[1]
        log.debug(
            "Undoing TA verification image dither offsets: "
            f"x={ta_dither_offset[0]}, y={ta_dither_offset[1]}"
            f" -> new position: x={x_center}, y={y_center} should be near reference point."
        )

        # Apply science exposure subarray offset
        log.debug(
            f"Applying science exposure subarray origin offsets: xstart={xstart}, ystart={ystart}"
        )
        x_center -= xstart - 1
        y_center -= ystart - 1

        # Apply dither offset from science exposure
        _, (offset_x, offset_y) = find_dither_position(result)
        log.debug(f"Applying science exposure dither offsets: x={offset_x}, y={offset_y}")
        x_center += offset_x
        y_center += offset_y

        # Apply filter offsets
        with dm.FilteroffsetModel(filteroffset_file) as filteroffset:
            col_offset, row_offset = retrieve_filter_offset(filteroffset, ta_filter)
            log.debug(f"Applying filter offsets: column={col_offset}, row={row_offset}")
            x_center += col_offset
            y_center += row_offset

        # Store final source position in science model and set completion flag
        result.source_xpos = x_center
        result.source_ypos = y_center
        result.meta.cal_step.targ_centroid = "COMPLETE"
        log.info("TA centering complete.")
        log.info(
            "Final source position in science data frame: "
            f"({result.source_xpos:.2f}, {result.source_ypos:.2f})"
        )

        # Reconstruct the container if needed
        result = self._rebuild_container(container, result)
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

    def _rebuild_container(self, container, updated_sci_model):
        """
        Rebuild the container with the updated science model.

        Parameters
        ----------
        container : ModelContainer or DataModel
            The original container if present. If None, the updated science model
            is returned as is.
        updated_sci_model : DataModel
            The updated science data model.

        Returns
        -------
        new_container : ModelContainer
            The rebuilt container with the updated science model.
        """
        if container is None:
            return updated_sci_model
        sci_idx = container.ind_asn_type("science")
        container[sci_idx[0]] = updated_sci_model
        return container
