"""Subtraction of background signal depending on the observing mode."""

import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.background.background_sub import background_sub
from jwst.background.background_sub_soss import subtract_soss_bkg
from jwst.background.background_sub_wfss import subtract_wfss_bkg
from jwst.stpipe import Step

__all__ = ["BackgroundStep"]

log = logging.getLogger(__name__)

WFSS_TYPES = ["NIS_WFSS", "NRC_GRISM", "NRC_WFSS"]


class BackgroundStep(Step):
    """Subtract background exposures from target exposures."""

    class_alias = "bkg_subtract"

    spec = """
        bkg_list = force_list(default=None)  # List of background files. Ignored for WFSS or if asn is provided
        save_combined_background = boolean(default=False)  # Save combined background image
        sigma = float(default=3.0)  # Clipping threshold
        maxiters = integer(default=None)  # Number of clipping iterations
        soss_source_percentile = float(default=35.0) # Threshold flux percentile to mask out source pixels
        soss_bkg_percentile = float_list(min=2, max=2, default=None) # Background percentiles to use; default is [25.0, 50.0]
        wfss_mmag_extract = float(default=None)  # WFSS minimum abmag to extract
        wfss_mask = string(default=None)  # WFSS source mask file
        wfss_maxiter = integer(default=5)  # WFSS iterative outlier rejection max iterations
        wfss_rms_stop = float(default=0)  # WFSS iterative outlier rejection RMS improvement threshold (percent)
        wfss_outlier_percent = float(default=1)  # WFSS outlier percentile to reject per iteration
    """  # noqa: E501

    # These reference files are only used for WFSS/GRISM or SOSS data.
    reference_file_types = ["bkg", "wavelengthrange"]

    # Define a suffix for optional saved output of the combined background
    bkg_suffix = "combinedbackground"

    def process(self, step_input, input_bkg_list=None):
        """
        Subtract designated background images from target exposures.

        Parameters
        ----------
        step_input : str, `~stdatamodels.jwst.datamodels.ImageModel` or \
                     `~stdatamodels.jwst.datamodels.IFUImageModel`
            Input target data model to which background subtraction is applied or asn file

        input_bkg_list : list, optional
            File name list of background exposures.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.ImageModel` or \
                 `~stdatamodels.jwst.datamodels.IFUImageModel`
            The background-subtracted target data model
        """
        asn = self.load_as_level2_asn(step_input)
        model, members_by_type = self._asn_get_data(asn)

        if model.meta.exposure.type in ["NIS_WFSS", "NRC_WFSS"]:
            # Get the reference file names
            bkg_name = self.get_reference_file(model, "bkg")
            wlrange_name = self.get_reference_file(model, "wavelengthrange")

            if bkg_name == "N/A":
                log.warning("No BKG reference file found. Skipping background subtraction.")
                model.meta.cal_step.bkg_subtract = "SKIPPED"
                return model

            log.info("Using BKG reference file %s", bkg_name)
            log.info("Using WavelengthRange reference file %s", wlrange_name)

            wfss_mask = None
            if self.wfss_mask is not None:
                with datamodels.ImageModel(self.wfss_mask) as mask_model:
                    if not mask_model.hasattr("mask"):
                        raise AttributeError(f"No 'mask' attribute found in {self.wfss_mask}")
                    wfss_mask = mask_model.mask.astype(bool)
                    if wfss_mask.shape != model.data.shape:
                        raise ValueError(
                            f"WFSS mask shape {wfss_mask.shape} does not match "
                            f"input data shape {model.data.shape}"
                        )

            # Do the background subtraction for WFSS/GRISM data
            rescaler_kwargs = {
                "p": self.wfss_outlier_percent,
                "maxiter": self.wfss_maxiter,
                "delta_rms_thresh": self.wfss_rms_stop / 100,
            }
            result = subtract_wfss_bkg(
                model,
                bkg_name,
                wlrange_name,
                self.wfss_mmag_extract,
                user_mask=wfss_mask,
                rescaler_kwargs=rescaler_kwargs,
            )
            if result is None:
                result = model
                result.meta.cal_step.bkg_subtract = "SKIPPED"
            elif not result.meta.cal_step.bkg_subtract:
                result.meta.cal_step.bkg_subtract = "COMPLETE"

        elif model.meta.exposure.type == "NIS_SOSS":
            # Fetch the background reference filename
            bkg_name = self.get_reference_file(model, "bkg")

            if bkg_name == "N/A":
                log.warning("No BKG reference file found. Skipping background subtraction.")
                model.meta.cal_step.bkg_subtract = "SKIPPED"
                return model

            log.info("Using BKG reference file %s", bkg_name)

            if self.soss_bkg_percentile is None:
                soss_bkg_percentile = [25.0, 50.0]
            else:
                soss_bkg_percentile = self.soss_bkg_percentile

            result = subtract_soss_bkg(
                model, bkg_name, self.soss_source_percentile, soss_bkg_percentile
            )
            if result is None:
                result = model
                result.meta.cal_step.bkg_subtract = "SKIPPED"
            else:
                result.meta.cal_step.bkg_subtract = "COMPLETE"

        else:
            # Get the background files to be subtracted
            if len(members_by_type["background"]) >= 1:
                bkg_list = members_by_type["background"]
            elif input_bkg_list is not None:
                if isinstance(input_bkg_list, str):
                    if "," in input_bkg_list:
                        bkg_list = input_bkg_list.split(sep=",")
                    else:
                        bkg_list = [input_bkg_list]
                else:
                    bkg_list = input_bkg_list
            else:
                bkg_list = self.bkg_list

            # Make sure to catch a trailing comma and ignore it
            if bkg_list is not None:
                bkg_list = [bg for bg in bkg_list if bg]

            # Make sure that the background list is not empty for this case,
            # or report and skip the step
            if bkg_list is None or len(bkg_list) == 0:
                log.warning("* No background list provided * Skipping step.")
                model.meta.cal_step.bkg_subtract = "SKIPPED"
                return model

            # check if input data is NRS_IFU
            tolerance = 1.0e-8
            do_sub = True
            if model.meta.instrument.name in ["NIRSPEC"]:
                # check if GWA_XTIL & GWA_YTIL values of source
                # background are the same. If not skip step
                input_xtilt = model.meta.instrument.gwa_xtilt
                input_ytilt = model.meta.instrument.gwa_ytilt
                for bkg_file in bkg_list:
                    with datamodels.open(bkg_file) as bkg_model:
                        bkg_xtilt = bkg_model.meta.instrument.gwa_xtilt
                        bkg_ytilt = bkg_model.meta.instrument.gwa_ytilt
                        if np.allclose(
                            (input_xtilt, input_ytilt),
                            (bkg_xtilt, bkg_ytilt),
                            atol=tolerance,
                            rtol=0,
                        ):
                            pass
                        else:
                            do_sub = False
                            break
            # Do the background subtraction
            if do_sub:
                bkg_model, result = background_sub(model, bkg_list, self.sigma, self.maxiters)
                result.meta.cal_step.bkg_subtract = "COMPLETE"
                if self.save_combined_background:
                    comb_bkg_path = self.save_model(bkg_model, suffix=self.bkg_suffix, force=True)
                    log.info(f"Combined background written to {comb_bkg_path}.")

            else:
                result = model
                result.meta.cal_step.bkg_subtract = "SKIPPED"
                log.warning("Skipping background subtraction")
                log.warning("GWA_XTIL and GWA_YTIL source values are not the same as bkg values")

        return result

    def _asn_get_data(self, asn):
        """
        Get the targets and catalog from an association.

        Parameters
        ----------
        asn : ``jwst.associations.lib.rules_level2_base.DMSLevel2bBase``
            Input association.

        Returns
        -------
        step_input : `~stdatamodels.jwst.datamodels.ImageModel` or \
                     `~stdatamodels.jwst.datamodels.IFUImageModel`
            Input target data model
        bkg_list : list
            File name list of background exposures
        """
        members_by_type = defaultdict(list)

        if len(asn["products"]) > 1:
            log.warning("Multiple products in input association. Using only the first one.")

        # Get the grism image and the catalog, direct image, and segmentation map
        exp_product = asn["products"][0]
        # Find all the member types in the product
        for member in exp_product["members"]:
            members_by_type[member["exptype"].lower()].append(member["expname"])

        # Get the science member. Technically there should only be one. Even if
        # there are more, we'll just get the first one found.
        science_member = members_by_type["science"]
        if len(science_member) != 1:
            log.warning(
                "Wrong number of science exposures found in {}".format(exp_product["name"])  # noqa: E501
            )
            log.warning("    Using only first one.")

        science_member = science_member[0]
        log.info("Working on input %s ...", science_member)

        # Open the datamodel and update it with the relevant info for the background step
        sci = self.prepare_output(science_member)
        exp_type = sci.meta.exposure.type
        if exp_type in WFSS_TYPES:
            try:
                sci.meta.source_catalog = Path(members_by_type["sourcecat"][0]).name
                log.info(f"Using sourcecat file {sci.meta.source_catalog}")
                sci.meta.segmentation_map = Path(members_by_type["segmap"][0]).name
                log.info(f"Using segmentation map {sci.meta.segmentation_map}")
                sci.meta.direct_image = Path(members_by_type["direct_image"][0]).name
                log.info(f"Using direct image {sci.meta.direct_image}")
            except IndexError:
                if sci.meta.source_catalog is None:
                    raise IndexError(
                        "No source catalog specified in association or datamodel."
                    ) from None

        return sci, members_by_type
