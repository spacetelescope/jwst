

from ..stpipe import Step
from . import badpix_selfcal
import numpy as np
from jwst import datamodels as dm
from jwst.master_background.master_background_step import split_container

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["BadpixSelfcalStep"]

class BadpixSelfcalStep(Step):
    """
    BadpixSelfcalStep: 
    """

    class_alias = "badpix_selfcal"
    bkg_suffix = "badpix_selfcal"

    spec = """
    flagfrac = float(default=0.001)  #fraction of pixels to flag on each of low and high end
    kernel_size = integer(default=15)  #size of kernel for median filter
    skip = boolean(default=True)
    """

    def process(self, input, bkg_list=None):
        """
        Flag CRs in the DQ array of a JWST exposure

        Parameters
        ----------
        input: JWST data model or association
            input science data to be corrected

        bkg_list: list of background ImageModels.
            When calwebb_spec2 is run, these are read from the association file
            inside calwebb_spec2.py. When run as a standalone step, these are either
            read directly from an input association file, or the user can provide
            this list manually

        Returns
        -------
        output: JWST data model or association
            data model with CRs flagged
        """
        with dm.open(input) as input_data:

            if isinstance(input_data, dm.ModelContainer):
                
                sci_models, bkg_list_asn = split_container(input_data)
                if bkg_list is None:
                    bkg_list = bkg_list_asn
                else:
                    log.warning("bkg_list provided as input, ignoring bkg_list from association file")
                
                # in calwebb_spec2 there should be only a single science exposure in an association
                input_sci = sci_models[0]

            elif isinstance(input_data, dm.IFUImageModel) or isinstance(input_data, dm.ImageModel):

                input_sci = input_data
                if bkg_list is None:
                    # true self-calibration on input data itself
                    bkg_list = [input_data]

            else:
                raise TypeError("Input data is not a ModelContainer, ImageModel, or IFUImageModel.\
                                Cannot continue.")
            
            # collapse background dithers into a single background model
            bkgd_3d = []
            for i, background_model in enumerate(bkg_list):
                bkgd_3d.append(background_model.data)
            bkgd = np.nanmin(np.asarray(bkgd_3d), axis=0)
            bad_indices = badpix_selfcal.badpix_selfcal(bkgd, self.flagfrac, self.kernel_size)

            # apply the flags to the science data
            input_sci = badpix_selfcal.apply_flags(input_sci, bad_indices)

            # apply the flags to the background data
            for i, background_model in enumerate(bkg_list):
                bkg_list[i] = badpix_selfcal.apply_flags(background_model, bad_indices)
            
            self.record_step_status(input_sci, "badpix_selfcal", success=True)
        return input_sci, bkg_list
