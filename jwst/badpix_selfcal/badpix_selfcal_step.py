

from ..stpipe import Step
from . import badpix_selfcal
import numpy as np
from jwst import datamodels as dm
from jwst.master_background.master_background_step import split_container

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

    def process(self, input):
        """
        Flag CRs in the DQ array of a JWST exposure

        Parameters
        ----------
        input: JWST data model or association
            input science data to be corrected

        Returns
        -------
        output: JWST data model or association
            data model with CRs flagged
        """
        with dm.open(input) as input_data:

            if isinstance(input_data, dm.ModelContainer):
                # handle associations; this should be the majority of cases where the step is used
                sci_models, background_models = split_container(input_data)

                # collapse background dithers into a single background model
                bkgd_list = []
                for i, background_model in enumerate(background_models):
                    bkgd_list.append(background_model.data)
                bkgd = np.nanmin(np.asarray(bkgd_list), axis=0)
                bad_indices = badpix_selfcal.badpix_selfcal(bkgd, self.flagfrac, self.kernel_size)

                for model in input_data:
                    # apply bad indices to both science and background models
                    badpix_selfcal.apply_flags(model, bad_indices)

            elif isinstance(input_data, dm.IFUImageModel) or isinstance(input_data, dm.ImageModel):
                # true self-calibration on input data itself
                bkgd = input_data.data
                bad_indices = badpix_selfcal.badpix_selfcal(bkgd, self.flagfrac, self.kernel_size)
                input_data = badpix_selfcal.apply_flags(input_data, bad_indices)

            else:
                raise TypeError("Input data is not a ModelContainer, ImageModel, or IFUImageModel.\
                                Cannot continue.")
            
            self.record_step_status(input_data, "badpix_selfcal", success=True)
        return input_data
