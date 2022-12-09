import logging
import numpy as np
from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class PixelReplacement:
    """Main class for performing pixel replacement.

    This class controls loading the input data model, selecting the
    method for pixel replacement, and executing each step. This class
    should provide modularization to allow for multiple options and possible
    future reference files.

    """

    # Shortcuts for DQ Flags
    DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
    REPLACED = datamodels.dqflags.pixel['RESERVED_4']

    # Shortcuts for dispersion direction for ease of reading
    HORIZONTAL = 1
    VERTICAL = 2

    default_suffix = 'pixrep'

    def __init__(self, input_model, **pars):
        """
        Initialize the class with input data model.

        Parameters
        ----------
        input_model : DataModel, str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        pars : dict, optional
            Optional parameters to modify how pixel replacement
            will execute.
        """
        self.input = input_model
        self.pars = dict()
        self.pars.update(pars)

        self.algorithm_dict = {
            'fit_profile': self.fit_profile,
        }

        self.algorithm = self.algorithm_dict[self.pars['algorithm']]

    def replace(self):
        """
        Process the input DataModel, unpack any model that holds
        more than one 2D spectrum, then apply selected algorithm
        to each 2D spectrum in input.
        """
        # ImageModel inputs (MIR_LRS-FIXEDSLIT)
        if isinstance(self.input, datamodels.ImageModel):
            self.output = self.algorithm(self.input)

        # MultiSlitModel inputs (WFSS, NRS_FIXEDSLIT, ?)
        elif isinstance(self.input, datamodels.MultiSlitModel):
            for i, slit in enumerate(self.input.slits):
                slit_model = datamodels.SlitModel(self.input.slits[i].instance)
                slit_replaced = self.algorithm(slit_model)
                self.input.slits[i] = slit_replaced

        # Unsure how to deal with IFU data - should this be run
        # on pre- or post- cube_build products?
        else:
            log.critical(f"For input of exposure type {self.input.exposure.type}\n"
                         f"Algorithm has not yet been implemented.")
            raise Exception

    def fit_profile(self, model):
        """
        Fit a profile to adjacent columns, scale profile to
        column with missing pixel(s), and find flux estimate
        from scaled profile.

        Parameters
        ----------
        model : DataModel
            Either the input to the pixel_replace step in the
            case of DataModels containing only one 2D spectrum,
            or a single 2D spectrum from the input DataModel
            containting multiple spectra (i.e. MultiSlitModel).
            Requires data and dq attributes.

        Returns
        -------
        model_replaced : DataModel
            DataModel with flagged bad pixels now flagged with
            TO-BE-DETERMINED and holding a flux value estimated
            from spatial profile, derived from adjacent columns.

        """
        dispaxis = model.meta.wcsinfo.dispersion_direction

        # Trucate array to region where good pixels exist
        good_pixels = np.where(~model.dq & 1)
        x_range = [np.min(good_pixels[0]), np.max(good_pixels[0]) + 1]
        y_range = [np.min(good_pixels[1]), np.max(good_pixels[1]) + 1]

        # Loop over axis of data array corresponding to cross-
        # dispersion direction by indexing data shape with
        # strange dispaxis argument.
        log.critical(f"Number of profiles: {model.data.shape[2 - dispaxis]}")
        for ind in range(model.data.shape[2 - dispaxis]):
            dq_slice = model.dq[self.custom_slice(dispaxis, ind)]
            n_bad = np.count_nonzero(dq_slice & self.DO_NOT_USE)
            if n_bad == len(dq_slice):
                log.debug(f"Slice {ind} contains no good pixels. Skipping replacement.")
            else:
                log.debug(f"Slice {ind} contains {len(dq_slice) - n_bad} good pixels.")



    def custom_slice(self, dispaxis, index):
        """
        Construct slice for ease of use with varying
        dispersion axis.

        Parameters
        ----------
        dispaxis : int
            Using module-defined HORIZONTAL=1,
            VERTICAL=2

        index : int or slice
            Index or indices of cross-dispersion
            vectors to slice

        Returns
        -------
        Tuple
            Slice constructed using np.s_
        """
        if dispaxis == self.HORIZONTAL:
            return np.s_[:, index]
        elif dispaxis == self.VERTICAL:
            return np.s_[index, :]
        else:
            raise Exception
