#! /usr/bin/env python

# Module for fringe fitting class
# Original python was by A. Greenbaum & A. Sivaramakrishnan

import logging
import numpy as np
from . import lg_model
from . import utils
from . import oifits

from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class FringeFitter:
    def __init__(self, instrument_data, **kwargs):
        """
        For the given information on the instrument and mask, calculate the
        fringe observables (visibilities and closure phases in the image plane.

        Parameters
        ----------
        instrument_data - jwst.ami.instrument_data.NIRISS object
            Information on the mask geometry (namely # holes), instrument,
            wavelength obs mode.

        kwargs options:
            oversample - model oversampling (also how fine to measure the centering)
            psf_offset - subpixel centering of your data, if known
            npix - number of data pixels to use. Default is the shape of the data frame.
            find_rotation - will find the best pupil rotation that matches the data
        """
        self.instrument_data = instrument_data

        # Options
        if "oversample" in kwargs:
            self.oversample = kwargs["oversample"]
        else:
            # default oversampling is 3
            self.oversample = 3

        if "find_rotation" in kwargs:
            # can be True/False or 1/0
            self.find_rotation = kwargs["find_rotation"]
        else:
            self.find_rotation = False

        if "psf_offset_ff" in kwargs:  # if so do not find center of image in data
            self.psf_offset_ff = kwargs["psf_offset_ff"]
        else:
            self.psf_offset_ff = None  # find center of image in data

        if "npix" in kwargs:
            self.npix = kwargs["npix"]
        else:
            self.npix = 'default'
        # Default: unweighted fit
        self.weighted = False    
        if "weighted" in kwargs:
            self.weighted = kwargs["weighted"]
        if self.weighted is True:
            log.info("leastsqnrm.weighted_operations() - weighted by Poisson variance")
        else:
            log.info("leastsqnrm.matrix_operations() - equally-weighted")

    def fit_fringes_all(self, input_model):
        """
        Extract the input data from input_model, and generate the best model to
        match the data (centering, scaling, rotation)
        May allow parallelization by integration (later)

        Parameters
        ----------
        input_model: instance Data Model
            DM object for input

        Returns
        -------
        output_model: AmiOIModel object
            AMI tables of median observables from LG algorithm fringe fitting in OIFITS format
        output_model_multi: AmiOIModel object
            AMI tables of observables for each integration from LG algorithm fringe fitting in OIFITS format
        lgfit:
            AMI cropped data, model, and residual data from LG algorithm fringe fitting
        """

        # scidata, dqmask are already centered around peak
        self.scidata, self.dqmask = self.instrument_data.read_data_model(input_model)

        # list for nrm objects for each slc
        self.nrm_list = []

        for slc in range(self.instrument_data.nwav):
            self.nrm_list.append(self.fit_fringes_single_integration(slc))


        # Now save final output model(s) of all slices, averaged slices to AmiOiModels
        # averaged
        oifits_model = oifits.RawOifits(self)
        output_model = oifits_model.make_oifits()

        # multi-integration
        oifits_model_multi = oifits.RawOifits(self,method='multi')
        output_model_multi = oifits_model_multi.make_oifits()

        # Save cropped/centered data, model, residual in AmiLgFitModel
        lgfit = self.make_lgfitmodel()

        return output_model, output_model_multi, lgfit
        
    def make_lgfitmodel(self):
        """
        Populate the LGFitModel with the output of the fringe fitting
        (LG algorithm)

        Parameters
        ----------

        Returns
        -------
        m: AmiLgFitModel object
            LG analysis centered data, fit, residual, and model info
        """
        nslices = len(self.nrm_list)
        # 3d arrays of centered data, models, and residuals (data - model)
        ctrd_arr = np.zeros((nslices,self.scidata.shape[1],self.scidata.shape[2]))
        n_ctrd_arr = np.zeros((nslices,self.scidata.shape[1],self.scidata.shape[2]))
        model_arr = np.zeros((nslices,self.scidata.shape[1],self.scidata.shape[2]))
        n_model_arr = np.zeros((nslices,self.scidata.shape[1],self.scidata.shape[2]))
        resid_arr = np.zeros((nslices,self.scidata.shape[1],self.scidata.shape[2]))
        n_resid_arr = np.zeros((nslices,self.scidata.shape[1],self.scidata.shape[2]))
        # Model parameters
        solns_arr = np.zeros((nslices,44))

        for i,nrmslc in enumerate(self.nrm_list):
            datapeak = nrmslc.reference.max()
            ctrd_arr[i,:,:] = nrmslc.reference
            n_ctrd_arr[i,:,:] = nrmslc.reference/datapeak
            model_arr[i,:,:] = nrmslc.modelpsf
            n_model_arr[i,:,:] = nrmslc.modelpsf/datapeak
            resid_arr[i,:,:] = nrmslc.residual
            n_resid_arr[i,:,:] = nrmslc.residual/datapeak
            solns_arr[i,:] = nrmslc.soln

        # Populate datamodel
        m = datamodels.AmiLgFitModel()
        m.centered_image = ctrd_arr
        m.norm_centered_image = n_ctrd_arr
        m.fit_image = model_arr
        m.norm_fit_image = n_model_arr
        m.resid_image = resid_arr
        m.norm_resid_image = n_resid_arr
        m.solns_table = np.recarray(solns_arr.shape[0], dtype=[('coeffs', 'f8', solns_arr.shape[1])], buf=solns_arr)

        return m


    def fit_fringes_single_integration(self, slc):
        """
        Generate the best model to
        match a single slice

        Parameters
        ----------
        slc: numpy array
            2D slice of data

        Returns
        -------
        nrm: NrmModel object
            Model with best fit results
        
        Notes
        -----
        After nrm.fit_image is called, these attributes are stored in nrm object:

        -----------------------------------------------------------------------------
        soln            --- resulting sin/cos coefficients from least squares fitting
        fringephase     --- baseline phases in radians
        fringeamp       --- baseline amplitudes (flux normalized)
        redundant_cps   --- closure phases in radians
        redundant_cas   --- closure amplitudes
        residual        --- fit residuals [data - model solution]
        cond            --- matrix condition for inversion
        fringepistons   --- zero-mean piston opd in radians on each hole (eigenphases)
        -----------------------------------------------------------------------------
        
        """

        nrm = lg_model.NrmModel(mask=self.instrument_data.mask,
                                pixscale=self.instrument_data.pscale_rad,
                                holeshape=self.instrument_data.holeshape,
                                affine2d=self.instrument_data.affine2d,
                                over=self.oversample)

        nrm.bandpass = self.instrument_data.wls[0]

        if self.npix == 'default':
            self.npix = self.scidata[slc,:, :].shape[0]

        self.ctrd = self.scidata[slc]
        self.dqslice = self.dqmask[slc]

        nrm.reference = self.ctrd  # self.ctrd is the cropped image centered on the brightest pixel

        if self.psf_offset_ff is None:
            # returned values have offsets x-y flipped:
            # Finding centroids the Fourier way assumes no bad pixels case - Fourier domain mean slope
            centroid = utils.find_centroid(self.ctrd, self.instrument_data.threshold)  # offsets from brightest pixel ctr
            # use flipped centroids to update centroid of image for JWST - check parity for GPI, Vizier,...
            # pixel coordinates: - note the flip of [0] and [1] to match DS9 view
            nrm.xpos = centroid[1]  # flip 0 and 1 to convert
            nrm.ypos = centroid[0]  # flip 0 and 1
            nrm.psf_offset = nrm.xpos, nrm.ypos  # renamed .bestcenter to .psf_offset
        else:
            nrm.psf_offset = self.psf_offset_ff  # user-provided psf_offsetoffsets from array center are here.

        nrm.make_model(fov=self.ctrd.shape[0], 
                       bandpass=nrm.bandpass,
                       over=self.oversample,
                       psf_offset=nrm.psf_offset,
                       pixscale=nrm.pixel)

        nrm.fit_image(self.ctrd, 
                      modelin=nrm.model,
                      psf_offset=nrm.psf_offset,
                      dqm=self.dqslice,
                      weighted=self.weighted)

        nrm.create_modelpsf()
        # model now stored as nrm.modelpsf, also nrm.residual
        self.nrm = nrm # this gets updated with each slice
        return nrm # to fit_fringes_all, where the output model will be created from list of nrm objects
