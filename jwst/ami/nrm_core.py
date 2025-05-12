import logging
import numpy as np
from . import lg_model
from . import utils
from . import oifits

from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class FringeFitter:
    """
    Fit fringes to get interferometric observables for the data.

    For the given information on the instrument and mask, calculate the
    fringe observables (visibilities and closure phases in the image plane.
    Original python was by A. Greenbaum & A. Sivaramakrishnan
    """

    def __init__(
        self,
        instrument_data,
        oversample=3,
        psf_offset_ff=None,
        npix="default",
        weighted=False,
    ):
        """
        Initialize the FringeFitter object.

        Parameters
        ----------
        instrument_data : jwst.ami.instrument_data.NIRISS object
            Information on the mask geometry (namely # holes), instrument,
            wavelength obs mode.
        oversample : int, optional
            Model oversampling (also how fine to measure the centering). Default is 3.
        psf_offset_ff : float, optional
            Subpixel centering of your data, if known. Default is None.
        npix : int, optional
            Number of data pixels to use. Default is to use the shape of the data frame.
        weighted : bool, optional
            If True, use Poisson variance for weighting, otherwise do not apply
            any weighting. Default is False.
        """
        self.instrument_data = instrument_data

        self.oversample = oversample
        self.psf_offset_ff = psf_offset_ff
        self.npix = npix
        self.weighted = weighted

        if self.weighted:
            log.info("leastsqnrm.weighted_operations() - weighted by Poisson variance")
        else:
            log.info("leastsqnrm.matrix_operations() - equally-weighted")

    def fit_fringes_all(self, input_model):
        """
        Generate the best model to match the data (centering, scaling, rotation).

        Parameters
        ----------
        input_model : instance Data Model
            DM object for input

        Returns
        -------
        output_model : AmiOIModel object
            AMI tables of median observables from LG algorithm fringe fitting in OIFITS format
        output_model_multi : AmiOIModel object
            AMI tables of observables for each integration
            from LG algorithm fringe fitting in OIFITS format
        lgfit : AmiLgFitModel object
            AMI cropped data, model, and residual data from LG algorithm fringe fitting

        Notes
        -----
        May allow parallelization by integration (later)
        """
        # scidata, dqmask are already centered around peak
        self.scidata, self.dqmask = self.instrument_data.read_data_model(input_model)
        self.instrument_data.isz = self.scidata.shape[1]

        # Initialize the output oifts models
        oifits_model = oifits.RawOifits(self.instrument_data)
        oifits_model.initialize_obsarrays()
        oifits_model_multi = oifits.RawOifits(self.instrument_data, method="multi")
        oifits_model_multi.initialize_obsarrays()

        # initialize the output lgfitmodel product arrays
        nslices = self.instrument_data.nslices
        # 3d arrays of centered data, models, and residuals (data - model)
        ctrd_arr = np.zeros((nslices, self.scidata.shape[1], self.scidata.shape[2]))
        n_ctrd_arr = np.zeros((nslices, self.scidata.shape[1], self.scidata.shape[2]))
        model_arr = np.zeros((nslices, self.scidata.shape[1], self.scidata.shape[2]))
        n_model_arr = np.zeros((nslices, self.scidata.shape[1], self.scidata.shape[2]))
        resid_arr = np.zeros((nslices, self.scidata.shape[1], self.scidata.shape[2]))
        n_resid_arr = np.zeros((nslices, self.scidata.shape[1], self.scidata.shape[2]))
        # Model parameters
        solns_arr = np.zeros((nslices, 44))

        for slc in range(nslices):
            log.info(f"Fitting fringes for iteration {slc} of {nslices}")
            nrmslc = self.fit_fringes_single_integration(slc)

            # populate the solutions of the lgfit model
            datapeak = nrmslc.reference.max()
            ctrd_arr[slc, :, :] = nrmslc.reference
            n_ctrd_arr[slc, :, :] = nrmslc.reference / datapeak
            model_arr[slc, :, :] = nrmslc.modelpsf
            n_model_arr[slc, :, :] = nrmslc.modelpsf / datapeak
            resid_arr[slc, :, :] = nrmslc.residual
            n_resid_arr[slc, :, :] = nrmslc.residual / datapeak
            solns_arr[slc, :] = nrmslc.soln

            # populate the oifits models
            oifits_model.populate_obsarray(slc, nrmslc)
            oifits_model_multi.populate_obsarray(slc, nrmslc)

        # Populate the LGFitModel with the output of the fringe fitting (LG algorithm).
        lgfit = datamodels.AmiLgFitModel()
        lgfit.centered_image = ctrd_arr
        lgfit.norm_centered_image = n_ctrd_arr
        lgfit.fit_image = model_arr
        lgfit.norm_fit_image = n_model_arr
        lgfit.resid_image = resid_arr
        lgfit.norm_resid_image = n_resid_arr
        lgfit.solns_table = np.recarray(
            solns_arr.shape[0],
            dtype=[("coeffs", "f8", solns_arr.shape[1])],
            buf=solns_arr,
        )

        # Now save final output model(s) of all slices, averaged slices to AmiOiModels
        output_model = oifits_model.make_oifits()
        output_model_multi = oifits_model_multi.make_oifits()

        return output_model, output_model_multi, lgfit

    def fit_fringes_single_integration(self, slc):
        """
        Generate the best model to match a single slice.

        Parameters
        ----------
        slc : int
            Index of the iteration to fit (0 to nslc-1).

        Returns
        -------
        nrm : LgModel object
            Model with best fit results for the given slice.

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
        nrm = lg_model.LgModel(
            self.instrument_data.nrm_model,
            bandpass=self.instrument_data.wls[0],
            mask=self.instrument_data.mask.maskname,
            pixscale=self.instrument_data.pscale_rad,
            holeshape=self.instrument_data.holeshape,
            affine2d=self.instrument_data.affine2d,
            over=self.oversample,
        )

        if self.npix == "default":
            self.npix = self.scidata[slc, :, :].shape[0]

        ctrd = self.scidata[slc]
        dqslice = self.dqmask[slc]

        nrm.reference = ctrd  # self.ctrd is the cropped image centered on the brightest pixel

        if self.psf_offset_ff is None:
            # returned values have offsets x-y flipped:
            # Finding centroids the Fourier way assumes no bad pixels case:
            # Fourier domain mean slope
            centroid = utils.find_centroid(ctrd)
            # centroid represents offsets from brightest pixel ctr
            # use flipped centroids to update centroid of image for JWST:
            # check parity for GPI, Vizier,...
            # pixel coordinates: - note the flip of [0] and [1] to match DS9 view
            nrm.xpos = centroid[1]  # flip 0 and 1 to convert
            nrm.ypos = centroid[0]  # flip 0 and 1
            nrm.psf_offset = nrm.xpos, nrm.ypos  # renamed .bestcenter to .psf_offset
        else:
            nrm.psf_offset = (
                self.psf_offset_ff
            )  # user-provided psf_offsetoffsets from array center are here.

        model = nrm.make_model(
            fov=ctrd.shape[0],
            psf_offset=nrm.psf_offset,
        )

        nrm.fit_image(
            ctrd,
            model_in=model,
            dqm=dqslice,
            weighted=self.weighted,
        )

        nrm.create_modelpsf()
        return nrm  # to fit_fringes_all, where output model is created from list of nrm objects
