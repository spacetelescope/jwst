#! /usr/bin/env python

# Original python by A. Greenbaum & A. Sivaramakrishnan 

"""
Contains 

FringeFitter - fit fringe phases and amplitudes to data in the image plane

Calibrate - Calibrate target data with calibrator data <<< DELETE

DELETE THIS:
BinaryAnalyze - Detection, mcmc modeling, visualization tools
                Much of this section is based on tools in the pymask code
                by F. Martinache, B. Pope, and A. Cheetham
                We especially thank A. Cheetham for help and advising for
                developing the analysis tools in this package. 

    LG++ anand@stsci.edu nrm_core changes:
        Removed use of 'centering' parameter, switch to psf_offsets, meant to be uniformly ds9-compatible 
        offsets from array center (regardless of even/odd array sizes).

            nrm_core.fit_image(): refslice UNTESTED w/ new utils.centroid()
            nrm_core.fit_image(): hold_centering UNTESTED w/ new utils.centroid()

"""

import os, sys, time
import numpy as np
from astropy.io import fits
from scipy.special import comb
from scipy.stats import sem, mstats
import pickle as pickle

#import nrm_analysis.oifits as oifits #### 080320 NEEDED ??

from nrm_analysis.fringefitting.LG_Model import NRM_Model
from nrm_analysis.misctools import utils  # AS LG++
from nrm_analysis.modeling.binarymodel import model_cp_uv, model_allvis_uv, model_v2_uv, model_t3amp_uv
from nrm_analysis.modeling.multimodel import model_bispec_uv

from multiprocessing import Pool

class FringeFitter:
    def __init__(self, instrument_data, **kwargs):
        """
        Fit fringes in the image plane

        Takes an instance of the appropriate instrument class
        Various options can be set

        kwarg options:
        oversample - model oversampling (also how fine to measure the centering)
        psf_offset - If you already know the subpixel centering of your data, give it here
                     (not recommended except when debugging with perfectly know image placement))
        savedir - Where do you want to save the new files to? Default is working directory.
        #datadir - Where is your data? Default is working directory.
        npix - How many pixels of your data do you want to use? 
               Default is the shape of a data [slice or frame].  Typically odd?
        debug - will plot the FT of your data next to the FT of a reference PSF.
                Needs poppy package to run
        verbose_save - saves more than the standard files
        interactive - default True, prompts user to overwrite/create fresh directory.  
                      False will overwrite files where necessary.
        find_rotation - will find the best pupil rotation that matches the data

        main method:
        * fit_fringes


        """
        self.instrument_data = instrument_data

        # Options
        if "oversample" in kwargs:
            self.oversample = kwargs["oversample"]
        else:
            #default oversampling is 3
            self.oversample = 3

        if "find_rotation" in kwargs:
            # can be True/False or 1/0
            self.find_rotation = kwargs["find_rotation"]
        else:
            self.find_rotation = False

        if "psf_offset_ff" in kwargs: # if so do not find center of image in data
            self.psf_offset_ff = kwargs["psf_offset_ff"]         
        else:
            self.psf_offset_ff = None #find center of image in data

        if "savedir" in kwargs:
            self.savedir = kwargs["savedir"]
        else:
            self.savedir = os.getcwd()

        if "npix" in kwargs:
            self.npix = kwargs["npix"]
        else:
            self.npix = 'default'

        if "debug" in kwargs:
            self.debug=kwargs["debug"]
        else:
            self.debug=False

        if "verbose_save" in kwargs:
            self.verbose_save = kwargs["verbose_save"]
        else:
            self.verbose_save = False

        if 'interactive' in kwargs:
            self.interactive = kwargs['interactive']
        else:
            self.interactive = True

        if "save_txt_only" in kwargs:
            self.save_txt_only = kwargs["save_txt_only"]
        else:
            self.save_txt_only = False

        # Create directories if they don't already exit
        try:
            os.mkdir(self.savedir)
        except:
            if self.interactive is True:
                ans = input()
                if ans == "y":
                    pass
                elif ans == "n":
                    sys.exit("use alternative save directory with kwarg 'savedir' when calling FringeFitter")
                else:
                    sys.exit("Invalid answer. Stopping.")
            else:
                pass

    # May 2017 J Sahlmann updates: parallelized fringe-fitting!

    def fit_fringes(self, fns, threads = 0):
        if type(fns) == str:
            fns = [fns, ]

        # Can get fringes for images in parallel
        store_dict = [{"object":self, "file":                 fn,"id":jj} \
                        for jj,fn in enumerate(fns)]

        t2 = time.time()
        for jj, fn in enumerate(fns):
            fit_fringes_parallel({"object":self, "file":                  fn,\
                                  "id":jj},threads)
        t3 = time.time()


    def save_output(self, slc, nrm):
        # cropped & centered PSF
        self.datapeak = self.ctrd.max()
        #TBD: Keep only n_*.fits files after testing is done and before doing ImPlaneIA delivery
        if self.save_txt_only==False:
            fits.PrimaryHDU(data=self.ctrd, \
                    header=self.scihdr).writeto(self.savedir+\
                    self.sub_dir_str+"/centered_{0}.fits".format(slc), \
                    overwrite=True)
            fits.PrimaryHDU(data=self.ctrd/self.datapeak, \
                    header=self.scihdr).writeto(self.savedir+\
                    self.sub_dir_str+"/n_centered_{0}.fits".format(slc), \
                    overwrite=True)

            model, modelhdu = nrm.plot_model(fits_true=1)
            # save to fits files
            fits.PrimaryHDU(data=nrm.residual).writeto(self.savedir+\
                        self.sub_dir_str+"/residual_{0:02d}.fits".format(slc), \
                        overwrite=True)
            fits.PrimaryHDU(data=nrm.residual/self.datapeak).writeto(self.savedir+\
                        self.sub_dir_str+"/n_residual_{0:02d}.fits".format(slc), \
                        overwrite=True)
            modelhdu.writeto(self.savedir+\
                        self.sub_dir_str+"/modelsolution_{0:02d}.fits".format(slc),\
                        overwrite=True)
            fits.PrimaryHDU(data=model/self.datapeak, \
                    header=modelhdu.header).writeto(self.savedir+\
                    self.sub_dir_str+"/n_modelsolution_{0:02d}.fits".format(slc), \
                    overwrite=True)

        # default save to text files
        np.savetxt(self.savedir+self.sub_dir_str + \
                   "/solutions_{0:02d}.txt".format(slc), nrm.soln)
        np.savetxt(self.savedir+self.sub_dir_str + \
                   "/phases_{0:02d}.txt".format(slc), nrm.fringephase)
        np.savetxt(self.savedir+self.sub_dir_str + \
                   "/amplitudes_{0:02d}.txt".format(slc), nrm.fringeamp)
        np.savetxt(self.savedir+self.sub_dir_str + \
                   "/CPs_{0:02d}.txt".format(slc), nrm.redundant_cps)
        np.savetxt(self.savedir+self.sub_dir_str + \
                   "/CAs_{0:02d}.txt".format(slc), nrm.redundant_cas)
        np.savetxt(self.savedir+self.sub_dir_str + \
                  "/fringepistons_{0:02d}.txt".format(slc), nrm.fringepistons)

        # write info that oifits wants only when writing out first slice.
        # this will relevant to all slices... so no slice number here.
        if slc == 0:
            pfn = self.savedir+self.sub_dir_str + "/info4oif_dict.pkl"
            pfd = open(pfn, 'wb')
            pickle.dump(self.instrument_data.info4oif_dict, pfd)
            pfd.close()


        # optional save outputs
        if self.verbose_save:
            np.savetxt(self.savedir+self.sub_dir_str+\
                       "/condition_{0:02d}.txt".format(slc), nrm.cond)
            np.savetxt(self.savedir+self.sub_dir_str+\
                       "/flux_{0:02d}.txt".format(slc), nrm.flux)
          
        if nrm.linfit_result is not None:          
            # save linearfit results to pickle file
            myPickleFile = os.path.join(self.savedir+self.sub_dir_str,"linearfit_result_{0:02d}.pkl".format(slc))
            with open( myPickleFile , "wb" ) as f:
                pickle.dump((nrm.linfit_result), f) 
                       

    def save_auto_figs(self, slc, nrm):
        
        # rotation
        if self.find_rotation==True:
            plt.figure()
            plt.plot(nrm.rots, nrm.corrs)
            plt.vlines(nrm.rot_measured, nrm.corrs[0],
                        nrm.corrs[-1], linestyles='--', color='r')
            plt.text(nrm.rots[1], nrm.corrs[1], 
                     "best fit at {0}".format(nrm.rot_measured))
            plt.savefig(self.savedir+self.sub_dir_str+\
                        "/rotationcorrelation_{0:02d}.png".format(slc))

def fit_fringes_parallel(args, threads):
    self = args['object']
    filename = args['file']
    id_tag = args['id']
    self.prihdr, self.scihdr, self.scidata = self.instrument_data.read_data(filename)
    self.sub_dir_str = self.instrument_data.sub_dir_str
    try:
        os.mkdir(self.savedir+self.sub_dir_str)
    except:
        pass

    store_dict = [{"object":self, "slc":slc} for slc in \
                  range(self.instrument_data.nwav)]

    if threads>0:
        pool = Pool(processes=threads)
        pool.map(fit_fringes_single_integration, store_dict)
        pool.close()
        pool.join()

    else:
        for slc in range(self.instrument_data.nwav):
            fit_fringes_single_integration({"object":self, "slc":slc})

def fit_fringes_single_integration(args):
    self = args["object"]
    slc = args["slc"]
    id_tag = args["slc"]
    nrm = NRM_Model(mask=self.instrument_data.mask,
                    pixscale=self.instrument_data.pscale_rad,
                    holeshape=self.instrument_data.holeshape,
                    affine2d=self.instrument_data.affine2d,
                    over = self.oversample)

    nrm.bandpass = self.instrument_data.wls[slc]

    if self.npix == 'default':
        self.npix = self.scidata[slc,:,:].shape[0]

    DBG = False # AS testing gross psf orientatioun while getting to LG++ beta release 2018 09
    if DBG:
        nrm.simulate(fov=self.npix, bandpass=self.instrument_data.wls[slc], over=self.oversample)
        fits.PrimaryHDU(data=nrm.psf).writeto(self.savedir + "perfect.fits", overwrite=True)

    # New or modified in LG++
    # center the image on its peak pixel:
    # AS subtract 1 from "r" below  for testing >1/2 pixel offsets
    # AG 03-2019 -- is above comment still relevant?
    
    if self.instrument_data.arrname=="NIRC2_9NRM":
        self.ctrd = utils.center_imagepeak(self.scidata[slc, :,:], 
                        r = (self.npix -1)//2 - 2, cntrimg=False)  
    elif self.instrument_data.arrname=="gpi_g10s40":
        self.ctrd = utils.center_imagepeak(self.scidata[slc, :,:], 
                        r = (self.npix -1)//2 - 2, cntrimg=True)  
    else:
        self.ctrd = utils.center_imagepeak(self.scidata[slc, :,:])  
    

    nrm.reference = self.ctrd  # self. is the cropped image centered on the brightest pixel
    if self.psf_offset_ff is None:
        # returned values have offsets x-y flipped:
        # Finding centroids the Fourier way assumes no bad pixels case - Fourier domain mean slope
        centroid = utils.find_centroid(self.ctrd, self.instrument_data.threshold) # offsets from brightest pixel ctr
        # use flipped centroids to update centroid of image for JWST - check parity for GPI, Vizier,...
        # pixel coordinates: - note the flip of [0] and [1] to match DS9 view
        image_center = utils.centerpoint(self.ctrd.shape) + np.array((centroid[1], centroid[0])) # info only, unused

        nrm.xpos = centroid[1]  # flip 0 and 1 to convert
        nrm.ypos = centroid[0]  # flip 0 and 1
        nrm.psf_offset = nrm.xpos, nrm.ypos  # renamed .bestcenter to .psf_offset

    else:
        nrm.psf_offset = self.psf_offset_ff # user-provided psf_offset python-style offsets from array center are here.


    nrm.make_model(fov = self.ctrd.shape[0], bandpass=nrm.bandpass, 
                   over=self.oversample,
                   psf_offset=nrm.psf_offset,  
                   pixscale=nrm.pixel)
    nrm.fit_image(self.ctrd, modelin=nrm.model, psf_offset=nrm.psf_offset)

    """
    Attributes now stored in nrm object:

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

    if self.debug==True:
        import matplotlib.pyplot as plt
        import poppy.matrixDFT as mft
        dataft = mft.matrix_dft(self.ctrd, 256, 512)
        refft = mft.matrix_dft(self.refpsf, 256, 512)
        plt.figure()
        plt.title("Data")
        plt.imshow(np.sqrt(abs(dataft)), cmap = "bone")
        plt.figure()
        plt.title("Reference")
        plt.imshow(np.sqrt(abs(refft)), cmap="bone")
        plt.show()
    
    self.save_output(slc, nrm)
    self.nrm = nrm # store  extracted values here
    return None


def cp_binary_model(params, constant, priors, spectrum_model, uvcoords, cp, cperr, stat="loglike"):
    # really want to be able to give this guy some general oi_data and have bm() sort it out.
    # Need to figure out how to add in the priors

    ##################################################
    # HOW DO I TUNE THIS DEPENDING ON MY OBSERVATIONS? - need a keyword or something, need help.
    # data = self.cp, self.cperr#, self.v2, self.v2err
    ##################################################

    # priors, i.e. bounds here
    # constant
    # uvcoords
    # cp
    # cperr

    for i in range(len(params)):
        if (params[i] < priors[i][0] or params[i] > priors[i][1]):  
            return -np.inf

    if spectrum_model == None:

        # Model from params
        #model_cp = model_cp_uv(self.uvcoords, params['con'], params['sep'], \
        #                   params['pa'], 1.0/self.constant['wavl'])
        model_cp = model_cp_uv(uvcoords, params[0], params[1], \
                            params[2], 1.0/constant['wavl'])
    elif spectrum_model == 'slope':
        # params needs 'con_start' starting contrast and 'slope,' sep & pa constant?
        wav_step = constant['wavl'][1] - constant['wavl'][0]
        # contrast model is con_start + slope*delta_lambda
        #contrast = params[0] + params[1]*wav_step
        band_diff = constant['wavl'][-1] - constant['wavl'][0]
        contrast = np.linspace(params[0], params[0]+(params[1]*band_diff), \
                               num = len(constant['wavl']))
        # Model from params
        model_cp = model_cp_uv(uvcoords, contrast, params[2], \
                            params[3], 1.0/constant['wavl'])
    elif spectrum_model == 'free' :
        # Model from params - params is contrast array nwav long, sep & pa constant
        model_cp = model_cp_uv(uvcoords, params, constant['sep'], \
                            constant['pa'], 1.0/constant['wavl'])
    else:
        sys.exit("Invalid spectrum model")

    chi2stat = reduced_chi2(cp, cperr, model_cp, 34.0)
    ll = logl(cp, cperr, model_cp)
    if stat == "loglike":
        return ll
    elif stat == "chi2":
        return chi2stat

def v2_binary_model(params, constant, priors, spectrum_model, uvcoords, v2, v2err, stat="loglike"):
    # really want to be able to give this guy some general oi_data and have bm() sort it out.
    # Need to figure out how to add in the priors

    ##################################################
    # HOW DO I TUNE THIS DEPENDING ON MY OBSERVATIONS? - need a keyword or something, need help.
    # data = self.cp, self.cperr#, self.v2, self.v2err
    ##################################################

    # priors, i.e. bounds here
    # constant
    # uvcoords
    # cp
    # cperr

    for i in range(len(params)):
        if (params[i] < priors[i][0] or params[i] > priors[i][1]):  
            return -np.inf

    if spectrum_model == None:

        # Model from params
        #model_cp = model_cp_uv(self.uvcoords, params['con'], params['sep'], \
        #                   params['pa'], 1.0/self.constant['wavl'])
        model_v2 = model_v2_uv(uvcoords, params[0], params[1], \
                            params[2], 1.0/constant['wavl'])
    elif spectrum_model == 'slope':
        # params needs 'con_start' starting contrast and 'slope,' sep & pa constant?
        wav_step = constant['wavl'][1] - constant['wavl'][0]
        # contrast model is con_start + slope*delta_lambda
        contrast = params[0] + params[1]*wav_step
        # Model from params
        model_v2 = model_v2_uv(uvcoords, contrast, params[2], \
                            params[3], 1.0/constant['wavl'])
    elif spectrum_model == 'free' :
        # Model from params - params is contrast array nwav long, sep & pa constant
        model_v2 = model_v2_uv(uvcoords, params, constant['sep'], \
                            constant['pa'], 1.0/constant['wavl'])
    else:
        sys.exit("Invalid spectrum model")

    chi2stat = reduced_chi2(v2, v2err, model_v2, 34.0)
    ll = logl(v2, v2err, model_v2)
    if stat == "loglike":
        return ll
    elif stat == "chi2":
        return chi2stat

def allvis_binary_model(params, constant, priors, spectrum_model, uvcoords, \
                        cp, cperr, vis, viserr, stat="loglike", dof = 1):
    """
    For now, writing a new function to do this with both cps and bisectrum visibilities
    """
    for i in range(len(params)):
        if (params[i] < priors[i][0] or params[i] > priors[i][1]):  
            return -np.inf

    if spectrum_model == None:

        # Model from params
        model_cp = model_cp_uv(uvcoords, params[0], params[1], \
                           params[2], 1.0/constant['wavl'])
        model_vis = model_t3amp_uv(uvcoords, params[0], params[1], \
                           params[2], 1.0/constant['wavl'])
        #model_cp, model_vis = model_allvis_uv(uvcoords, uvcoords_vis, params[0], params[1], \
        #                    params[2], 1.0/constant['wavl'])

        model = np.zeros((model_cp.shape[0] + model_vis.shape[0], model_cp.shape[1]))
        model[:model_cp.shape[0],...] = model_cp
        model[model_cp.shape[0]:,...] = model_vis
    else:
        sys.exit("Invalid spectrum model")

    allvisobs = np.zeros((cp.shape[0]+vis.shape[0], cp.shape[1]))
    allvisobserr = np.zeros((cp.shape[0]+vis.shape[0], cp.shape[1]))
    allvisobs[:cp.shape[0], ...] = cp
    allvisobs[cp.shape[0]:, ...] = vis
    allvisobserr[:cp.shape[0], ...] = cperr
    allvisobserr[cp.shape[0]:, ...] = viserr

    if stat == "loglike":
        ll = logl(allvisobs, allvisobserr, model)
        return ll
    elif stat == "chi2":
        chi2stat = reduced_chi2(allvisobs, allvisobserr, model, dof)
        return chi2stat

def cp_multi_model(params, constant, priors, spectrum_model, uvcoords, \
                        cp, cperr, stat="loglike", dof = 1):
    """
    For now, writing a new function to do this with both cps and bisectrum visibilities
    """
    for i in range(len(params)):
        if (params[i] < priors[i%3][0] or params[i] > priors[i%3][1]):  
            return -np.inf
    # chop up the params by type input: [c1, s1, p1, c2, s2, p2, ...]
    nsource = len(params) / 3
    cons = params[0::3]
    seps = params[1::3]
    pas = params[2::3]
    if spectrum_model == None:

        # Model from params
        model_cp, model_vis = model_bispec_uv(uvcoords, cons, seps, \
                           pas, 1.0/constant['wavl'])
        #model_cp = np.zeros(cp.shape)
        #for q in range(nsource):
        #    model_cp += model_cp_uv(uvcoords, cons[q], seps[q], pas[q], 1.0/constant['wavl'])

        #model = np.zeros((model_cp.shape[0] + model_vis.shape[0], model_cp.shape[1]))
        #model[:model_cp.shape[0],...] = model_cp
        #model[model_cp.shape[0]:,...] = model_vis
    else:
        sys.exit("Invalid spectrum model")

    """
    allvisobs = np.zeros((cp.shape[0]+vis.shape[0], cp.shape[1]))
    allvisobserr = np.zeros((cp.shape[0]+vis.shape[0], cp.shape[1]))
    allvisobs[:cp.shape[0], ...] = cp
    allvisobs[cp.shape[0]:, ...] = vis
    allvisobserr[:cp.shape[0], ...] = cperr
    allvisobserr[cp.shape[0]:, ...] = viserr
    """

    if stat == "loglike":
        ll = logl(cp, cperr, model_cp)
        return ll
    elif stat == "chi2":
        chi2stat = reduced_chi2(cp, cperr, model_cp, dof)
        return chi2stat

def bispec_multi_model(params, constant, priors, spectrum_model, uvcoords, \
                        cp, cperr, vis, viserr, stat="loglike", dof = 1):
    """
    For now, writing a new function to do this with both cps and bisectrum visibilities
    """
    for i in range(len(params)):
        if (params[i] < priors[i%3][0] or params[i] > priors[i%3][1]):  
            return -np.inf

    # chop up the params by type input: [c1, s1, p1, c2, s2, p2, ...]
    nsource = len(params) / 3
    cons = params[0::3]
    seps = params[1::3]
    pas = params[2::3]
    if spectrum_model == None:

        # Model from params
        model_cp, model_vis = model_bispec_uv(uvcoords, cons, seps, \
                           pas, 1.0/constant['wavl'])

        model = np.zeros((model_cp.shape[0] + model_vis.shape[0], model_cp.shape[1]))
        model[:model_cp.shape[0],...] = model_cp
        model[model_cp.shape[0]:,...] = model_vis
    else:
        sys.exit("Invalid spectrum model")

    allvisobs = np.zeros((cp.shape[0]+vis.shape[0], cp.shape[1]))
    allvisobserr = np.zeros((cp.shape[0]+vis.shape[0], cp.shape[1]))
    allvisobs[:cp.shape[0], ...] = cp
    allvisobs[cp.shape[0]:, ...] = vis
    allvisobserr[:cp.shape[0], ...] = cperr
    allvisobserr[cp.shape[0]:, ...] = viserr

    if stat == "loglike":
        ll = logl(allvisobs, allvisobserr, model)
        return ll
    elif stat == "chi2":
        chi2stat = reduced_chi2(allvisobs, allvisoberr, model, dof)
        return chi2stat

def get_data(self):
    # Move this function out, pass values to the object
    try:
        self.oifdata = oifits.open(self.oifitsfn)

    try:
        self.avparang = self.oifdata.avparang
        self.parang_range = self.oifdata.parang_range

    print("Apparently we opened ", self.oifitsfn, "successfully")

    self.telescope = list(self.oifdata.wavelength.keys())[0]
    self.ncp = len(self.oifdata.t3)
    self.nbl = len(self.oifdata.vis2)
    self.wavls = self.oifdata.wavelength[self.telescope].eff_wave
    self.eff_band = self.oifdata.wavelength[self.telescope].eff_band
    self.nwav = len(self.wavls)
    self.uvcoords = np.zeros((2, 3, self.ncp))
    self.uvcoords_vis = np.zeros((2, self.nbl))

    # Now collect fringe observables and coordinates
    self.cp = np.zeros((self.ncp, self.nwav))
    self.cperr = np.zeros((self.ncp, self.nwav))
    self.t3amp = np.zeros((self.ncp, self.nwav))
    self.t3amperr = np.zeros((self.ncp, self.nwav))
    self.v2 = np.zeros((self.nbl, self.nwav))
    self.v2err = np.zeros((self.nbl, self.nwav))
    self.pha = np.zeros((self.nbl, self.nwav))
    self.phaerr = np.zeros((self.nbl, self.nwav))

    #Now, if extra_error is specified and there is a wavelength axis, we scale the 
    # estimated error with wavelength - User specifies error at shorted wavl?
    #self.extra_error = self.extra_error*np.ones(self.nwav)*self.wavls[0] / (self.wavls)

    for ii in range(self.ncp):
        self.cp[ii, :] = self.oifdata.t3[ii].t3phi
        self.cperr[ii, :] = np.sqrt(self.oifdata.t3[ii].t3phierr**2 + self.extra_error**2)
        self.t3amp[ii, :] = self.oifdata.t3[ii].t3amp
        self.t3amperr[ii, :] = np.sqrt(self.oifdata.t3[ii].t3amperr**2 + self.extra_error**2)
        self.uvcoords[0,:,ii] = self.oifdata.t3[ii].u1coord, self.oifdata.t3[ii].u2coord,\
                    -(self.oifdata.t3[ii].u1coord+self.oifdata.t3[ii].u2coord)
        self.uvcoords[1, :,ii] = self.oifdata.t3[ii].v1coord, self.oifdata.t3[ii].v2coord,\
                    -(self.oifdata.t3[ii].v1coord+self.oifdata.t3[ii].v2coord)

    for jj in range(self.nbl):
        self.v2[jj, :] = self.oifdata.vis2[jj].vis2data
        self.v2err[jj, :] = self.oifdata.vis2[jj].vis2err
        try:
            self.pha[jj, :] = self.oifdata.vis[jj].vispha
            self.phaerr[jj, :] = self.oifdata.vis[jj].visphaerr
            self.cv = np.sqrt(self.v2)*np.exp(-1j*self.pha)
        except:
            pass
        self.uvcoords_vis[0,jj] = self.oifdata.vis2[jj].ucoord
        self.uvcoords_vis[1,jj] = self.oifdata.vis2[jj].vcoord

    # hack right now to take care of 0 values, set to some limit, 0.001 right now
    floor = 0.00001
    self.cperr[self.cperr<floor] = self.cperr[self.cperr!=0.0].mean()
    self.phaerr[self.phaerr<floor] = self.phaerr[self.phaerr!=0.0].mean()
    self.v2err[self.v2err<floor] = self.v2err[self.v2err!=0.0].mean()
    
    # replicate the uv coordinates over the wavelength axis
    self.uvcoords = np.tile(self.uvcoords, (self.nwav, 1, 1, 1))
    self.uvcoords_vis = np.tile(self.uvcoords_vis, (self.nwav, 1, 1))
    # Now uvcoords is shape (nwav, 2, 3, ncps)
    # So we move nwav axis to the end:
    self.uvcoords = np.rollaxis(self.uvcoords, 0, 4)
    self.uvcoords_vis = np.rollaxis(self.uvcoords_vis, 0, 3)


def detec_calc_loop(dictlist):
    # ndetected should have shape (nsep, ncon, nang) -- the first 3 dimensions of the cp model
    simcps = dictlist['model'].copy()
    simcps += dictlist['randerrors']

    chi2bin_m_chi2null = np.sum( ((dictlist['model'] - simcps)**2 - (simcps**2)) /(dictlist['dataerrors']**2), axis=(-1,-2))

    detected = chi2bin_m_chi2null<0.0

    return detected

def detec_calc_loop_all(dictlist):
    # ndetected should have shape (nsep, ncon, nang) -- the first 3 dimensions of the cp model
    simcps = dictlist['model'].copy()
    simcps += dictlist['randerrors']

    nullmodel = np.zeros(simcps.shape)
    nullmodel[:,:,:,simcps.shape[3]/2:, :] = 1.0
    chi2bin_m_chi2null = np.sum( ((dictlist['model'] - simcps)**2 - ((simcps-nullmodel)**2)) /(dictlist['dataerrors']**2), axis=(-1,-2))

    detected = chi2bin_m_chi2null<0.0

    return detected


def logl(data, err, model):
    """
    Likelihood given data, errors, and the model values
    These are all shape (nobservable, nwav)
    """
    #for ii in range(len(model)):
    #   #ll += -0.5*np.log(2*np.pi)*data[2*ii].size + np.sum(-np.log(data[2*ii+1]**2)
    #return -0.5*np.log(2*np.pi) - np.sum(np.log(err)) - np.sum((model - data)**2/(2*data**2))
    #return -0.5*np.log(2*np.pi)*data.size + np.sum(-np.log(err**2) - 0.5*((model - data)/err)**2)
    chi2 = np.nansum(((data-model)/err)**2)
    loglike = -chi2/2
    #return np.sum(-np.log(err**2) - 0.5*((model - data)/err)**2)
    return loglike


def logl_cov(flatdata, invcovmat, flatmodel):
    """
    flatdata and flatmodel must be the same shape & len=1
    """
    v_i = (flatdata - flatmodel).reshape(flatdata.shape[0], 1)
    chi2 = np.dot(flatdata - flatmodel, np.dot(incovmat, (flatdata-flatmodel).T))
    loglike = -chi2/2
    return loglike

def reduced_chi2(data, err, model, dof=1.0):
    return (1/float(dof))*np.sum(((model - data)**2)/(err**2), axis=(-1,-2))


def assemble_cov_mat(self):
    meancps = np.mean(self.cp, axis=0)
    flat_mean_sub_cp = (self.cp - meancps[None,:]).flatten
    covmat = flat_mean_sub_cp[None,:]*flat_mean_sub_cp[:,None]
    return covmat

def chi2_grid_loop(args):
    # Model from data, err, uvcoords, params, wavls
    p0, p1, p2 = args['params']
    modelcps = np.rollaxis(model_cp_uv(args['uvcoords'], p0, p1, p2, 1/args['wavls']), 0, -1)
    chi2 = np.nansum( (modelcps - args['data'])**2 / args['error']**2, axis = (-1,-2))/ args["dof"]
    return chi2


def chi2_grid_loop_all(args):
    # Model from data, err, uvcoords, params, wavls
    p0, p1, p2 = args['params']
    modelcps = np.rollaxis(model_cp_uv(args['uvcoords'], p0, p1, p2, 1/args['wavls']), 0, -1)
    modelt3 = np.rollaxis(model_t3amp_uv(args['uvcoords'], p0, p1, p2, 1/args['wavls']), 0, -1)
    model = np.concatenate((modelcps, modelt3), axis=2)
    chi2 = np.nansum( (model - args['data'])**2 / args['error']**2, axis = (-1,-2))/ args["dof"]
    return chi2




