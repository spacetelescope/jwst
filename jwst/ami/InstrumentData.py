#
#  Module for defining data format, wavelength info, an mask geometry for these
#   instruments: GPI NRM, NIRISS AMI, VISIR SAM
#

import logging
import numpy as np
from astropy.io import fits
import os
import time

from .misctools.mask_definitions import NRM_mask_definitions
from .misctools import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

um = 1.0e-6   # need ?

# utility routines for InstrumentData classes


def set_cvsupport_threshold(instr, k, v):
    """
    Short Summary
    ------------
    Set threshold for where 'splodge' data in CV space contains signal

    Parameters
    ----------
    instr: InstrumentData
        InstrumentData instance

    k:

    v:

    Returns
    -------

    instr:
    thresh: Threshold for the absolute value of the FT(interferogram).
            Normalize abs(CV = FT(a)) for unity peak, and define the support
            of "good" CV when this is above threshold
    """

    instr.cvsupport_threshold[k] = v


class GPI:

    def __init__(self, reffile, **kwargs):
        """
        Short Summary
        ------------
        Initialize GPI class

        Parameters
        ----------
        reffile: <TYPE>
            one or a list of reference fits files for gpi-pipeline reduced
            containing useful header info

        < kwargs>

       ??  optionally: 'gpifilterpath'
                - Point to a directory which contains GPI filter files
                  code will read in the relevant file and pick out a
                  sample of wavelengths and transmissions that span the
                  list, to be used later to generate the model.
        """

        # only one NRM on GPI:
        self.arrname = "gpi_g10s40"
        self.pscale_mas = 14.1667  # 14.27 looks like a better match March 2019
        self.pscale_rad = utils.mas2rad(self.pscale_mas)
        self.mask = NRM_mask_definitions(maskname=self.arrname)
        self.mask.ctrs = np.array(self.mask.ctrs)
        # Hard code -1.5 deg rotation in data (April 2016)
        # (can be moved to NRM_mask_definitions later)
        self.mask.ctrs = utils.rotate2dccw(self.mask.ctrs, -3.7 * np.pi / 180.)
        # Add in hole/baseline properties ?
        self.holeshape = "circ"
        affine2d = kwargs.get('affine2d', None)
        if affine2d is None:
            self.affine2d = utils.Affine2d(mx=1.0, my=1.0,
                                           sx=0.0, sy=0.0,
                                           xo=0.0, yo=0.0, name="Ideal")
        else:
            self.affine2d = affine2d

        # Get info from reference file
        self.hdr0 = []
        self.hdr1 = []
        self.refdata = []
        if type(reffile) == str:
            reffile = [reffile, ]
        for fn in reffile:
            reffits = fits.open(fn)
            self.hdr0.append(reffits[0].header)
            self.hdr1.append(reffits[1].header)
            self.refdata.append(reffits[1].data)
            reffits.close()
        # instrument settings
        self.mode = self.hdr0[0]["DISPERSR"]
        self.obsmode = self.hdr0[0]["OBSMODE"]
        self.band = self.obsmode[-1]  # K1 is two letters
        self.ref_imgs_dir = "refimgs_" + self.band + "/"

        # finding centroid from phase slope only considered cv_phase data
        # when cv_abs data exceeds this cvsupport_threshold.
        # Absolute value of cv data normalized to unity maximum
        # for the threshold application.
        # Data reduction gurus: tweak the threshold value with experience...
        self.cvsupport_threshold = {"Y": 0.02, "J": 0.02, "H": 0.02, "1": 0.02,
                                    "2": 0.02}
        self.threshold = self.cvsupport_threshold[self.band]

        # Special mode for collapsed data
        if self.hdr1[0]["NAXIS3"] == 1:
            # This is just a way handle data that is manually collapsed.
            # Not a standard data format for GPI.
            self.mode = "WOLLASTON_FAKEOUT"

        # wavelength info: spect mode or pol more
        if "PRISM" in self.mode:
            # GPI's spectral mode
            self.nwav = self.hdr1[0]["NAXIS3"]
            self.wls = np.linspace(self.hdr1[0]["CRVAL3"],
                                   self.hdr1[0]["CRVAL3"]
                                   + self.hdr1[0]['CD3_3'] * self.nwav,
                                   self.nwav) * um
            self.eff_band = um * np.ones(self.nwav) * (self.wls[-1] - self.wls[0]) / self.nwav
        elif "WOLLASTON" in self.mode:
            # GPI's pol mode. Will define this for the DIFFERENTIAL
            # VISIBILITIES
            # diff vis: two channels 0/45 and 22/67

            self.nwav = 2

            # Define the bands in case we use a tophat filter
            band_ctrs = {"Y": (1.14+0.95)*um/2., "J": (1.35+1.12)*um/2.,
                         "H": (1.80+1.50)*um/2., "1": (2.19+1.9)*um/2.,
                         "2": (2.4+2.13)*um/2.0}
            band_wdth = {"Y": (1.14-0.95)*um, "J": (1.35-1.12)*um,
                         "H": (1.80-1.50)*um,
                         "1": (2.19-1.9)*um, "2": (2.4-2.13)*um}
            wghts = np.ones(15)
            wavls = np.linspace(band_ctrs[self.band]-band_wdth[self.band]/2.0,
                                band_ctrs[self.band]+band_wdth[self.band]/2.0,
                                num=15)

            if 'gpifilterpath' in kwargs:
                if self.band == "Y":
                    filterfile = kwargs["gpifilterpath"] + "GPI-filter-Y.fits"
                    cutoff = 0.7
                if self.band == "J":
                    filterfile = kwargs["gpifilterpath"] + "GPI-filter-J.fits"
                    cutoff = 0.7
                if self.band == "H":
                    filterfile = kwargs["gpifilterpath"] + "GPI-filter-H.fits"
                    cutoff = 0.7
                if self.band == "1":
                    filterfile = kwargs["gpifilterpath"] + "GPI-filter-K1.fits"
                    cutoff = 0.94
                if self.band == "2":
                    filterfile = kwargs["gpifilterpath"] + "GPI-filter-K2.fits"
                    cutoff = 0.94

                # Read in gpi filter file
                fitsfilter = fits.open(filterfile)[1].data
                wavls = []
                wghts = []
                # Sample the filter file so the filter is only 50 elements long
                skip = len(fitsfilter[0][0]) / 50
                for ii in range(len(fitsfilter[0][0])/skip):
                    if fitsfilter[0][1][skip*ii] > cutoff:
                        wavls.append(fitsfilter[0][0][skip*ii]*1.0e-6)
                        wghts.append(fitsfilter[0][1][skip*ii])

            lam_c = band_ctrs[self.band]
            lam_w = band_wdth[self.band]

            transmission = np.array([[wghts[f], wavls[f]]
                                    for f in range(len(wghts))])
            if "FAKEOUT" in self.mode:
                self.nwav = 1
                self.wls = [transmission, ]
                self.eff_band = np.array([lam_w, ])
            else:
                self.wls = [transmission, transmission]
                self.eff_band = np.array([lam_w, lam_w])
        else:
            log.critical("Check your reference file header."
                         + "Keywork DISPERSR='{0}' not understood".
                         format(self.mode))

        # For OIFits structure
        self.wavextension = (self.wls, self.eff_band)
        if "FAKEOUT" in self.mode:
            self.wavextension = ([lam_c, ], [lam_w, ])

        # Observation info
        self.telname = "GEMINI"
        self.ra, self.dec = self.hdr0[0]["RA"], self.hdr0[0]["DEC"]

        try:
            self.date = self.hdr0[0]["DATE"]
        except Exception:
            self.date = self.hdr0[0]["DATE-OBS"]

        self.month = self.date[-5:-3]
        self.day = self.date[-2:]
        self.year = self.date[:4]
        # self.parang = self.hdr0["PAR_ANG"]
        # AVPARANG added Aug 2 2016
        self.parangs = []
        self.itime = []
        self.crpa = []
        for ii in range(len(reffile)):
            self.parangs.append(self.hdr1[ii]["AVPARANG"])
            self.itime.append(self.hdr1[ii]["ITIME"])
            if "CRPA" in self.hdr0[ii]:
                self.crpa.append(self.hdr0[ii]["CRPA"])
        self.avparang = np.mean(self.parangs)
        if len(self.crpa) > 0:
            self.avcassang = np.mean(self.crpa)
        else:
            self.avcassang = 0.0
        self.parang_range = abs(self.parangs[-1] - self.parangs[0])
        self.totalinttime = np.sum(self.itime)
        try:
            self.pa = self.hdr0[0]["PA"]
        except Exception:
            self.pa = 0.0
        try:
            self.objname = self.hdr0[0]["OBJECT"]
        except Exception:
            self.objname = "Unknown"

        # Look for additional keyword arguments ?
    def read_data(self, fn):

        """
        Short Summary
        ------------
        ??

        Parameters
        ----------
        fn< ??:
           >>

        Returns
        -------
        ???

        """

        fitsfile = fits.open(fn)
        sci = fitsfile[1].data
        hdr = fitsfile[1].header
        fitsfile.close()

        if 'distorcorr' in fn:
            self.sub_dir_str = fn[-32:-21]
        else:
            self.sub_dir_str = fn[-21:-10]

        return sci, hdr


class VISIR:
    def __init__(self, objname="obj", band="11.3", src="A0V",
                 affine2d=None):
        """
        Short Summary
        ------------
        Initialize VISIR class

        Parameters
        ----------
        objname: string
            name of object observed

        band: ??
           ??

        src: ???
            if pysynphot is installed, can provide a guess at the
            stellar spectrum

        affine2d:  ??
            ???

        """
        self.band = band

        self.objname = objname

        self.arrname = "visir_sam"
        self.pscale_mas = 45
        self.pscale_rad = utils.mas2rad(self.pscale_mas)
        self.mask = NRM_mask_definitions(maskname=self.arrname)
        self.mask.ctrs = np.array(self.mask.ctrs)
        self.holeshape = "hex"

        if affine2d is None:
            self.affine2d = utils.Affine2d(mx=1.0, my=1.0,
                                           sx=0.0, sy=0.0,
                                           xo=0.0, yo=0.0, name="Ideal")
        else:
            self.affine2d = affine2d

        # tophat filter
        # this can be swapped with an actual filter file
        if self.band == "11.3":
            self.lam_c = 11.3*1e-6  # 11.3 microns
            self.lam_w = 0.6/11.3  # 0.6 micron bandpass
        elif self.band == "10.5":
            self.lam_c = 10.6*1e-6  # 11.3 microns
            self.lam_w = 0.1/10.5  # 0.6 micron bandpass
        else:
            raise ValueError("options for band are '11.3' or '10.5' \n{0} not supported".format(band))

        self.filt = utils.tophatfilter(self.lam_c, self.lam_w, npoints=10)
        try:
            self.wls = [utils.combine_transmission(self.filt, src), ]
        except Exception:
            self.wls = [self.filt, ]
        # self.wavextension = (self.lam_c, self.lam_w)
        # self.wavextension = (self.lam_c*np.ones(self.nexp), \
        #                     self.lam_w*np.ones(self.nexp))
        self.wavextension = ([self.lam_c, ], [self.lam_w, ])
        self.nwav = 1

        # finding centroid from phase slope only considered cv_phase data when
        # cv_abs data exceeds this.
        # absolute value of cv data normalized to unity maximum for the
        # threshold application.
        self.cvsupport_threshold = {"10.5": 0.02, "11.3": 0.02}
        self.threshold = self.cvsupport_threshold[self.band]

        self.ref_imgs_dir = "refimgs/"

        # Observation info - I don't know yet how JWST data headers will be
        # structured
        self.telname = "VLT"
        try:
            self.ra, self.dec = self.hdr0["RA"], self.hdr0["DEC"]
        except Exception:
            self.ra, self.dec = 00, 00
        try:
            self.date = self.hdr0["DATE"]
            self.month = self.date[-5:-3]
            self.day = self.date[-2:]
            self.year = self.date[:4]
        except Exception:
            lt = time.localtime()
            self.date = "{0}{1:02d}{2:02d}".format(lt[0], lt[1], lt[2])
            self.month = lt[1]
            self.day = lt[2]
            self.year = lt[0]
        try:
            self.parang = self.hdr0["PAR_ANG"]
        except Exception:
            self.parang = 00
        try:
            self.pa = self.hdr0["PA"]
        except Exception:
            self.pa = 00
        try:
            self.itime = self.hdr1["ITIME"]
        except Exception:
            self.itime = 00

    def read_data(self, fn):
        """
        Short Summary
        ------------
        ??

        Parameters
        ----------
        fn: ?
            ??

        Returns
        -------
        ??

        """

        # for datacube of exposures, need to read as 3D (nexp, npix, npix)
        fitsfile = fits.open(fn)
        scidata = fitsfile[0].data
        hdr = fitsfile[0].header
        # self.sub_dir_str = self.filt+"_"+objname
        self.sub_dir_str = '/' + fn.split('/')[-1].replace('.fits', '')
        # self.nexp = scidata.shape[0]
        # rewrite wavextension to be same length as nexp
        if len(scidata.shape) == 3:
            self.nwav = scidata.shape[0]
            [self.wls.append(self.wls[0]) for f in range(self.nwav-1)]
            return scidata, hdr
        elif len(scidata.shape) == 2:
            return np.array([scidata, ]), hdr
        else:
            log.critical("invalid data dimensions for NIRISS. Should have dimensionality of 2 or 3.")

        return scidata, hdr


class NIRISS:
    def __init__(self, filt, 
                       objname="obj", 
                       src="A0V", 
                       out_dir='', 
                       chooseholes=None, 
                       affine2d=None, 
                       bandpass=None,
                       **kwargs):
        """
        Short Summary
        ------------
        Initialize NIRISS class
        Either user has webbpsf and filter file can be read, or this will use a tophat and give a warning
        chooseholes: None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

        Parameters
        ----------
        filt: ?
 
        objname: ? 

        src: ?
 
        out_dir: ?

        chooseholes: ?

        affine2d: ?

        bandpass: None or [(wt,wlen),(wt,wlen),...].  Monochromatic would be e.g. [(1.0, 4.3e-6)]
                  Explicit bandpass arg will replace *all* niriss filter-specific variables with 
                  the given bandpass, so you can simulate 21cm psfs through something called "F430M"!
        """

        self.chooseholes = chooseholes
        self.objname = objname
        self.filt = filt

        self.lam_bin = {"F277W": 50, "F380M": 20, "F430M":40, "F480M":30}  # 12 waves in f430 - data analysis
                                                                           
        # use 150 for 3 waves ax f430m 
        self.lam_c = {"F277W": 2.77e-6,  # central wavelength (SI)
                      "F380M": 3.8e-6, 
                      "F430M": 4.28521033106325E-06,
                      "F480M": 4.8e-6}
        self.lam_w = {"F277W":0.2, "F380M": 0.1, "F430M": 0.0436, "F480M": 0.08} # fractional filter width 
        
        try:
            self.throughput = utils.trim_webbpsf_filter(self.filt, specbin=self.lam_bin[self.filt])
        except Exception:
            self.throughput = utils.tophatfilter(self.lam_c[self.filt], self.lam_w[self.filt], npoints=11)

        # Nominal
        """self.lam_c = {"F277W": 2.77e-6,  # central wavelength (SI)
                      "F380M": 3.8e-6, 
                      "F430M": 4.28521033106325E-06,
                      "F480M": 4.8e-6}
        self.lam_w = {"F277W":0.2, "F380M": 0.1, "F430M": 0.0436, "F480M": 0.08} # fractional filter width """

        # update nominal filter parameters with those of the filter read in and used in the analysis...
        # Weighted mean wavelength in meters, etc, etc "central wavelength" for the filter:
        from scipy.integrate import simps 
        self.lam_c[self.filt] = (self.throughput[:,1] * self.throughput[:,0]).sum() / self.throughput[:,0].sum() 
        area = simps(self.throughput[:,0], self.throughput[:,1])
        ew = area / self.throughput[:,0].max() # equivalent width
        beta = ew/self.lam_c[self.filt] # fractional bandpass
        self.lam_w[self.filt] = beta

        if bandpass is not None:
            bandpass = np.array(bandpass)  # type simplification
            wt = bandpass[:,0]
            wl = bandpass[:,1]
            cw = (wl*wt).sum()/wt.sum() # Weighted mean wavelength in meters "central wavelength"
            area = simps(wt, wl)
            ew = area / wt.max() # equivalent width
            beta = ew/cw # fractional bandpass
            self.lam_c = {"F277W":cw, "F380M": cw, "F430M": cw, "F480M": cw,}
            self.lam_w = {"F277W": beta, "F380M": beta, "F430M": beta, "F480M": beta} 
            self.throughput = bandpass

        try:
            self.wls = [utils.combine_transmission(self.throughput, src), ]
        except Exception:
            self.wls = [self.throughput, ]

        # Wavelength info for NIRISS bands F277W, F380M, F430M, or F480M
        self.wavextension = ([self.lam_c[self.filt],], [self.lam_w[self.filt],])
        self.nwav = 1

        # only one NRM on JWST:
        self.telname = "NIRISS"
        self.instrument = "NIRISS"
        self.arrname = "jwst_g7s6c"
        self.holeshape = "hex"
        self.mask = NRM_mask_definitions(maskname=self.arrname, chooseholes=chooseholes, holeshape=self.holeshape)

        # save affine deformation of pupil object or create a no-deformation object. 
        # We apply this when sampling the PSF, not to the pupil geometry.
        # This will set a default Ideal or a measured rotation, for example,
        # and include pixel scale changes due to pupil distortion.
        # Separating detector tilt pixel scale effects from pupil distortion effects is 
        # yet to be determined... see comments in Affine class definition.
        # AS AZG 2018 08 15 Ann Arbor
        if affine2d is None:
            self.affine2d = utils.Affine2d(mx=1.0,my=1.0, 
                                           sx=0.0,sy=0.0, 
                                           xo=0.0,yo=0.0, name="Ideal")
        else:
            self.affine2d = affine2d

        # finding centroid from phase slope only considered cv_phase data 
        # when cv_abs data exceeds this cvsupport_threshold.  
        # Absolute value of cv data normalized to unity maximum
        # for the threshold application.
        # Data reduction gurus: tweak the threshold value with experience...
        # Gurus: tweak cvsupport with use...
        self.cvsupport_threshold = {"F277W":0.02, "F380M": 0.02, "F430M": 0.02, "F480M": 0.02}

        self.threshold = self.cvsupport_threshold[filt]

        self.ref_imgs_dir = os.path.join(out_dir,"refimgs_"+self.filt+"/")

    def set_pscale(self, pscalex_deg=None, pscaley_deg=None):
        """
        Short Summary
        ------------
        Override pixel scale in header

        Parameters
        ----------
        pscalex_deg:

        pscaley_deg: 
       
     
        Returns
        -------
 
        """
        if pscalex_deg is not None:
            self.pscalex_deg = pscalex_deg
        if pscaley_deg is not None:
            self.pscaley_deg = pscaley_deg
        self.pscale_mas = 0.5 * (pscalex_deg + pscaley_deg) * (60*60*1000)
        self.pscale_rad = utils.mas2rad(self.pscale_mas)

    def read_data(self, fn, mode="slice"):
        """
        Short Summary
        ------------

        Parameters
        ----------
        fn:

        mode:
     
        Returns
        -------

        """
        # mode options are slice or UTR
        # for single slice data, need to read as 3D (1, npix, npix)
        # for utr data, need to read as 3D (ngroup, npix, npix)
        fitsfile = fits.open(fn)
        scidata = fitsfile[1].data
        prihdr = fitsfile[0].header
        scihdr = fitsfile[1].header
        self.updatewithheaderinfo(prihdr, scihdr) # mirage header or MAST header
        self.sub_dir_str = '/' + fn.split('/')[-1].replace('.fits', '')
        if len(scidata.shape) == 3:
            self.nwav = scidata.shape[0]
            [self.wls.append(self.wls[0]) for f in range(self.nwav-1)]
            return prihdr, scihdr, scidata
        elif len(scidata.shape) == 2:
            return prihdr, scihdr, scidata
        else:
            log.critical("invalid data dimensions for NIRISS. Should have dimensionality of 2 or 3.")

    def updatewithheaderinfo(self, ph, sh):
        """ 
        Short Summary
        ------------

        Parameters
        ----------
        ph: primary header

        sh: science header
     
        Returns
        -------

        """

        # The info4oif_dict will get pickled to disk when we write txt files of results.
        # That way we don't drag in objects like InstrumentData into code that reads text results
        # and writes oifits files - a simple built-in dictionary is the only object used in this transfer.
        info4oif_dict = {}
        self.telname = "JWST"; info4oif_dict['telname'] = self.telname
        info4oif_dict['telname'] = self.telname

        info4oif_dict['filt'] = self.filt
        info4oif_dict['lam_c'] = self.lam_c
        info4oif_dict['lam_w'] = self.lam_w
        info4oif_dict['lam_bin'] = self.lam_bin

        self.objname = ph["TARGNAME"]; info4oif_dict['objname'] = self.objname
        self.ra = ph["TARG_RA"]; info4oif_dict['ra'] = self.ra
        self.dec = ph["TARG_DEC"]; info4oif_dict['dec'] = self.dec

        # / axis 1 DS9 coordinate of the reference pixel (always POS1)
        # / axis 2 DS9 coordinate of the reference pixel (always POS1)
        self.crpix1 = sh["CRPIX1"]; info4oif_dict['crpix1'] = self.crpix1
        self.crpix2 = sh["CRPIX2"]; info4oif_dict['crpix2'] = self.crpix2
        # need Paul Goudfrooij's table for actual crval[1,2] for true pointing to detector pixel coords (DS9)

        self.instrument = ph["INSTRUME"]; info4oif_dict['instrument'] = self.instrument
        self.pupil = ph["PUPIL"]; info4oif_dict['pupil'] = self.pupil
        # "ImPlaneIA internal mask name" - same as maskname
        self.arrname = "jwst_g7s6c"; info4oif_dict['arrname'] = self.arrname   # ???unclean hardcoding

        # if data was generated on the average pixel scale of the header
        # then this is the right value that gets read in, and used in fringe fitting
        pscalex_deg = sh["CDELT1"]
        pscaley_deg = sh["CDELT2"]
        #
        # To use ami_sim's eg 65.6 mas/pixel scale we hardcode it here...
        pscalex_deg = 65.6 / (1000 * 60 * 60)
        pscaley_deg = 65.6 / (1000 * 60 * 60)
        info4oif_dict['pscalex_deg'] = pscalex_deg
        info4oif_dict['pscaley_deg'] = pscaley_deg
        # Whatever we did set is averaged for isotropic pixel scale here
        self.pscale_mas = 0.5 * (pscalex_deg + pscaley_deg) * (60*60*1000); info4oif_dict['pscale_mas'] = self.pscale_mas
        self.pscale_rad = utils.mas2rad(self.pscale_mas); info4oif_dict['pscale_rad'] = self.pscale_rad
        self.mask = NRM_mask_definitions(maskname=self.arrname, chooseholes=self.chooseholes,
                                         holeshape=self.holeshape)

        str = ph["DATE-OBS"]
        self.year = str[:4]; info4oif_dict['year'] = self.year
        self.month = str[5:7]; info4oif_dict['month'] = self.month
        self.day = str[8:10]; info4oif_dict['day'] = self.day
        self.parangh = sh["ROLL_REF"]; info4oif_dict['parangh'] = self.parangh
        self.pa = sh["PA_V3"]; info4oif_dict['pa'] = self.pa

        # An INTegration is NGROUPS "frames", not relevant here but context info.
        # 2d => "cal" file combines all INTegrations (ramps)
        # 3d=> "calints" file is a cube of all INTegrations (ramps)
        if sh["NAXIS"] == 2:
            # all INTegrations or 'ramps'
            self.itime = ph["EFFINTTM"] * ph["NINTS"]; info4oif_dict['itime'] = self.itime
        elif sh["NAXIS"] == 3:
            # each slice is one INTegration or 'ramp'
            self.itime = ph["EFFINTTM"]; info4oif_dict['itime'] = self.itime

        info4oif_dict['ctrs'] = self.mask.ctrs
        info4oif_dict['hdia'] = self.mask.hdia
        self.info4oif_dict = info4oif_dict # save it when writing extracted observables txt

    def reset_nwav(self, nwav):
        """
        Short Summary
        ------------
        rather than calling InstrumentData in the niriss example just to reset just call this routine

        Parameters
        ----------
        nwav:
     
        Returns
        -------
        """ 
        self.nwav = nwav


class NIRC2:
    def __init__(self, reffile, **kwargs):
        """
        Short Summary
        ------------
        Initialize NIRC2 class

        Parameters
        ----------
        reffile:     

        objname: string w/name of object observed

        src: if pysynphot is installed, can provide a guess at the stellar spectrum

        IFU simulation option set IFU = True
        """

        if "IFU" in kwargs:
            if kwargs["IFU"] is True:
                self.mode = "PRISM"
            else:
                self.mode = "BROADBAND"
        else:
            self.mode = "BROADBAND"
        if "src" in kwargs:
            src = kwargs["src"]
        else:
            pass

        affine2d = kwargs.get('affine2d', None)

        if affine2d is None:
            self.affine2d = utils.Affine2d(mx=1.0,my=1.0, 
                                           sx=0.0,sy=0.0, 
                                           xo=0.0,yo=0.0, name="Ideal")
        else:
            self.affine2d = affine2d

        self.arrname = "NIRC2_9NRM"
        self.pscale_mas = 9.952 # mas
        self.pscale_rad = utils.mas2rad(self.pscale_mas)
        self.mask = NRM_mask_definitions(maskname=self.arrname)
        self.mask.ctrs = np.array(self.mask.ctrs)
        # Hard code -1.5 deg rotation in data (April 2016)
        # (can be moved to NRM_mask_definitions later)
        self.mask.ctrs = utils.rotate2dccw(self.mask.ctrs, 1.0*np.pi/180.0)
        # Add in hole/baseline properties ?
        self.holeshape = "circ"

        self.threshold = 0.02

        # Get info from reference file 
        self.hdr = []

        if type(reffile) == str:
            reffile = [reffile, ]
        for fn in reffile:
            reffits = fits.open(fn)
            self.hdr.append(reffits[0].header)
            reffits.close()
        # instrument settings
        self.band = self.hdr[0]["FWINAME"]

        self.objname = self.hdr[0]["OBJECT"]

        # tophat filter
        # this can be swapped with an actual filter file
        # band_ctrs = {"Kp":1.633*um/2.0,"Lp":3.1*um}
        band_ctrs = {"J":1.248*um, "H":1.633*um, "CH4_short":1.5923*um,"Kp":2.2*um,"Lp":3.1*um}
        band_wdth = {"J":0.163*um, "H":0.296*um, "CH4_short":(0.1257)*um,"Kp":(0.3)*um, "Lp":(4.126 - 3.426)*um}

        lam_c = band_ctrs[self.band]
        lam_w = band_wdth[self.band]

        # wavelength info: spect mode or pol more
        if "PRISM" in self.mode:
            # GPI's spectral mode
            self.nwav = 36 # self.hdr1[0]["NAXIS3"]
            self.wls = np.linspace(lam_c - lam_w/2.0, lam_c+lam_w/2.0, num=36)*1e-6
            self.eff_band = um*np.ones(self.nwav)*(self.wls[-1] - self.wls[0])/self.nwav
            # For OIFits structure
            self.wavextension = (self.wls, self.eff_band)
        elif "BROADBAND" in self.mode:
            # Copied from GPI's pol mode.

            self.nwav = 1

            # Define the bands in case we use a tophat filter
            wghts = np.ones(11)
            wavls = np.linspace(lam_c-lam_w/2.0,
                                lam_c+lam_w/2.0, num=len(wghts))

            transmission = np.array([[wghts[f], wavls[f]] for f in range(len(wghts))])
            self.wls = [transmission, ]
            # self.wls = [np.sum(wghts*wavls) /float(len(wavls)), ]
            self.eff_band = np.array([lam_w, ])
            # For OIFits structure
            self.wavextension = ([lam_c,], [lam_w,])

        try:
            self.wls = [utils.combine_transmission(transmission, src), ]
        except Exception:
            self.wls = [transmission, ]

        self.nwav = 1

        # finding centroid from phase slope only considered cv_phase data when cv_abs data exceeds this.  
        # absolute value of cv data normalized to unity maximum for the threshold application.
        self.cvsupport_threshold = {"threshold":0.02} # Gurus: tweak with use...

        self.ref_imgs_dir = "refimgs/"

        # Observation info
        self.telname = "Keck"
        try:
            self.ra, self.dec = self.hdr[0]["RA"], self.hdr[0]["DEC"]
        except Exception:
            self.ra, self.dec = 00, 00
        try:
            self.date = self.hdr[0]["DATE"]
            self.month = self.date[8:10]
            self.day = self.date[5:7]
            self.year = self.date[:4]
        except Exception:
            lt = time.localtime()
            self.date = "{0}{1:02d}{2:02d}".format(lt[0],lt[1],lt[2])
            self.month = lt[1]
            self.day = lt[2]
            self.year = lt[0]

        self.parangs = []
        self.itime = []
        self.crpa = []
        self.rotposns = []
        self.instangs = []
        self.derotangs = []
        for ii in range(len(reffile)):
            self.parangs.append(self.hdr[ii]["PARANG"])
            self.rotposns.append(self.hdr[ii]["ROTPOSN"])
            self.instangs.append(self.hdr[ii]["INSTANGL"])
            # From Tom Esposito:
            # PARANG + ROTPOSN - INSTANGL - 0.262 
            self.derotangs.append(self.hdr[ii]["PARANG"] + self.hdr[ii]["ROTPOSN"]
                                  - self.hdr[ii]["INSTANGL"] - 0.262)
            self.itime.append(self.hdr[ii]["ITIME"])
            if "CRPA" in self.hdr[ii]:
                self.crpa.append(self.hdr[ii]["CRPA"])

        self.avderotang = np.mean(self.derotangs)
        self.avparang = np.mean(self.parangs)
        if len(self.crpa) > 0:
            self.avcassang = np.mean(self.crpa)
        else:
            self.avcassang = 0.0
        self.parang_range = abs(self.parangs[-1] - self.parangs[0])
        self.totalinttime = np.sum(self.itime)

        try:
            self.pa = self.hdr0["PA"]
        except Exception:
            self.pa = 00

        self.parang_range = abs(self.parangs[-1] - self.parangs[0])
        self.totalinttime = np.sum(self.itime)
        self.ref_imgs_dir = "refimgs/"

    def read_data(self, fn):
        """
        Short Summary
        ------------
        ??

        Parameters
        ----------
        fn: ?
            ??     

        Returns
        -------
        ??     

        """

        fitsfile = fits.open(fn)
        sci = fitsfile[0].data
        hdr = fitsfile[0].header
        fitsfile.close()
        # fitshdr = fitsfile[0].header
        self.sub_dir_str = fn.split("/")[-1][:-5]

        if len(sci.shape) == 3:
            if self.mode == "PRISM":
                return sci, hdr
            elif self.mode == "BROADBAND":
                self.nwav = sci.shape[0]
                [self.wls.append(self.wls[0]) for f in range(self.nwav-1)]
                return sci, hdr
        elif len(sci.shape) == 2:
            if self.mode == "BROADBAND":
                return np.array([sci,]), hdr
        else:
            log.critical("invalid data dimensions for NIRC2. Should have dimensionality of 2 or 3.")

        return sci, hdr
