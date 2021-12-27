""" Dictionary of basic instrument parameters
"""

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class InstrumentInfo():

    def __init__(self):
        """ Dictionary of basic instrument parameters

        These parameters are filled in from the cube reference file
        and  MIRI resolution reference file
        """

        self.multich_wavelength = []
        self.multich_sroi = []
        self.multich_wroi = []
        self.multich_power = []
        self.multich_softrad = []
        self.multich_scalerad = []

        self.prism_wavelength = []
        self.prism_sroi = []
        self.prism_wroi = []
        self.prism_power = []
        self.prism_softrad = []
        self.prism_scalerad = []

        self.med_wavelength = []
        self.med_sroi = []
        self.med_wroi = []
        self.med_power = []
        self.med_softrad = []
        self.med_scalerad = []

        self.high_wavelength = []
        self.high_sroi = []
        self.high_wroi = []
        self.high_power = []
        self.high_softrad = []
        self.high_scalerad = []
# _______________________________________________________________________
        # This is basic information on the MIRI channels
        # information that will not change:
        # number of slices, starting slice number,
        # ending # slice number and default scales
        self.Info = {}
# _____________________________________________________________________
# channel 1 parameters
        self.Info['1'] = {}
        self.Info['1']['nslices'] = 21
        self.Info['1']['start_slice'] = 101
        self.Info['1']['end_slice'] = 121
        self.Info['1']['xstart'] = 0
        self.Info['1']['xend'] = 512

        self.Info['1']['short'] = {}
        self.Info['1']['short']['ascale'] = None
        self.Info['1']['short']['bscale'] = None
        self.Info['1']['short']['wscale'] = None
        self.Info['1']['short']['sroi'] = None
        self.Info['1']['short']['wroi'] = None
        self.Info['1']['short']['wavemin'] = None
        self.Info['1']['short']['wavemax'] = None
        self.Info['1']['short']['softrad'] = None
        self.Info['1']['short']['scalerad'] = None
        self.Info['1']['short']['msm_power'] = None

        self.Info['1']['medium'] = {}
        self.Info['1']['medium']['ascale'] = None
        self.Info['1']['medium']['bscale'] = None
        self.Info['1']['medium']['wscale'] = None
        self.Info['1']['medium']['sroi'] = None
        self.Info['1']['medium']['wroi'] = None
        self.Info['1']['medium']['wavemin'] = None
        self.Info['1']['medium']['wavemax'] = None
        self.Info['1']['medium']['softrad'] = None
        self.Info['1']['medium']['scalerad'] = None
        self.Info['1']['medium']['msm_power'] = None

        self.Info['1']['long'] = {}
        self.Info['1']['long']['ascale'] = None
        self.Info['1']['long']['bscale'] = None
        self.Info['1']['long']['wscale'] = None
        self.Info['1']['long']['sroi'] = None
        self.Info['1']['long']['wroi'] = None
        self.Info['1']['long']['wavemin'] = None
        self.Info['1']['long']['wavemax'] = None
        self.Info['1']['long']['softrad'] = None
        self.Info['1']['long']['scalerad'] = None
        self.Info['1']['long']['msm_power'] = None
# _______________________________________________________________________
# channel 2 parameters
        self.Info['2'] = {}
        self.Info['2']['nslices'] = 17
        self.Info['2']['start_slice'] = 201
        self.Info['2']['end_slice'] = 217
        self.Info['2']['xstart'] = 513
        self.Info['2']['xend'] = 1031

        self.Info['2']['short'] = {}
        self.Info['2']['short']['ascale'] = None
        self.Info['2']['short']['bscale'] = None
        self.Info['2']['short']['wscale'] = None
        self.Info['2']['short']['sroi'] = None
        self.Info['2']['short']['wroi'] = None
        self.Info['2']['short']['wavemin'] = None
        self.Info['2']['short']['wavemax'] = None
        self.Info['2']['short']['softrad'] = None
        self.Info['2']['short']['scalerad'] = None
        self.Info['2']['short']['msm_power'] = None

        self.Info['2']['medium'] = {}
        self.Info['2']['medium']['ascale'] = None
        self.Info['2']['medium']['bscale'] = None
        self.Info['2']['medium']['wscale'] = None
        self.Info['2']['medium']['sroi'] = None
        self.Info['2']['medium']['wroi'] = None
        self.Info['2']['medium']['wavemin'] = None
        self.Info['2']['medium']['wavemax'] = None
        self.Info['2']['medium']['softrad'] = None
        self.Info['2']['medium']['scalerad'] = None
        self.Info['2']['medium']['msm_power'] = None

        self.Info['2']['long'] = {}
        self.Info['2']['long']['ascale'] = None
        self.Info['2']['long']['bscale'] = None
        self.Info['2']['long']['wscale'] = None
        self.Info['2']['long']['sroi'] = None
        self.Info['2']['long']['wroi'] = None
        self.Info['2']['long']['wavemin'] = None
        self.Info['2']['long']['wavemax'] = None
        self.Info['2']['long']['softrad'] = None
        self.Info['2']['long']['scalerad'] = None
        self.Info['2']['long']['msm_power'] = None
# ________________________________________________________________
# channel 3 parameters
        self.Info['3'] = {}
        self.Info['3']['nslices'] = 16
        self.Info['3']['start_slice'] = 301
        self.Info['3']['end_slice'] = 316
        self.Info['3']['xstart'] = 513
        self.Info['3']['xend'] = 1031

        self.Info['3']['short'] = {}
        self.Info['3']['short']['ascale'] = None
        self.Info['3']['short']['bscale'] = None
        self.Info['3']['short']['wscale'] = None
        self.Info['3']['short']['sroi'] = None
        self.Info['3']['short']['wroi'] = None
        self.Info['3']['short']['wavemin'] = None
        self.Info['3']['short']['wavemax'] = None
        self.Info['3']['short']['softrad'] = None
        self.Info['3']['short']['scalerad'] = None
        self.Info['3']['short']['msm_power'] = None

        self.Info['3']['medium'] = {}
        self.Info['3']['medium']['ascale'] = None
        self.Info['3']['medium']['bscale'] = None
        self.Info['3']['medium']['wscale'] = None
        self.Info['3']['medium']['sroi'] = None
        self.Info['3']['medium']['wroi'] = None
        self.Info['3']['medium']['wavemin'] = None
        self.Info['3']['medium']['wavemax'] = None
        self.Info['3']['medium']['softrad'] = None
        self.Info['3']['medium']['scalerad'] = None
        self.Info['3']['medium']['msm_power'] = None

        self.Info['3']['long'] = {}
        self.Info['3']['long']['ascale'] = None
        self.Info['3']['long']['bscale'] = None
        self.Info['3']['long']['wscale'] = None
        self.Info['3']['long']['sroi'] = None
        self.Info['3']['long']['wroi'] = None
        self.Info['3']['long']['wavemin'] = None
        self.Info['3']['long']['wavemax'] = None
        self.Info['3']['long']['softrad'] = None
        self.Info['3']['long']['scalerad'] = None
        self.Info['3']['long']['msm_power'] = None
# _________________________________________________________________
# channel 4 parameters
        self.Info['4'] = {}
        self.Info['4']['nslices'] = 12
        self.Info['4']['start_slice'] = 401
        self.Info['4']['end_slice'] = 412
        self.Info['4']['xstart'] = 0
        self.Info['4']['xend'] = 512

        self.Info['4']['short'] = {}
        self.Info['4']['short']['ascale'] = None
        self.Info['4']['short']['bscale'] = None
        self.Info['4']['short']['wscale'] = None
        self.Info['4']['short']['sroi'] = None
        self.Info['4']['short']['wroi'] = None
        self.Info['4']['short']['wavemin'] = None
        self.Info['4']['short']['wavemax'] = None
        self.Info['4']['short']['softrad'] = None
        self.Info['4']['short']['scalerad'] = None
        self.Info['4']['short']['msm_power'] = None

        self.Info['4']['medium'] = {}
        self.Info['4']['medium']['ascale'] = None
        self.Info['4']['medium']['bscale'] = None
        self.Info['4']['medium']['wscale'] = None
        self.Info['4']['medium']['sroi'] = None
        self.Info['4']['medium']['wroi'] = None
        self.Info['4']['medium']['wavemin'] = None
        self.Info['4']['medium']['wavemax'] = None
        self.Info['4']['medium']['softrad'] = None
        self.Info['4']['medium']['scalerad'] = None
        self.Info['4']['medium']['msm_power'] = None

        self.Info['4']['long'] = {}
        self.Info['4']['long']['ascale'] = None
        self.Info['4']['long']['bscale'] = None
        self.Info['4']['long']['wscale'] = None
        self.Info['4']['long']['wroi'] = None
        self.Info['4']['long']['sroi'] = None
        self.Info['4']['long']['wavemin'] = None
        self.Info['4']['long']['wavemax'] = None
        self.Info['4']['long']['softrad'] = None
        self.Info['4']['long']['scalerad'] = None
        self.Info['4']['long']['msm_power'] = None

        self.Info['1']['short-medium'] = {}
        self.Info['1']['short-medium']['ascale'] = None
        self.Info['1']['short-medium']['bscale'] = None
        self.Info['1']['short-medium']['wscale'] = None
        self.Info['1']['short-medium']['sroi'] = None
        self.Info['1']['short-medium']['wroi'] = None
        self.Info['1']['short-medium']['wavemin'] = None
        self.Info['1']['short-medium']['wavemax'] = None
        self.Info['1']['short-medium']['softrad'] = None
        self.Info['1']['short-medium']['scalerad'] = None
        self.Info['1']['short-medium']['msm_power'] = None

        self.Info['2']['short-medium'] = {}
        self.Info['2']['short-medium']['ascale'] = None
        self.Info['2']['short-medium']['bscale'] = None
        self.Info['2']['short-medium']['wscale'] = None
        self.Info['2']['short-medium']['sroi'] = None
        self.Info['2']['short-medium']['wroi'] = None
        self.Info['2']['short-medium']['wavemin'] = None
        self.Info['2']['short-medium']['wavemax'] = None
        self.Info['2']['short-medium']['softrad'] = None
        self.Info['2']['short-medium']['scalerad'] = None
        self.Info['2']['short-medium']['msm_power'] = None

        self.Info['3']['short-medium'] = {}
        self.Info['3']['short-medium']['ascale'] = None
        self.Info['3']['short-medium']['bscale'] = None
        self.Info['3']['short-medium']['wscale'] = None
        self.Info['3']['short-medium']['sroi'] = None
        self.Info['3']['short-medium']['wroi'] = None
        self.Info['3']['short-medium']['wavemin'] = None
        self.Info['3']['short-medium']['wavemax'] = None
        self.Info['3']['short-medium']['softrad'] = None
        self.Info['3']['short-medium']['scalerad'] = None
        self.Info['3']['short-medium']['msm_power'] = None

        self.Info['4']['short-medium'] = {}
        self.Info['4']['short-medium']['ascale'] = None
        self.Info['4']['short-medium']['bscale'] = None
        self.Info['4']['short-medium']['wscale'] = None
        self.Info['4']['short-medium']['sroi'] = None
        self.Info['4']['short-medium']['wroi'] = None
        self.Info['4']['short-medium']['wavemin'] = None
        self.Info['4']['short-medium']['wavemax'] = None
        self.Info['4']['short-medium']['softrad'] = None
        self.Info['4']['short-medium']['scalerad'] = None
        self.Info['4']['short-medium']['msm_power'] = None

        self.Info['1']['short-long'] = {}
        self.Info['1']['short-long']['ascale'] = None
        self.Info['1']['short-long']['bscale'] = None
        self.Info['1']['short-long']['wscale'] = None
        self.Info['1']['short-long']['sroi'] = None
        self.Info['1']['short-long']['wroi'] = None
        self.Info['1']['short-long']['wavemin'] = None
        self.Info['1']['short-long']['wavemax'] = None
        self.Info['1']['short-long']['softrad'] = None
        self.Info['1']['short-long']['scalerad'] = None
        self.Info['1']['short-long']['msm_power'] = None

        self.Info['2']['short-long'] = {}
        self.Info['2']['short-long']['ascale'] = None
        self.Info['2']['short-long']['bscale'] = None
        self.Info['2']['short-long']['wscale'] = None
        self.Info['2']['short-long']['sroi'] = None
        self.Info['2']['short-long']['wroi'] = None
        self.Info['2']['short-long']['wavemin'] = None
        self.Info['2']['short-long']['wavemax'] = None
        self.Info['2']['short-long']['softrad'] = None
        self.Info['2']['short-long']['scalerad'] = None
        self.Info['2']['short-long']['msm_power'] = None

        self.Info['3']['short-long'] = {}
        self.Info['3']['short-long']['ascale'] = None
        self.Info['3']['short-long']['bscale'] = None
        self.Info['3']['short-long']['wscale'] = None
        self.Info['3']['short-long']['sroi'] = None
        self.Info['3']['short-long']['wroi'] = None
        self.Info['3']['short-long']['wavemin'] = None
        self.Info['3']['short-long']['wavemax'] = None
        self.Info['3']['short-long']['softrad'] = None
        self.Info['3']['short-long']['scalerad'] = None
        self.Info['3']['short-long']['msm_power'] = None

        self.Info['4']['short-long'] = {}
        self.Info['4']['short-long']['ascale'] = None
        self.Info['4']['short-long']['bscale'] = None
        self.Info['4']['short-long']['wscale'] = None
        self.Info['4']['short-long']['sroi'] = None
        self.Info['4']['short-long']['wroi'] = None
        self.Info['4']['short-long']['wavemin'] = None
        self.Info['4']['short-long']['wavemax'] = None
        self.Info['4']['short-long']['softrad'] = None
        self.Info['4']['short-long']['scalerad'] = None
        self.Info['4']['short-long']['msm_power'] = None

        self.Info['1']['medium-short'] = {}
        self.Info['1']['medium-short']['ascale'] = None
        self.Info['1']['medium-short']['bscale'] = None
        self.Info['1']['medium-short']['wscale'] = None
        self.Info['1']['medium-short']['sroi'] = None
        self.Info['1']['medium-short']['wroi'] = None
        self.Info['1']['medium-short']['wavemin'] = None
        self.Info['1']['medium-short']['wavemax'] = None
        self.Info['1']['medium-short']['softrad'] = None
        self.Info['1']['medium-short']['scalerad'] = None
        self.Info['1']['medium-short']['msm_power'] = None

        self.Info['2']['medium-short'] = {}
        self.Info['2']['medium-short']['ascale'] = None
        self.Info['2']['medium-short']['bscale'] = None
        self.Info['2']['medium-short']['wscale'] = None
        self.Info['2']['medium-short']['sroi'] = None
        self.Info['2']['medium-short']['wroi'] = None
        self.Info['2']['medium-short']['wavemin'] = None
        self.Info['2']['medium-short']['wavemax'] = None
        self.Info['2']['medium-short']['softrad'] = None
        self.Info['2']['medium-short']['scalerad'] = None
        self.Info['2']['medium-short']['msm_power'] = None

        self.Info['3']['medium-short'] = {}
        self.Info['3']['medium-short']['ascale'] = None
        self.Info['3']['medium-short']['bscale'] = None
        self.Info['3']['medium-short']['wscale'] = None
        self.Info['3']['medium-short']['sroi'] = None
        self.Info['3']['medium-short']['wroi'] = None
        self.Info['3']['medium-short']['wavemin'] = None
        self.Info['3']['medium-short']['wavemax'] = None
        self.Info['3']['medium-short']['softrad'] = None
        self.Info['3']['medium-short']['scalerad'] = None
        self.Info['3']['medium-short']['msm_power'] = None

        self.Info['4']['medium-short'] = {}
        self.Info['4']['medium-short']['ascale'] = None
        self.Info['4']['medium-short']['bscale'] = None
        self.Info['4']['medium-short']['wscale'] = None
        self.Info['4']['medium-short']['sroi'] = None
        self.Info['4']['medium-short']['wroi'] = None
        self.Info['4']['medium-short']['wavemin'] = None
        self.Info['4']['medium-short']['wavemax'] = None
        self.Info['4']['medium-short']['softrad'] = None
        self.Info['4']['medium-short']['scalerad'] = None
        self.Info['4']['medium-short']['msm_power'] = None

        self.Info['1']['medium-long'] = {}
        self.Info['1']['medium-long']['ascale'] = None
        self.Info['1']['medium-long']['bscale'] = None
        self.Info['1']['medium-long']['wscale'] = None
        self.Info['1']['medium-long']['sroi'] = None
        self.Info['1']['medium-long']['wroi'] = None
        self.Info['1']['medium-long']['wavemin'] = None
        self.Info['1']['medium-long']['wavemax'] = None
        self.Info['1']['medium-long']['softrad'] = None
        self.Info['1']['medium-long']['scalerad'] = None
        self.Info['1']['medium-long']['msm_power'] = None

        self.Info['2']['medium-long'] = {}
        self.Info['2']['medium-long']['ascale'] = None
        self.Info['2']['medium-long']['bscale'] = None
        self.Info['2']['medium-long']['wscale'] = None
        self.Info['2']['medium-long']['sroi'] = None
        self.Info['2']['medium-long']['wroi'] = None
        self.Info['2']['medium-long']['wavemin'] = None
        self.Info['2']['medium-long']['wavemax'] = None
        self.Info['2']['medium-long']['softrad'] = None
        self.Info['2']['medium-long']['scalerad'] = None
        self.Info['2']['medium-long']['msm_power'] = None
        self.Info['2']['medium-long']['rp_wave_cutoff'] = None
        self.Info['2']['medium-long']['rp_a_low'] = None
        self.Info['2']['medium-long']['rp_b_low'] = None
        self.Info['2']['medium-long']['rp_c_low'] = None
        self.Info['2']['medium-long']['rp_a_high'] = None
        self.Info['2']['medium-long']['rp_b_high'] = None
        self.Info['2']['medium-long']['rp_c_high'] = None
        self.Info['2']['medium-long']['rp_a_ave'] = None
        self.Info['2']['medium-long']['rp_b_ave'] = None
        self.Info['2']['medium-long']['rp_c_ave'] = None

        self.Info['3']['medium-long'] = {}
        self.Info['3']['medium-long']['ascale'] = None
        self.Info['3']['medium-long']['bscale'] = None
        self.Info['3']['medium-long']['wscale'] = None
        self.Info['3']['medium-long']['sroi'] = None
        self.Info['3']['medium-long']['wroi'] = None
        self.Info['3']['medium-long']['wavemin'] = None
        self.Info['3']['medium-long']['wavemax'] = None
        self.Info['3']['medium-long']['softrad'] = None
        self.Info['3']['medium-long']['scalerad'] = None
        self.Info['3']['medium-long']['msm_power'] = None

        self.Info['4']['medium-long'] = {}
        self.Info['4']['medium-long']['ascale'] = None
        self.Info['4']['medium-long']['bscale'] = None
        self.Info['4']['medium-long']['wscale'] = None
        self.Info['4']['medium-long']['sroi'] = None
        self.Info['4']['medium-long']['wroi'] = None
        self.Info['4']['medium-long']['wavemin'] = None
        self.Info['4']['medium-long']['wavemax'] = None
        self.Info['4']['medium-long']['softrad'] = None
        self.Info['4']['medium-long']['scalerad'] = None
        self.Info['4']['medium-long']['msm_power'] = None

        self.Info['1']['long-short'] = {}
        self.Info['1']['long-short']['ascale'] = None
        self.Info['1']['long-short']['bscale'] = None
        self.Info['1']['long-short']['wscale'] = None
        self.Info['1']['long-short']['sroi'] = None
        self.Info['1']['long-short']['wroi'] = None
        self.Info['1']['long-short']['wavemin'] = None
        self.Info['1']['long-short']['wavemax'] = None
        self.Info['1']['long-short']['softrad'] = None
        self.Info['1']['long-short']['scalerad'] = None
        self.Info['1']['long-short']['msm_power'] = None

        self.Info['2']['long-short'] = {}
        self.Info['2']['long-short']['ascale'] = None
        self.Info['2']['long-short']['bscale'] = None
        self.Info['2']['long-short']['wscale'] = None
        self.Info['2']['long-short']['sroi'] = None
        self.Info['2']['long-short']['wroi'] = None
        self.Info['2']['long-short']['wavemin'] = None
        self.Info['2']['long-short']['wavemax'] = None
        self.Info['2']['long-short']['softrad'] = None
        self.Info['2']['long-short']['scalerad'] = None
        self.Info['2']['long-short']['msm_power'] = None

        self.Info['3']['long-short'] = {}
        self.Info['3']['long-short']['ascale'] = None
        self.Info['3']['long-short']['bscale'] = None
        self.Info['3']['long-short']['wscale'] = None
        self.Info['3']['long-short']['sroi'] = None
        self.Info['3']['long-short']['wroi'] = None
        self.Info['3']['long-short']['wavemin'] = None
        self.Info['3']['long-short']['wavemax'] = None
        self.Info['3']['long-short']['softrad'] = None
        self.Info['3']['long-short']['scalerad'] = None
        self.Info['3']['long-short']['msm_power'] = None

        self.Info['4']['long-short'] = {}
        self.Info['4']['long-short']['ascale'] = None
        self.Info['4']['long-short']['bscale'] = None
        self.Info['4']['long-short']['wscale'] = None
        self.Info['4']['long-short']['sroi'] = None
        self.Info['4']['long-short']['wroi'] = None
        self.Info['4']['long-short']['wavemin'] = None
        self.Info['4']['long-short']['wavemax'] = None
        self.Info['4']['long-short']['softrad'] = None
        self.Info['4']['long-short']['scalerad'] = None
        self.Info['4']['long-short']['msm_power'] = None

        self.Info['1']['long-medium'] = {}
        self.Info['1']['long-medium']['ascale'] = None
        self.Info['1']['long-medium']['bscale'] = None
        self.Info['1']['long-medium']['wscale'] = None
        self.Info['1']['long-medium']['sroi'] = None
        self.Info['1']['long-medium']['wroi'] = None
        self.Info['1']['long-medium']['wavemin'] = None
        self.Info['1']['long-medium']['wavemax'] = None
        self.Info['1']['long-medium']['softrad'] = None
        self.Info['1']['long-medium']['scalerad'] = None
        self.Info['1']['long-medium']['msm_power'] = None

        self.Info['2']['long-medium'] = {}
        self.Info['2']['long-medium']['ascale'] = None
        self.Info['2']['long-medium']['bscale'] = None
        self.Info['2']['long-medium']['wscale'] = None
        self.Info['2']['long-medium']['sroi'] = None
        self.Info['2']['long-medium']['wroi'] = None
        self.Info['2']['long-medium']['wavemin'] = None
        self.Info['2']['long-medium']['wavemax'] = None
        self.Info['2']['long-medium']['softrad'] = None
        self.Info['2']['long-medium']['scalerad'] = None
        self.Info['2']['long-medium']['msm_power'] = None

        self.Info['3']['long-medium'] = {}
        self.Info['3']['long-medium']['ascale'] = None
        self.Info['3']['long-medium']['bscale'] = None
        self.Info['3']['long-medium']['wscale'] = None
        self.Info['3']['long-medium']['sroi'] = None
        self.Info['3']['long-medium']['wroi'] = None
        self.Info['3']['long-medium']['wavemin'] = None
        self.Info['3']['long-medium']['wavemax'] = None
        self.Info['3']['long-medium']['softrad'] = None
        self.Info['3']['long-medium']['scalerad'] = None
        self.Info['3']['long-medium']['msm_power'] = None

        self.Info['4']['long-medium'] = {}
        self.Info['4']['long-medium']['ascale'] = None
        self.Info['4']['long-medium']['bscale'] = None
        self.Info['4']['long-medium']['wscale'] = None
        self.Info['4']['long-medium']['sroi'] = None
        self.Info['4']['long-medium']['wroi'] = None
        self.Info['4']['long-medium']['wavemin'] = None
        self.Info['4']['long-medium']['wavemax'] = None
        self.Info['4']['long-medium']['softrad'] = None
        self.Info['4']['long-medium']['scalerad'] = None
        self.Info['4']['long-medium']['msm_power'] = None

# ####################################################################
# NIRSPEC Parameters
        self.Info['prism'] = {}
        self.Info['prism']['clear'] = {}
        self.Info['prism']['clear']['nslices'] = 30
        self.Info['prism']['clear']['wscale'] = 0.005
        self.Info['prism']['clear']['ascale'] = 0.1
        self.Info['prism']['clear']['bscale'] = 0.1
        self.Info['prism']['clear']['wroi'] = None
        self.Info['prism']['clear']['sroi'] = None
        self.Info['prism']['clear']['wavemin'] = None
        self.Info['prism']['clear']['wavemax'] = None
        self.Info['prism']['clear']['softrad'] = None
        self.Info['prism']['clear']['msm_power'] = None
        self.Info['prism']['clear']['scalerad'] = None

        self.Info['prism']['opaque'] = {}
        self.Info['prism']['opaque']['nslices'] = 30
        self.Info['prism']['opaque']['wscale'] = 0.005
        self.Info['prism']['opaque']['ascale'] = 0.1
        self.Info['prism']['opaque']['bscale'] = 0.1

        self.Info['g140m'] = {}
        self.Info['g140m']['f070lp'] = {}
        self.Info['g140m']['f070lp']['nslices'] = 30
        self.Info['g140m']['f070lp']['wscale'] = 0.000636
        self.Info['g140m']['f070lp']['ascale'] = 0.1
        self.Info['g140m']['f070lp']['bscale'] = 0.1
        self.Info['g140m']['f070lp']['wroi'] = None
        self.Info['g140m']['f070lp']['sroi'] = None
        self.Info['g140m']['f070lp']['wavemin'] = None
        self.Info['g140m']['f070lp']['wavemax'] = None
        self.Info['g140m']['f070lp']['softrad'] = None
        self.Info['g140m']['f070lp']['msm_power'] = None
        self.Info['g140m']['f070lp']['scalerad'] = None

        self.Info['g140m']['f100lp'] = {}
        self.Info['g140m']['f100lp']['nslices'] = 30
        self.Info['g140m']['f100lp']['wscale'] = 0.000636
        self.Info['g140m']['f100lp']['ascale'] = 0.1
        self.Info['g140m']['f100lp']['bscale'] = 0.1
        self.Info['g140m']['f100lp']['wroi'] = None
        self.Info['g140m']['f100lp']['sroi'] = None
        self.Info['g140m']['f100lp']['wavemin'] = None
        self.Info['g140m']['f100lp']['wavemax'] = None
        self.Info['g140m']['f100lp']['softrad'] = None
        self.Info['g140m']['f100lp']['msm_power'] = None
        self.Info['g140m']['f100lp']['scalerad'] = None

        self.Info['g140m']['opaque'] = {}
        self.Info['g140m']['opaque']['nslices'] = 30
        self.Info['g140m']['opaque']['wscale'] = 0.000636
        self.Info['g140m']['opaque']['ascale'] = 0.1
        self.Info['g140m']['opaque']['bscale'] = 0.1

        self.Info['g235m'] = {}
        self.Info['g235m']['f170lp'] = {}
        self.Info['g235m']['f170lp']['nslices'] = 30
        self.Info['g235m']['f170lp']['wscale'] = 0.00106
        self.Info['g235m']['f170lp']['ascale'] = 0.1
        self.Info['g235m']['f170lp']['bscale'] = 0.1
        self.Info['g235m']['f170lp']['wroi'] = None
        self.Info['g235m']['f170lp']['sroi'] = None
        self.Info['g235m']['f170lp']['wavemin'] = None
        self.Info['g235m']['f170lp']['wavemax'] = None
        self.Info['g235m']['f170lp']['softrad'] = None
        self.Info['g235m']['f170lp']['msm_power'] = None
        self.Info['g235m']['f170lp']['scalerad'] = None

        self.Info['g235m']['opaque'] = {}
        self.Info['g235m']['opaque']['nslices'] = 30
        self.Info['g235m']['opaque']['wscale'] = 0.00106
        self.Info['g235m']['opaque']['ascale'] = 0.1
        self.Info['g235m']['opaque']['bscale'] = 0.1

        self.Info['g395m'] = {}
        self.Info['g395m']['f290lp'] = {}
        self.Info['g395m']['f290lp']['nslices'] = 30
        self.Info['g395m']['f290lp']['wscale'] = 0.00179
        self.Info['g395m']['f290lp']['ascale'] = 0.1
        self.Info['g395m']['f290lp']['bscale'] = 0.1
        self.Info['g395m']['f290lp']['wroi'] = None
        self.Info['g395m']['f290lp']['sroi'] = None
        self.Info['g395m']['f290lp']['wavemin'] = None
        self.Info['g395m']['f290lp']['wavemax'] = None
        self.Info['g395m']['f290lp']['softrad'] = None
        self.Info['g395m']['f290lp']['msm_power'] = None
        self.Info['g395m']['f290lp']['scalerad'] = None

        self.Info['g395m']['opaque'] = {}
        self.Info['g395m']['opaque']['nslices'] = 30
        self.Info['g395m']['opaque']['wscale'] = 0.00179
        self.Info['g395m']['opaque']['ascale'] = 0.1
        self.Info['g395m']['opaque']['bscale'] = 0.1

        self.Info['g140h'] = {}
        self.Info['g140h']['f070lp'] = {}
        self.Info['g140h']['f070lp']['nslices'] = 30
        self.Info['g140h']['f070lp']['wscale'] = 0.000235
        self.Info['g140h']['f070lp']['ascale'] = 0.1
        self.Info['g140h']['f070lp']['bscale'] = 0.1
        self.Info['g140h']['f070lp']['wroi'] = None
        self.Info['g140h']['f070lp']['sroi'] = None
        self.Info['g140h']['f070lp']['wavemin'] = None
        self.Info['g140h']['f070lp']['wavemax'] = None
        self.Info['g140h']['f070lp']['softrad'] = None
        self.Info['g140h']['f070lp']['msm_power'] = None
        self.Info['g140h']['f070lp']['scalerad'] = None

        self.Info['g140h']['f100lp'] = {}
        self.Info['g140h']['f100lp']['nslices'] = 30
        self.Info['g140h']['f100lp']['wscale'] = 0.000235
        self.Info['g140h']['f100lp']['ascale'] = 0.1
        self.Info['g140h']['f100lp']['bscale'] = 0.1
        self.Info['g140h']['f100lp']['wroi'] = None
        self.Info['g140h']['f100lp']['sroi'] = None
        self.Info['g140h']['f100lp']['wavemin'] = None
        self.Info['g140h']['f100lp']['wavemax'] = None
        self.Info['g140h']['f100lp']['softrad'] = None
        self.Info['g140h']['f100lp']['msm_power'] = None
        self.Info['g140h']['f100lp']['scalerad'] = None

        self.Info['g140h']['opaque'] = {}
        self.Info['g140h']['opaque']['nslices'] = 30
        self.Info['g140h']['opaque']['wscale'] = 0.000235
        self.Info['g140h']['opaque']['ascale'] = 0.1
        self.Info['g140h']['opaque']['bscale'] = 0.1

        self.Info['g235h'] = {}
        self.Info['g235h']['f170lp'] = {}
        self.Info['g235h']['f170lp']['nslices'] = 30
        self.Info['g235h']['f170lp']['wscale'] = 0.000396
        self.Info['g235h']['f170lp']['ascale'] = 0.1
        self.Info['g235h']['f170lp']['bscale'] = 0.1
        self.Info['g235h']['f170lp']['wroi'] = None
        self.Info['g235h']['f170lp']['sroi'] = None
        self.Info['g235h']['f170lp']['wavemin'] = None
        self.Info['g235h']['f170lp']['wavemax'] = None
        self.Info['g235h']['f170lp']['softrad'] = None
        self.Info['g235h']['f170lp']['msm_power'] = None
        self.Info['g235h']['f170lp']['scalerad'] = None

        self.Info['g235h']['opaque'] = {}
        self.Info['g235h']['opaque']['nslices'] = 30
        self.Info['g235h']['opaque']['wscale'] = 0.000396
        self.Info['g235h']['opaque']['ascale'] = 0.1
        self.Info['g235h']['opaque']['bscale'] = 0.1

        self.Info['g395h'] = {}
        self.Info['g395h']['f290lp'] = {}
        self.Info['g395h']['f290lp']['nslices'] = 30
        self.Info['g395h']['f290lp']['wscale'] = 0.000665
        self.Info['g395h']['f290lp']['ascale'] = 0.1
        self.Info['g395h']['f290lp']['bscale'] = 0.1
        self.Info['g395h']['f290lp']['wroi'] = None
        self.Info['g395h']['f290lp']['sroi'] = None
        self.Info['g395h']['f290lp']['wavemin'] = None
        self.Info['g395h']['f290lp']['wavemax'] = None
        self.Info['g395h']['f290lp']['softrad'] = None
        self.Info['g395h']['f290lp']['msm_power'] = None
        self.Info['g395h']['f290lp']['scalerad'] = None

        self.Info['g395h']['opaque'] = {}
        self.Info['g395h']['opaque']['nslices'] = 30
        self.Info['g395h']['opaque']['wscale'] = 0.000665
        self.Info['g395h']['opaque']['ascale'] = 0.1
        self.Info['g395h']['opaque']['bscale'] = 0.1

# ******************************************************************
# Functions

    def SetMultiChannelTable(self, wave, sroi, wroi, power, softrad):
        self.multich_wavelength.append(wave)
        self.multich_sroi.append(sroi)
        self.multich_wroi.append(wroi)
        self.multich_power.append(power)
        self.multich_softrad.append(softrad)
        self.multich_scalerad.append(None)

    def SetMultiChannelEMSMTable(self, wave, sroi, wroi, scalerad):
        self.multich_wavelength.append(wave)
        self.multich_sroi.append(sroi)
        self.multich_wroi.append(wroi)
        self.multich_scalerad.append(scalerad)
        self.multich_power.append(None)
        self.multich_softrad.append(None)

    def SetMultiChannelDrizTable(self, wave):
        self.multich_wavelength.append(wave)
        self.multich_sroi.append(None)
        self.multich_wroi.append(None)
        self.multich_scalerad.append(None)
        self.multich_power.append(None)
        self.multich_softrad.append(None)

    def SetPrismTable(self, wave, sroi, wroi, power, softrad):
        self.prism_wavelength.append(wave)
        self.prism_sroi.append(sroi)
        self.prism_wroi.append(wroi)
        self.prism_power.append(power)
        self.prism_softrad.append(softrad)
        self.prism_scalerad.append(None)

    def SetMedTable(self, wave, sroi, wroi, power, softrad):
        self.med_wavelength.append(wave)
        self.med_sroi.append(sroi)
        self.med_wroi.append(wroi)
        self.med_power.append(power)
        self.med_softrad.append(softrad)
        self.med_scalerad.append(None)

    def SetHighTable(self, wave, sroi, wroi, power, softrad):
        self.high_wavelength.append(wave)
        self.high_sroi.append(sroi)
        self.high_wroi.append(wroi)
        self.high_power.append(power)
        self.high_softrad.append(softrad)
        self.high_scalerad.append(None)

    def SetPrismEMSMTable(self, wave, sroi, wroi, scalerad):
        self.prism_wavelength.append(wave)
        self.prism_sroi.append(sroi)
        self.prism_wroi.append(wroi)
        self.prism_scalerad.append(scalerad)
        self.prism_power.append(None)
        self.prism_softrad.append(None)

    def SetMedEMSMTable(self, wave, sroi, wroi, scalerad):
        self.med_wavelength.append(wave)
        self.med_sroi.append(sroi)
        self.med_wroi.append(wroi)
        self.med_scalerad.append(scalerad)
        self.med_power.append(None)
        self.med_softrad.append(None)

    def SetHighEMSMTable(self, wave, sroi, wroi, scalerad):
        self.high_wavelength.append(wave)
        self.high_sroi.append(sroi)
        self.high_wroi.append(wroi)
        self.high_scalerad.append(scalerad)
        self.high_softrad.append(None)
        self.high_power.append(None)

    def SetXSliceLimits(self, x1, x2, parameter1):
        self.Info[parameter1]['xstart'] = x1
        self.Info[parameter1]['xend'] = x2

    def SetMSM(self, parameter1, parameter2, sroi, wroi, power, softrad):
        self.Info[parameter1][parameter2]['sroi'] = sroi
        self.Info[parameter1][parameter2]['wroi'] = wroi
        self.Info[parameter1][parameter2]['msm_power'] = power
        self.Info[parameter1][parameter2]['softrad'] = softrad
        self.Info[parameter1][parameter2]['scalerad'] = None

    def SetEMSM(self, parameter1, parameter2, sroi, wroi, scalerad):
        self.Info[parameter1][parameter2]['sroi'] = sroi
        self.Info[parameter1][parameter2]['wroi'] = wroi
        self.Info[parameter1][parameter2]['msm_power'] = None
        self.Info[parameter1][parameter2]['softrad'] = None
        self.Info[parameter1][parameter2]['scalerad'] = scalerad

    def SetSpatialSize(self, value, parameter1, parameter2=None):
        if parameter2 is None:
            self.Info[parameter1]['ascale'] = value
            self.Info[parameter1]['bscale'] = value
        else:
            self.Info[parameter1][parameter2]['ascale'] = value
            self.Info[parameter1][parameter2]['bscale'] = value

    def SetSpectralStep(self, value, parameter1, parameter2):
        self.Info[parameter1][parameter2]['wscale'] = value

    def SetWaveMin(self, value, parameter1, parameter2):
        self.Info[parameter1][parameter2]['wavemin'] = value

    def SetWaveMax(self, value, parameter1, parameter2):
        self.Info[parameter1][parameter2]['wavemax'] = value

    def SetSpatialROI(self, value, parameter1, parameter2):
        self.Info[parameter1][parameter2]['sroi'] = value

    def SetWaveROI(self, value, parameter1, parameter2):
        self.Info[parameter1][parameter2]['wroi'] = value

    def SetMSMPower(self, value, parameter1, parameter2):
        self.Info[parameter1][parameter2]['msm_power'] = value

    def SetSoftRad(self, value, parameter1, parameter2):
        self.Info[parameter1][parameter2]['softrad'] = value

    def SetScaleRad(self, value, parameter1, parameter2):
        self.Info[parameter1][parameter2]['scalerad'] = value

    def GetWaveRoi(self, parameter1, parameter2):
        roiw = self.Info[parameter1][parameter2]['wroi']
        return roiw

    def GetSpatialRoi(self, parameter1, parameter2):
        rois = self.Info[parameter1][parameter2]['sroi']
        return rois

    def GetWaveMin(self, parameter1, parameter2):
        wavemin = self.Info[parameter1][parameter2]['wavemin']
        return wavemin

    def GetWaveMax(self, parameter1, parameter2):
        wavemax = self.Info[parameter1][parameter2]['wavemax']
        return wavemax

    def GetMSMPower(self, parameter1, parameter2):
        weight_power = self.Info[parameter1][parameter2]['msm_power']
        return weight_power

    def GetSoftRad(self, parameter1, parameter2):
        softrad = self.Info[parameter1][parameter2]['softrad']
        return softrad

    def GetScaleRad(self, parameter1, parameter2):
        scalerad = self.Info[parameter1][parameter2]['scalerad']
        return scalerad

    def GetScale(self, parameter1, parameter2):
        scale = (self.Info[parameter1][parameter2]['ascale'],
                 self.Info[parameter1][parameter2]['bscale'],
                 self.Info[parameter1][parameter2]['wscale'])
        return scale

    def Get_multichannel_table(self, weighting):
        table = (self.multich_wavelength,
                 self.multich_sroi,
                 self.multich_wroi,
                 self.multich_power,
                 self.multich_softrad,
                 self.multich_scalerad)
        return table

    def Get_prism_table(self):
        table = (self.prism_wavelength,
                 self.prism_sroi,
                 self.prism_wroi,
                 self.prism_power,
                 self.prism_softrad,
                 self.prism_scalerad)
        return table

    def Get_med_table(self):
        table = (self.med_wavelength,
                 self.med_sroi,
                 self.med_wroi,
                 self.med_power,
                 self.med_softrad,
                 self.med_scalerad)
        return table

    def Get_high_table(self):
        table = (self.high_wavelength,
                 self.high_sroi,
                 self.high_wroi,
                 self.high_power,
                 self.high_softrad,
                 self.high_scalerad)
        return table

    def GetMIRISliceEndPts(self, parameter1):
        slice_xstart = self.Info[parameter1]['xstart']
        slice_xend = self.Info[parameter1]['xend']
        return slice_xstart, slice_xend

    def GetStartSlice(self, parameter1):
        sliceno = self.Info[parameter1]['start_slice']
        return sliceno

    def GetEndSlice(self, parameter1):
        sliceno = self.Info[parameter1]['end_slice']
        return sliceno

    def GetNSlice(self, parameter1):
        numslice = self.Info[parameter1]['nslices']
        return numslice
