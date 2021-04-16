#! /usr/bin/env python
#
# mc_3d.py - 2D Monte Carlo for adding cosmic rays to JWST data cube.
#
#
# History: 09/22/09 - first version
#          01/12/10 - Massimo R. pointed out that the flux in the SP spectrum I had been using includes
#                     flares, resulting in a too high SP flux. To correct this, I've adjusted the fluxes
#                     to be 4.2/cm2/sec (GCR) and 0.6/cm2/sec (SP) to sum to his suggested total value of
#                     4.8, although this is (almost surely) not the correct (i.e., non-flare) SP spectrum.
#          04/06/10 - changing to have arbitrary impact directions and energies for incident cosmic rays,
#                     which are sampled from Massimo Robberto's CR library
#          05/19/10 - fixed bug which appears for non-square input arrays
#          08/17/10 - fixed minor bugs; added optional command line parameter 'solar_env' for solar environment type
#                     (SUNMIN, SUNMAX, FLARES); changed mandatory command line parameter 'cr_lib_dir' to be directory
#                     of installed fits files for CR templates; flux value now calculated from solar_env
#          10/22/10 - force output file of signal plus cosmic rays to be float64, irrespective of input file
#
#
# Required Input:
#  'data_cube_file' : data cube file created by create_cube
#  'cr_lib_dir': directory of installed fits files for CR templates to sample from
#  'solar_env': solar environment - allowed values are 'SUNMIN', 'SUNMAX', FLARES'; default = 'SUNMIN'
#
# Optional Input:
#  'verb' : level of verbosity for print statements; verb=2 is extremely verbose, verb=1 is moderate,
#           verb=0 is minimal; default = 0
#  'm_flux' : multiplicative factor for flux; value given is multiplied by the expected fluxes for cosmic
#             rays ; default = 1.0
#  'old_cr_file': file of previously generated GCR/SPs to add to the input data cube; default = None (which
#                 causes new cosmic rays to be generated)
#
# Arrays output to files:
#  'new_cube.fits': input data cube with GCRs and SPs added
#  'cr_and_sp.fits': cube of cosmic rays added
#
# linux command line usage examples:
#   input data cube "acube.fits", verbosity level = 2, flux = 10* 'expected flux' for SOLARMIN, sample from CR templates
#       in directory "cr_data/", generate new CRs:
#   ./mc_3d.py acube.fits "cr_data/" SUNMIN 10. 2
#
#   same but use previously generated cosmic rays in old_cr_and_sp.fits (flux multiplier is specified but will be ignored):
#   ./mc_3d.py acube.fits "cr_data/" SUNMIN 10. 2  "old_cr_and_sp.fits"
#
# pyraf usage examples:
#   input data cube "acube.fits", default verbosity level = 1, default m_flux = 1.*'expected flux', default solar_env SOLARMIN,
#       sample from CR templates in directory "cr_data/", generate new CRs:
#   aa= mc_3d.cr_add_3d("acube.fits", "cr_data/")
# --> status = aa.add_crs_all_slices()
#
#   input data cube "acube.fits", verbosity level = 2, solar_env = SUNMAX, sample from CR templates in directory "cr_data/",
#        use previously generated cosmic rays in old_cr_and_sp.fits (flux multiplier is specified but will be ignored):
# --> aa= mc_3d.cr_add_3d("acube.fits", "cr_data/","SUNMAX", 1., 2,"old_cr_and_sp.fits")
# --> status = aa.add_crs_all_slices()
#

import math
import random
import sys
import time

from astropy.io import fits
import numpy as np

ERROR_RETURN = 2

MAX_CR = 5  # max number of GCRs allowed per pixel in a given integration

E_BAND = 0.9 # eV per ionization: approximate value, get better one later
EXP_TIME = 10.6  # exposure time in seconds

class cr_add_3d:
    """ Add arbitrarily angled incident cosmic rays to input data cube
    """

    def __init__(self, data_cube_file, cr_lib_dir, solar_env='SUNMIN', m_flux=1.0, verb=0, old_cr_file=None):
        """ Constructor
        @param data_cube_file: input data cube to add GCR/SPs to
        @type data_cube_file: string
        @param cr_lib_dir: type of CR library CR fluxes
        @type cr_lib_dir: string
        @param solar_env: type of solar environment
        @type solar_env: string
        @param m_flux: multiplicative factor for flux
        @type m_flux: float
        @param old_cr_file: file of previously generatedd GCR/SPs to add to the input data cube
        @type old_cr_file: string
        @param verb: verbosity level
        @type verb: integer
        """

        if (solar_env == 'SUNMIN'):
            cr_flux = 4.8983E-8 # expected flux in number/(micron^2)/sec for 100 mil Al shield, solar minimum
        elif (solar_env == 'SUNMAX'):
            cr_flux = 1.7783E-8 # expected flux in number/(micron^2)/sec for 100 mil Al shield, solar maximum
        elif (solar_env == 'FLARES'):
            cr_flux = 3046.83E-8 # expected flux in number/(micron^2)/sec for 100 mil Al shield, solar flare
        else:
            print(' Invalid solar environment specified, so will use solar min ')
            solar_env = 'SUNMIN'
            cr_flux = 4.8983E-8 # expected flux in number/(micron^2)/sec for 100 mil Al shield, solar minimum

        self._data_cube_file = data_cube_file
        self._cr_lib_dir = cr_lib_dir + 'CRs_MCD5.5_' + solar_env
        self._verb = int(verb)
        self._old_cr_file = old_cr_file
        self._solar_env = solar_env
        self._cr_flux = float(m_flux) * float(cr_flux)

        if (self._verb > 0):
            print('Input parameters:')
            print('   input data cube file = ', data_cube_file)
            print('   m_flux = ', m_flux)
            print('   verb = ', verb)
            print('   cr_flux = ', self._cr_flux, '/micron^2/sec, (recommended',
                'value for this solar activity level with 100 mil Al shielding) ')
            print('   cr_lib_dir = ', cr_lib_dir)
            print('   solar_env = ', solar_env)
            print('   old CR file = ', old_cr_file)
            print('  ')

    def add_crs_all_slices(self):
        """ Generate and add GCR/SPs for all slices in cube
        """

        verb = self._verb
        self.fh_data_cube = fits.open(self._data_cube_file, mode='update')
        prihdr = self.fh_data_cube[0].header

        total_slices = self.fh_data_cube[0].data.shape[0]
        asize_1 = self.fh_data_cube[0].data.shape[1]
        asize_2 = self.fh_data_cube[0].data.shape[2]

        self._asize_1 = asize_1; self._asize_2 = asize_2

        if (verb > 0):
            print('   ')
            print(' Data cube info: ', self.fh_data_cube.info())
            print(' Array dimensions :', asize_1, asize_2)
            print(' Number of slices :', total_slices)
            try:
                print(' READMODE specified in the input header : ', prihdr['READMODE'])
            except Exception:
                print('WARNING: unable to determine readmode from header ')

        # read 'NAVGIMAG'] from input header to later output to output header
        try:
            self._navg = self.fh_data_cube[0].header['NAVGIMAG']
            if (verb > 0):
                print(' From the input header keyword NAVGIMAG, navg = ', self._navg)
                print('  ')
        except Exception:
            print('WARNING: unable to access NAVGIMAG from primary header of',
                'data cube, so assuming that the number of samples averaged',
                'over is 1 - which may not be what you want ! ')
            self._navg = 1

        self.fh_data_cube[0].data = (self.fh_data_cube[0].data).astype(np.float64)

        # generate new GCRs and SPs
        if self._old_cr_file is None:
            if (verb > 0): print(' Generating GCRs and SPs...')

            try:
                self._instrume = prihdr['INSTRUME']
                print(' This data cube is for instrument : ', self._instrume)
            except Exception:
                print('WARNING: unable to determine the instrument from header, will default to NIRCAM ')
                self._instrume = 'NIRCAM'

            if (self._instrume == 'NIRCAM'):
                self._electron_per_adu = 1.3
                self._pix_length = 18.
            elif (self._instrume == 'MIRI'):
                self._electron_per_adu = 2.2
                self._pix_length = 25.
            else:
                print('Fatal ERROR - unsupported instrument')

            self._half_pix_length = self._pix_length / 2.
            self._pix_area = self._pix_length**2.  # area in micron2

            self._grand_total_gcr = 0

            self._tot_gcr_and_sp = np.zeros((asize_2, asize_1), dtype=np.float64)

            # data cube of created crs and sps :
            self.cr_and_sp_only_cube = np.zeros((total_slices, asize_2, asize_1), dtype=np.float64)

            # generate number of cosmic rays
            exp_gcr = self._cr_flux * self._pix_area * EXP_TIME  # expected number of cosmic rays
            self._cumul_pois_g = calc_poisson(exp_gcr) # generate distribution of number of CR

            if (verb > 0):
                print('Instrument: ', self._instrume)
                print(' cosmic ray flux = ', self._cr_flux, '/micron^2/sec')
                print('Expected number of  cosmic rays per pixel = ', exp_gcr)
                print('Pixel length = ', self._pix_length, ' microns ')
                print('Pixel area = ', self._pix_area, 'microns^2')
                print('Exposure time = ', EXP_TIME, 'seconds')
                print('  ')

            for ii_slice in range(total_slices):
                if (verb > 0):
                    print('   ')
                    print('  ... beginning this slice:  ', ii_slice)
                self._new_gcr_and_sp = np.zeros((asize_2, asize_1), dtype=np.float64)
                self._slice = ii_slice
                self.add_crs_this_slice()  # add all GCRs for this slice

            self.write_file(self.fh_data_cube[0].data[:, :, :], 'new_cube.fits')
            self.write_file(self.cr_and_sp_only_cube[:, :, :], 'cr_and_sp.fits')

            if (verb > 0):
                print(' Wrote data cube of sources plus CRs to new_cube.fits')
                print(' Wrote cube of created CRs to cr_and_sp.fits')
                print(' Total number of GCRs summed over all slices = ', self._grand_total_gcr)

        else:  #  self._old_cr_file <> None  so will use previously generated CR/SP file
            self.fh_old_cr = fits.open(self._old_cr_file)

            if (verb > 0):
                print(' Will be using previously generated GCRs and SPs from the specified file: ')
                print(' ... with info: ', self.fh_old_cr.info())

            for ii_slice in range(1, total_slices):     # add cr slice 0 to data slice 1, etc
                self._slice = ii_slice

                slice_old_cr = self.fh_old_cr[0].data[self._slice - 1, :, :]  # note -1
                if (verb > 1):
                    print(' processing slice: ', slice_old_cr)

            #  need to add this old cr slice to all later data slices
                for jj_slice in range(ii_slice, total_slices):
                    self.fh_data_cube[0].data[jj_slice, :, :] += slice_old_cr

            self.write_file(self.fh_data_cube[0].data[:, :, :], 'new_cube.fits')
            print(' Wrote data cube of sources plus GCRs and SPs to new_cube.fits')



    def add_crs_this_slice(self):
        """ Add all GCR/SPs to the current slice
        """

        verb = self._verb

        if (verb > 1): print(' Start of add_crs_this_slice() for slice number: ', self._slice)
        try:
            pix_val = self.fh_data_cube[0].data[self._slice, :, :]  # do all pixels in this slice
            input_data = pix_val.copy()

            if (self._slice > 0): # add cr/sp for all but 0th slice
                tot_nz_gcr_and_sp = np.where(self._tot_gcr_and_sp > 0.0)
                pix_val += self._tot_gcr_and_sp

                if (verb > 0):
                    print(' the number of pixels affected by gcr and sp for slice', self._slice, ' = ', len(tot_nz_gcr_and_sp[0]))
                if (verb > 1):
                    print(' the amplitudes of CR/SP added for this slice: ', self._tot_gcr_and_sp[tot_nz_gcr_and_sp])


        except Exception as errmess:
            print('Fatal ERROR in add_crs_this_slice() : ', errmess)
            sys.exit(ERROR_RETURN)

        try:
            self.fh_data_cube[0].header
        except Exception:
            print('Fatal ERROR : unable to open primary header of data cube')
            sys.exit(ERROR_RETURN)

        #  generate number of  crs in this pixel
        total_gcr = 0  # total number of cosmic rays hitting array
        tot_pix_done = 0  # number of pixels looped over

        asize_1 = self._asize_1; asize_2 = self._asize_2

        for ii_pix in range(0, self._asize_1):
            for jj_pix in range(0, self._asize_2):

                tot_pix_done += 1

                if (verb > 1):
                    print('Pixel [', ii_pix, ',', jj_pix, '] initially has value:', pix_val[ii_pix, jj_pix])

                # generate CRs
                got_num = 0 # 0 for not yet generated number of  crs, reset to 1 when generated
                prob = random.uniform(0.0, 1.0)
                ii_cr = 0

                while (ii_cr < MAX_CR) and (got_num == 0):
                    if (prob < self._cumul_pois_g[ii_cr]):
                        num_gcr = ii_cr
                        got_num = 1
                    ii_cr += 1

                if (verb > 1 and num_gcr > 0): print(' The number of cosmic rays incident on this pixel will be ', num_gcr)

                # for each cr for this pixel, determine energy deposited in primary and secondary affected pixels
                for each_cr in range(num_gcr):
                    total_gcr += 1
                    self._grand_total_gcr += 1

                    # Get energy and template from library for CR affecting this primary pixel (ii_pix, jj_pix), and
                    #    through interpixel capacitance coupling, neighboring pixels. First generate random file.
                    file_num = int(10 * random.uniform(0.0, 1.0))

                    # use 'raw' (not interpixel coupling) version of files
                    fname = self._cr_lib_dir + "_0" + str(file_num) + ".fits"

                    if (verb > 1):
                        print('  ... a cosmic ray will be grabbed from file : ', fname)

                    fh = fits.open(fname)
                    cr_cube_data = fh[1].data
                    cr_cube_energy = fh[3].data

                    slice_num = int(1000 * random.uniform(0.0, 1.0))  # 1000 slices in library file

                    if (verb > 1):
                        print(' ... the energy deposited by this cosmic ray :', cr_cube_energy[slice_num], ' MeV')

                    cr_slice_data = cr_cube_data[slice_num, :, :]

                    # center pixel of CR template
                    x_ctr_pix = int(cr_slice_data.shape[0] / 2)
                    y_ctr_pix = int(cr_slice_data.shape[0] / 2)

                    if (verb > 1):
                        print('  The CR template at its center is cr_slice_data'
                            '[y_ctr_pix, x_ctr_pix] =', cr_slice_data[y_ctr_pix, x_ctr_pix])
                        print('  The CR max deposited energy  =', cr_slice_data.max())

                    wh_nz = np.where(cr_slice_data > 0.0)
                    num_pix_affected = len(cr_slice_data[wh_nz])

                    if (verb > 1):
                        print(' The number of secondary pixels potentially affected by this cosmic ray', num_pix_affected)
                        print(' ... and the energies deposited in these neighboring pixels : ', cr_slice_data[wh_nz])

                    num_2nd_aff = 0 # number of 2ndary pixels affected (including edge effects)
                    for qq in range(0, num_pix_affected - 0): # loop over affected secondary pixels
                        if (verb > 1):
                            print('   ')
                            print('  For secondary potentially affected pixel: ', qq)

                        x_tpl = wh_nz[0][qq]  # template-centered coordinates
                        y_tpl = wh_nz[1][qq]  # template-centered coordinates

                        if (verb > 1):
                            print(' The template-centered coordinates for this',
                                'secondary affected pixel are (', x_tpl, ',', y_tpl, ')')

                    # the secondary pixel coordinates to which energy is added are xval, yval, which are :
                        xval = ii_pix - x_ctr_pix + x_tpl
                        yval = jj_pix - y_ctr_pix + y_tpl

                        # start check that xval, yval are within array
                        beyond_edge = False  # not beyond edge
                        if ((xval >= asize_1) or (yval >= asize_2) or (xval < 0) or (yval < 0)):
                            if (verb > 1):
                                print('  this 2ndary pixel is beyond the array edge so this 2ndary pixel does not exist')
                            beyond_edge = True

                        if (beyond_edge == False):  # include this 2ndary pixel
                            num_2nd_aff += 1 # number of 2ndary pixels affected (including edge effects)
                            if (verb > 1):
                                print(' the pixel coordinates to which energy is added are xval, yval, which are :', xval, yval)
                                print('    will include this secondary pixel, just incremented num_2nd_aff to be', num_2nd_aff)
                            delta_e_mev = cr_slice_data[wh_nz][qq]  # energy deposited in this pixel in MeV

                            # convert this delta_e_mev from MeV to ADU before adding it to pix_val which is in ADU
                            delta_e_adu = (delta_e_mev * 1E6 / E_BAND) / (self._electron_per_adu) # ... in ADU

                            if (verb > 1):
                                print('  the energy deposited by this GCR into',
                                    'this pixels is', delta_e_mev, ' MeV, which is',
                                    delta_e_adu, ' ADU')
                                print('  ')

                            # add energy deposit to arrays for this pixel
                            self._new_gcr_and_sp[yval, xval] += delta_e_adu
                            # to compare to detections in a later program
                            self.cr_and_sp_only_cube[self._slice, yval, xval] += delta_e_adu
                        else:
                            if (verb > 1): print(' will *NOT* include this 2ndary pixel - beyond edge')

                    if (verb > 1): print(' for this CR each_cr = ', each_cr,
                        ' the number of secondary  pixels affected (including',
                            'edge effects) = num_2nd_aff =', num_2nd_aff)


        if (verb > 1):
            print(' just before adding the latest GCR and SP, the min, mean,',
                'and max (over the slice) of the cumulative total of the',
                'amplitudes of GCR and SP for the current slice: ',
                self._tot_gcr_and_sp.min(), self._tot_gcr_and_sp.mean(),
                self._tot_gcr_and_sp.max())

        self._tot_gcr_and_sp += self._new_gcr_and_sp  #  add all pixels of this slice to total

        if (verb > 1):
            print(' just after adding the latest GCR and SP, the min, mean,',
                'and max (over the slice) of the cumulative total of the',
                'amplitudes of GCR and SP for the current slice: ',
                self._tot_gcr_and_sp.min(), self._tot_gcr_and_sp.mean(),
                self._tot_gcr_and_sp.max())
            print('         ')
            print(' original pix_val in ADU = ', input_data)
            print(' final pix_val in ADU : ', pix_val)

        if (verb > 1):
            final_minus_original_data = pix_val - input_data
            print(' final minus original in ADU for this slice (slice number',
                self._slice, '): ', final_minus_original_data)
            print(' ... final minus original has min, mean, max : ',
                final_minus_original_data.min(), final_minus_original_data.mean(),
                final_minus_original_data.max())

        if (verb > 1):
            print(' the total number of gcr hitting this slice (slice number', self._slice, '): ', total_gcr)
            print(' the final array values for this slice : ', self.fh_data_cube[0].data[self._slice, :, :])

        if (verb > 1): print(' end of add crs for this slice')


    def write_file(self, data, output_fname):
        """ Write data cube to file, with relevant keywords
        @param data: data cube created
        @type data: floats
        @param output_fname: name of output file
        @type output_fname: string
        """
        fitsobj = fits.HDUList()
        hdu = fits.PrimaryHDU()

        prihdr = hdu.header
        prihdr['cr_flux'] = (self._cr_flux, 'cr_flux, /(micron^2*sec)')
        prihdr['NAVGIMAG'] = (self._navg, 'images averaged from create_cube')

        hdu.data = data
        fitsobj.append(hdu)
        fitsobj.writeto(output_fname, overwrite=True)
        fitsobj.close()


def calc_poisson(exp_num):
    """  Calculate cumulative poisson distribution for number of GCR/SPs given expected number
    @param exp_num: expected number of incident particles (GCR or SP) for this pixel
    @type exp_num: integer
    @return: normalized cumulative poisson distribution
    @rtype: float array
    """

    pois = np.zeros(MAX_CR, dtype=np.float64)  # only allow a max of MAX_CR particles/pixel/integration
    cumul_pois = np.zeros(MAX_CR, dtype=np.float64)

    for ii in range(0, MAX_CR):
        if ii > 0:
            f1 = lookup_factorial(ii)
        else:
            f1 = 1.
        pois[ii] = (exp_num**ii) * math.e**(-exp_num) / f1
        if ii > 0:
            cumul_pois[ii] += cumul_pois[ii - 1] + pois[ii]
        else:
            cumul_pois[ii] += pois[ii]

    return cumul_pois


def lookup_factorial(kk):

    """ lookup factorial table for small arguments
        @param kk: input value
        @type kk: integer
        @return: lookup factorial of input
        @rtype: integer
        """
    if kk == 0: return 1
    elif kk == 1: return 1
    elif kk == 2: return 2
    elif kk == 3: return 6
    elif kk == 4: return 24
    elif kk == 5: return 120
    elif kk == 6: return 720
    elif kk == 7: return 5040
    elif kk == 8: return 40320
    elif kk == 9: return 362880
    elif kk == 10: return 3628800
    else:
        print(' argument unexpectedly large in factorial() ')
        sys.exit(ERROR_RETURN)

if __name__ == "__main__":
    """ Get datacube, optional parameters, and add GCR/SP to datacube
    """
    usage = "USAGE: mc_3d.py datacube cr_lib_dir [solar_env [m_flux] [old_cr_file] [[[verb]]]] "

    if len(sys.argv) < 2:
        print("syntax: mc_3d.py datacube cr_lib_dir [solar_env [m_flux] [verb] [[[old_cr_file]]]] ")
        sys.exit(ERROR_RETURN)

    data_cube_file = sys.argv[1]
    cr_lib_dir = sys.argv[2]

    if (len(sys.argv) > 3):
        solar_env = sys.argv[3]
    else:
        solar_env = "SUNMIN"

    if (len(sys.argv) > 4):
        m_flux = sys.argv[4]
    else:
        m_flux = 1.0

    if (len(sys.argv) > 5):
        verb = int(sys.argv[5])
    else:
        verb = int(1)

    if (len(sys.argv) > 6):
        old_cr_file = sys.argv[6]
    else:
        old_cr_file = None

    try:
        tstart0 = time.time()

        if (verb > 0):
            print('Starting the program mc_3d: ')
            print('The current time (start) = ', time.asctime())
        aa = cr_add_3d(data_cube_file, cr_lib_dir, solar_env=solar_env, m_flux=m_flux, verb=verb, old_cr_file=old_cr_file)
        status = aa.add_crs_all_slices()

        tstop = time.time()

        if (verb > 0):
            print('.... completing the program mc_3d ')
            print('The current time (stop) = ', time.asctime())

    except Exception as errmess:
        print(' Fatal ERROR in main(): ', errmess)
