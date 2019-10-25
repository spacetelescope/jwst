#! /usr/bin/env python
#
# mc_2d.py - simplified 2D Monte Carlo for adding cosmic rays and solar particles to JWST data cube.
#                   In this version GCRs and SPs are normally incident to the centers of pixels, so each
#                   GCR/SP only affects a single pixel
#
# History: 09/22/09 - first version
#          01/12/10 - Massimo R. pointed out that the flux in the SP spectrum I had been using includes
#                     flares, resulting in a too high SP flux. To correct this, I've adjusted the fluxes
#                     to be 4.2/cm2/sec (GCR) and 0.6/cm2/sec (SP) to sum to his suggested total value of
#                     4.8, although this is (almost surely) not the correct (i.e., non-flare) SP spectrum.
#
#  Energy spectum for SP is from "The Radiation Environment for the NGST" (September 2000)
#  Energy spectum for GCR is from figure 4-2 of "JWST Mission Operations Document" (JWST-OPS-002018, Revision A)
#  Energy deposited (stopping power) based on figure 2-7 (p32) of Segre's "Nuclei and Particles", 2nd edition
#
# Required Input:
#  'data_cube_file' : data cube file created by create_cube
#
# Optional Input:
#  'verb' : level of verbosity for print statements; verb=2 is extremely verbose, verb=1 is moderate,
#           verb=0 is minimal; default = 0
#  'm_flux' : multiplicative factor for flux; value given is multiplied by the expected fluxes for galactic cosmic
#             rays and solar particles; default = 1.0
#  'old_cr_file': file of previously generated GCR/SPs to add to the input data cube; default = None (which
#                 causes new GCR/SPs to be generated)
#
# Arrays output to files:
#  'new_cube.fits': input data cube with GCRs and SPs added
#  'cr_and_sp.fits': cube of GCRs and SPs added
#
# linux command line usage for a verbosity level = 2 and flux = 10* 'expected flux', generate new CRs:
# ./mc_2d.py output/021710a_regtest_30_noisy_10times_CR/acube.fits 2 10.
#
# linux command line usage for a verbosity level = 2 and flux = 10* 'expected flux', use file of previously generated CRs:
# ./mc_2d.py output/021710a_regtest_30_noisy_10times_CR/acube.fits 2 10. given_cr_and_sp.fits
#
# pyraf usage
# --> aa= mc_2d.cr_add_2d("output/021710a_regtest_30_noisy_10times_CR/acube.fits",2,10)
# --> status = aa.add_crs_all_slices()
#

import energy_dists
import math
import random
import sys
import time

from astropy.io import fits
import numpy as np

global GCR_FLUX
global SP_FLUX

ERROR_RETURN = 2

HALF_PI = math.pi / 2.0
PI = math.pi

MAX_CR = 5  # max number of GCRs allowed per pixel in a given integration
MAX_SP = 5  # max number of SPs allowed per pixel in a given integration

PIX_HEIGHT = 8. # height of pixel in microns

E_BAND = 0.9 # eV per ionization: approximate value, get better one later ;
EXP_TIME = 10.6  # exposure time in seconds; change to get from header ?

GCR_FLUX = 4.2E-8 # expected galactic cosmic ray flux in number/(micron^2)/sec
SP_FLUX = 0.6E-8 # expected solar particle flux in number/(micron^2)/sec

class cr_add_2d:
    """ Add normally incident cosmic rays and solar particles to input data cube
    """

    def __init__(self, data_cube_file, verb=0, m_flux=1.0, old_cr_file=None):
        """ Constructor
        @param data_cube_file: input data cube to add GCR/SPs to
        @type data_cube_file: string
        @param verb: verbosity level
        @type verb: integer
        @param m_flux: multiplicative factor for flux
        @type m_flux: float
        @param old_cr_file: file of previously generatedd GCR/SPs to add to the input data cube
        @type old_cr_file: string
        """

        self._data_cube_file = data_cube_file
        self._verb = int(verb)
        self._gcr_flux = float(m_flux) * float(GCR_FLUX)
        self._sp_flux = float(m_flux) * float(SP_FLUX)
        self._old_cr_file = old_cr_file

        if (self._verb > 0):
            print('Input parameters:')
            print('   input data cube file = ', data_cube_file)
            print('   m_flux = ', m_flux)
            print('   verb = ', verb)
            print('   gcr_flux = ', self._gcr_flux, '/micron^2/sec')
            print('   sp_flux = ', self._sp_flux, '/micron^2/sec')
            print('   old CR file = ', old_cr_file)
            print('  ')
            print('The energy distribution (energy in MeV, cumulative prob) for the galactic cosmic rays is: ')
            print(energy_dists.gcr_list)
            print('The energy distribution (energy in MeV, cumulative prob) for the solar particles is: ')
            print(energy_dists.sp_list)
            print('  ')

    def add_crs_all_slices(self):
        """ Generate and add GCR/SPs for all slices in cube
        """

        verb = self._verb
        self.fh_data_cube = fits.open(self._data_cube_file, mode='update')
        prihdr = self.fh_data_cube[0].header

        total_slices = self.fh_data_cube[0].data.shape[0]

        asize_2 = self.fh_data_cube[0].data.shape[1]
        asize_1 = self.fh_data_cube[0].data.shape[2]
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
            print('Fatal ERROR : unable to access NAVGIMAG fromx primary header of data cube')
            sys.exit(ERROR_RETURN)

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
            self._grand_total_sp = 0

            self._tot_gcr_and_sp = np.zeros((asize_1, asize_2), dtype=np.float32) # cumulative total for all slices

            # data cube of created crs and sps :
            self.cr_and_sp_only_cube = np.zeros((total_slices, asize_1, asize_2), dtype=np.float32)

            # generate number of galactic cosmic rays
            exp_gcr = self._gcr_flux * self._pix_area * EXP_TIME  # expected number of galactic cosmic rays
            self._cumul_pois_g = calc_poisson(exp_gcr) # generate distribution of number of GCR

            # generate number of solar particles
            exp_sp = self._sp_flux * self._pix_area * EXP_TIME  # expected number of solar particles
            self._cumul_pois_s = calc_poisson(exp_sp) # generate distribution of number of SP

            if (verb > 0):
                print('Instrument: ', self._instrume)
                print('Solar particle flux = ', self._sp_flux, '/micron^2/sec')
                print('Galactic cosmic ray flux = ', self._gcr_flux, '/micron^2/sec')
                print('Expected number of solar particles per pixel = ', exp_sp)
                print('Expected number of galactic cosmic rays per pixel = ', exp_gcr)
                print('Pixel length = ', self._pix_length, ' microns ; pixel height = ', PIX_HEIGHT, ' microns')
                print('Pixel area = ', self._pix_area, 'microns^2')
                print('Exposure time = ', EXP_TIME, 'seconds')
                print('  ')

            if (verb > 1):
                print('  The normalized cumulative probability distribution for solar particles : ', self._cumul_pois_s)
                print('  The normalized cumulative probability distribution for galactic cosmic rays : ', self._cumul_pois_g)


            for ii_slice in range(total_slices):
                if (verb > 1):
                    print('   ')
                    print('  ... beginning this slice:  ', ii_slice)
                self._new_gcr_and_sp = np.zeros((asize_1, asize_2), dtype=np.float32)  # .... for the current slice
                self._slice = ii_slice
                self.add_crs_this_slice()  # add all GCR/SPs for this slice

            self.write_file(self.fh_data_cube[0].data[:, :, :], 'new_cube.fits')
            self.write_file(self.cr_and_sp_only_cube[:, :, :], 'cr_and_sp.fits')

            if (verb > 0):
                print(' Wrote data cube of sources plus CRs to new_cube.fits')
                print(' Wrote cube of created CRs to cr_and_sp.fits')
                print(' Total number of GCRs summed over all slices = ', self._grand_total_gcr)
                print(' Total number of SPs summed over all slices = ', self._grand_total_sp)

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
                pix_val += self._tot_gcr_and_sp   # add the gcrs and sps from ALL previous slices

                if (verb > 0):
                    print('  number of pixels affected by gcr and sp for slice', self._slice, '  = ', len(tot_nz_gcr_and_sp[0]))
                if (verb > 1):
                    print('  amplitudes of CR/SP added for this slice: ', self._tot_gcr_and_sp[tot_nz_gcr_and_sp])


        except Exception as errmess:
            print('Fatal ERROR in add_crs_this_slice() : ', errmess)
            sys.exit(ERROR_RETURN)

        try:
            self.fh_data_cube[0].header
        except Exception:
            print('Fatal ERROR : unable to open primary header of data cube')
            sys.exit(ERROR_RETURN)

        #  generate number of galactic crs in this pixel
        total_gcr = 0  # total number of cosmic rays hitting array
        total_sp = 0  # total number of solar particles hitting array
        tot_pix_done = 0  # number of pixels looped over

        # force CR to be normally incident at center of pixel
        self._in_x = self._half_pix_length
        self._in_y = self._half_pix_length
        self._in_z = PIX_HEIGHT

        for ii_pix in range(0, self._asize_1):
            for jj_pix in range(0, self._asize_2):
                tot_pix_done += 1

                if (verb > 1):
                    print(' Pixel [', ii_pix, ',', jj_pix, '] initially has value:', pix_val[ii_pix, jj_pix])

                # generate and trace GCRs
                got_num = 0 # 0 for not yet generated number of galactic crs, reset to 1 when generated
                prob = random.uniform(0.0, 1.0)
                ii_cr = 0

                while (ii_cr < MAX_CR) and (got_num == 0):
                    if (prob < self._cumul_pois_g[ii_cr]):
                        num_gcr = ii_cr
                        got_num = 1
                    ii_cr += 1

                for each_cr in range(num_gcr): # for each galactic cr for this pixel, do ray tracing.....
                    total_gcr += 1
                    self._grand_total_gcr += 1

                    # determine which energy bin is to be assigned to ; energy is in MeV
                    energy = self.energy_dist('gcr')

                    if (verb > 1): print(' There are ', num_gcr, ' galactic cosmic rays for this pixel ...')

                    delta_e_mev = self.calculate_energy_dep(ii_pix, jj_pix, energy) # ... is energy deposited in Mev

                    # convert this delta_e_mev from MeV to ADU before adding it to pix_val which is in ADU
                    delta_e_adu = (delta_e_mev * 1E6 / E_BAND) / (self._electron_per_adu) # ... in ADU

                    if (verb > 1):
                        print('  the energy deposited by this GCR into this pixels is', delta_e_adu, ' ADU')
                        print('  ')

                    # add energy deposit to arrays for this pixel
                    self._new_gcr_and_sp[jj_pix, ii_pix] += delta_e_adu
                    # to compare to detections in a later program
                    self.cr_and_sp_only_cube[self._slice, jj_pix, ii_pix] += delta_e_adu

                # generate and trace SPs
                got_num = 0 # 0 for not yet generated number of solar particles, reset to 1 when generated
                prob = random.uniform(0.0, 1.0)
                ii_sp = 0

                while (ii_sp < MAX_SP) and (got_num == 0):
                    if (prob < self._cumul_pois_s[ii_sp]):
                        num_sp = ii_sp
                        got_num = 1
                    ii_sp += 1

                for each_sp in range(num_sp): # for each solar particle for this pixel, do ray tracing.....
                    total_sp += 1
                    self._grand_total_sp += 1

                # determine which energy bin is to be assigned to ; energy is in MeV
                    energy = self.energy_dist('sp')

                    if (verb > 1): print(' There are ', num_sp, ' solar particles for this pixel ...')

                    delta_e_mev = self.calculate_energy_dep(ii_pix, jj_pix, energy) # .. is energy deposited in MeV

                # convert this delta_e from MeV to ADU before adding it to pix_val which is in ADU
                    delta_e_adu = (delta_e_mev * 1E6 / E_BAND) / (self._electron_per_adu) # ... in ADU

                    if (verb > 1):
                        print('  the energy deposited by this SP into this pixels is', delta_e_adu, ' ADU')
                        print('   ')

                    # add energy deposit to arrays for this pixel
                    self._new_gcr_and_sp[jj_pix, ii_pix] += delta_e_adu
                    self.cr_and_sp_only_cube[self._slice, jj_pix, ii_pix] += delta_e_adu

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
            print(' the total number of solar particles hitting this slice: ', total_sp)
            print(' the final array values for this slice : ', self.fh_data_cube[0].data[self._slice, :, :])

        if (verb > 1): print(' end of add crs for this slice')


    def calculate_energy_dep(self, ii_pix, jj_pix, energy):
        """ Calculate the energy deposited in the pixel for this GCR/SP
        @param ii_pix: x-coordinate of pixel within array
        @type ii_pix: integer
        @param jj_pix: y-coordinate of pixel within array
        @type jj_pix: integer
        @param energy: energy of incident particle in MeV
        @type energy: float
        """

        self._de_by_dx = calc_stopping_power(energy)

        # propagate ray thru pixel it is incident upon ....
        len_pix = self.ray_trace(ii_pix, jj_pix, energy)

        delta_e_mev = self._de_by_dx * len_pix   # this pixel's energy increase in MeV

        if (self._verb > 1):
            print(' Ray trace() completed for this CR/SP ')
            print(' output coordinates out_x = ', self._out_x, ' self._out_y = ', self._out_y)
            print(' de_by_dx = ', self._de_by_dx, ' which gives this pixels energy increase in MeV = ', delta_e_mev)
            print(' The number of pixels affected by this cosmic ray = 1 ') # forced

        return delta_e_mev  # energy deposited in this pixel in MeV


    def ray_trace(self, ii_pix, jj_pix, energy):
        """  Ray trace through pixel to calculate path length of incident particle within pixel
        (simplified case since the particles are normally incident at the center of the pixel)
        @param ii_pix: x-coordinate of pixel within array
        @type ii_pix: integer
        @param jj_pix: y-coordinate of pixel within array
        @type jj_pix: integer
        @param energy: energy of incident particle in MeV
        @type energy: float
        @return: path length of incident particle within pixel
        @rtype: float
        """

        verb = self._verb
        self._out_x = self._in_x; self._out_y = self._in_y; self._out_z = 0.0

        len_pix = self._in_z # length of interaction

        if (verb > 1):
            print(' Start of ray_trace for pixel with array coords: x = ', ii_pix, ', y = ', jj_pix)
            print('  - within pixel, impact coords are in_x = ', self._in_x, ', in_y = ', self._in_y, ', in_z = ', self._in_z)
            print('  - CR/SP energy =  ', energy, ' MeV')

        return len_pix


    def energy_dist(self, dist_type):
        """ sample from GCR or SP energy distribution and return particle energy in MeV
        @param dist_type: energy distribution type ('gcr', or 'sp')
        @type dist_type: string
        @return: particle energy in MeV
        @rtype: float
        """
        verb = self._verb

        if (dist_type == 'gcr'):
            part_list = energy_dists.gcr_list  # each list entry contains : index, energy, probability
        else:
            part_list = energy_dists.sp_list  # each list entry contains : index, energy, probability

        len_part_list = len(part_list)
        part_prob = random.uniform(0.0, 1.0)

        for ii_l in range(len_part_list):
            if (part_prob < part_list[ii_l][1]):
                if (ii_l == 0):
                    return part_list[0][0]
                elif (ii_l == len_part_list):
                    return part_list[len_part_list - 1][0]
                else:
                    return part_list[ii_l - 1][0]

        if (verb > 1):
            print(' For a cr/sp of type ', dist_type, ' the generated probability = ', part_prob)

        return part_list[len_part_list - 1][0] # is in MeV



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
        prihdr['GCR_FLUX'] = (self._gcr_flux, 'GCR_FLUX, /(micron^2*sec)')
        prihdr['SP_FLUX'] = (self._sp_flux, 'SP_FLUX, /(micron^2*sec)')
        prihdr['NAVGIMAG'] = (self._navg, 'images averaged from create_cube')

        hdu.data = data
        fitsobj.append(hdu)
        fitsobj.writeto(output_fname, overwrite=True)
        fitsobj.close()

# Calculate stopping power in MeV/micron for this energy in Silicon; removed density dependence
# Derived from figure 2.7 of Segre, rescaled to apply to Silicon. This applies to both NIRCAM and MIRI
def calc_stopping_power(energy):
    """  Calculate stopping power in MeV/micron for this incident energy in Silicon
    @param energy: energy of incidemt particle in MeV
    @type energy: float
    @return: stopping power in MeV/micron
    @rtype: float
    """
    if (energy < 1.):
        de_by_dx = 0.060
    elif (energy > 1000.):
        de_by_dx = 0.060 * (1000. ** (-0.684))
    else:
        de_by_dx = 0.060 * (energy ** (-0.684))

    return de_by_dx


def calc_poisson(exp_num):
    """  Calculate cumulative poisson distribution for number of GCR/SPs given expected number
    @param exp_num: expected number of incident particles (GCR or SP) for this pixel
    @type exp_num: integer
    @return: normalized cumulative poisson distribution
    @rtype: float array
    """

    pois = np.zeros(MAX_CR, dtype=np.float32)  # only allow a max of MAX_CR particles/pixel/integration
    cumul_pois = np.zeros(MAX_CR, dtype=np.float32)

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
    usage = "usage: ./mc_2d datacube [verb [m_flux] [[old_cr_file]]] "

    if len(sys.argv) < 2:
        print("syntax: ./mc_2d datacube [verb [m_flux] [[old_cr_file]]] ")
        sys.exit(ERROR_RETURN)

    data_cube_file = sys.argv[1]

    if (len(sys.argv) > 2):
        verb = int(sys.argv[2])
    else:
        verb = int(1)

    if (len(sys.argv) > 3):
        m_flux = sys.argv[3]
    else:
        m_flux = 1.0

    if (len(sys.argv) > 4):
        old_cr_file = sys.argv[4]
    else:
        old_cr_file = None

    try:
        tstart0 = time.time()

        if (verb > 0):
            print('Starting the program mc_2d: ')
            print('The current time (start) = ', time.asctime())
        aa = cr_add_2d(data_cube_file, verb=verb, m_flux=m_flux, old_cr_file=old_cr_file)
        status = aa.add_crs_all_slices()

        tstop = time.time()

        if (verb > 0):
            print('.... completing the program mc_2d ')
            print('The current time (stop) = ', time.asctime())

    except Exception as errmess:
        print(' Fatal ERROR in main(): ', errmess)
