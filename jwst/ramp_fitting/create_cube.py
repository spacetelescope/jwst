#! /usr/bin/env python
#
# Author: Dave Grumm
# Program: create_cube.py
# Purpose: routine to create data cube for a specified readout mode for JWST data
# History: 03/24/10 - first version
#          08/17/10 - changed datatype to float64 for compatibility with other programs
#
# The 9 NIRCAM readout modes are:
#  DEEP8, DEEP2, MEDIUM8, MEDIUM2, SHALLOW4, SHALLOW2, BRIGHT2, BRIGHT1, RAPID

# Required Input:
#   'readmode' : read mode
#  'nread' : number of reads
#  'asize' : length of input array
#  'sample_image' : file containing image - must have the same length as asize
#
# Optional Input:
#  'noisy' : whether or not the data cube will be noisy; default = False
#  'verb' : level of verbosity for print statements; verb=2 is extremely verbose,
#           verb=1 is moderate, verb=0 is minimal; default = 0
#
# array output to files:
#    'acube.fits' : the created cube whose dimensions are n_ind_reads * asize * asize,
#                   where n_ind_reads = nframes*nread + (nread-1)*nskip
#                   is the number of individual reads and asize is the array
#                   length input on the command line. The values of nframes.
#                   ngroup, and nskip are determined from the specified readmode
#
# linux usage example:
#  ./create_cube.py 'deep8' 4 100 sim_100.fits True 1
#  ... which specifies the values:
#  readmode = 'deep8', nread = 4, array size = 100, image to sample = 'sim_100.fits', noisy = True, verbosity = 1

import random
import sys
import time

from astropy.io import fits
import numpy as np

ERROR_RETURN = 2
DELTA_T = 10.6 #  time per single readout in seconds
ELECTRON_PER_ADU = 1.3 # for NIRCAM, is 2.2 for MIRI
READ_NOISE = 10. # electrons

class create_cube:
    """ Create a noisy or noiseless data cube.
    """
    def __init__(self, mode, nread, asize, sample_image, noisy=False, verb=0):

        mode = mode.upper()

        # check input mode and set parameters accordingly
        if (mode == 'DEEP8'): ngroup = 20; nframes = 8; nskip = 12; inst = 'NIRCAM'
        elif (mode == 'DEEP2'): ngroup = 20; nframes = 2; nskip = 18; inst = 'NIRCAM'
        elif (mode == 'MEDIUM8'): ngroup = 10; nframes = 8; nskip = 2; inst = 'NIRCAM'
        elif (mode == 'MEDIUM2'): ngroup = 10; nframes = 2; nskip = 8; inst = 'NIRCAM'
        elif (mode == 'SHALLOW4'): ngroup = 10; nframes = 4; nskip = 1; inst = 'NIRCAM'
        elif (mode == 'BRIGHT2'): ngroup = 10; nframes = 2; nskip = 1; inst = 'NIRCAM'
        elif (mode == 'BRIGHT1'): ngroup = 10; nframes = 1; nskip = 1; inst = 'NIRCAM'
        elif (mode == 'RAPID'): ngroup = 10; nframes = 1; nskip = 0; inst = 'NIRCAM'
        elif (mode == 'NRSSLOW'): ngroup = 0; nframes = 4; nskip = 0; inst = 'NIRSPEC'
        elif (mode == 'NRSFAST'): ngroup = 0; nframes = 1; nskip = 0; inst = 'NIRSPEC'
        elif (mode == 'TFISLOW'): ngroup = 0; nframes = 4; nskip = 0; inst = 'TFI'
        elif (mode == 'TFIFAST'): ngroup = 0; nframes = 1; nskip = 0; inst = 'TFI'
        elif (mode == 'MIRISLOW'): ngroup = 10; nframes = 1; nskip = 0; inst = 'MIRI'
        elif (mode == 'MIRIFAST'): ngroup = 0; nframes = 1; nskip = 0; inst = 'MIRI'
        else:
            print('Fatal ERROR:  unsupported mode  ', mode)
            sys.exit(ERROR_RETURN)

        print('The specified input parameters are :')
        print('  inst : ', inst)
        print('  mode : ', mode)
        print('  ngroup : ', ngroup)
        print('  nframes : ', nframes)
        print('  nskip : ', nskip)
        print('  nread : ', nread)
        print('  noisy : ', noisy)
        print('  sample_image : ', sample_image)
        print('  verb = ', verb)

        # do some parameter type checking
        [mode, nread] = check_pars(mode, nread, asize, sample_image, ngroup, noisy)

        self._nread = int(nread)
        self._verb = verb
        self._inst = inst
        self._mode = mode
        self._ngroup = ngroup
        self._nframes = nframes
        self._nskip = nskip
        self._read_range = np.arange(ngroup) + 1
        self._noisy = noisy
        self._read_noise = READ_NOISE / ELECTRON_PER_ADU
        self._asize_1 = int(asize); self._asize_2 = int(asize)
        self._sample_image = sample_image

    def make_cube(self):
        """  add stuff.......................
        """
        inst = self._inst
        mode = self._mode
        ngroup = self._ngroup
        nframes = self._nframes
        nskip = self._nskip
        nread = self._nread
        noisy = self._noisy
        verb = int(self._verb)
        asize_1 = self._asize_1; asize_2 = self._asize_2
        sample_image = self._sample_image

        if (verb > 0): print('Start of make_cube... ')

        print('The parameters to be used are : ')
        print(' inst : ', inst)
        print(' mode : ', mode)
        print(' ngroup : ', ngroup)
        print(' nframes : ', nframes)
        print(' nskip : ', nskip)
        print(' nread : ', nread)
        print(' noisy : ', noisy)
        print(' sample_image : ', sample_image)
        print(' verb = ', verb)
        print('  ')
        print('The readnoise for each pixel is  ', self._read_noise)
    #   print 'The requested readout mode for ', inst, 'is ', mode,' which has ngroup = ', ngroup,','
    #   print '    nframes = ', nframes,', skip = ' ,  nskip, ', nread = ' , nread

        if noisy == 'True':
            print(' The datacube create will be noisy ')
        else:
            print(' The datacube create will be noiseless')

        total_reads = self._nframes

    # calculate the total integration time
        n_ind_reads = nframes * nread + (nread - 1) * nskip
        t_int = DELTA_T * n_ind_reads

        print('The number of reads per group is nframes = ', total_reads,
            ',  the total number of individual reads = ', n_ind_reads,
            'and the total integration time = ', t_int)

    # calculate the integration time per single read (total divided by nread)
        t_read = t_int / float(nread)
        print(' The integration time per single read = ', t_read)

###        acube = np.zeros((n_ind_reads, asize_2, asize_1), dtype = np.float32)  # < 080210
        acube = np.zeros((n_ind_reads, asize_2, asize_1), dtype=np.float64)  # try 080210

        print(' The output cube will have dimensions:', n_ind_reads, asize_2, asize_1)

        if sample_image is None:
            print('Fatal ERROR: no sample_image has been specified ', end=' ')
            sys.exit(ERROR_RETURN)
        else:  # open and read image
            fh_reset = fits.open(sample_image)
            reset_data = fh_reset[0].data
            print(' The shape of the sample_image : ', reset_data.shape)

        sum_abs_total_noise = 0.0

        for ii_slice in range(0, n_ind_reads):   # 1st create noiseless cube, looping over all slices
            acube[ii_slice, :, :] = reset_data * (ii_slice + 1)

            if (verb > 0): print('for slice = ', ii_slice, ' noiseless cube has',
                'min, mean, max, std = ', acube[ii_slice, :, :].min(),
                acube[ii_slice, :, :].mean(), acube[ii_slice, :, :].max(),
                acube[ii_slice, :, :].std())


        if (noisy == 'True'):  # add noise if requested
            for xx_out in range(0, asize_1):
                if ((verb > 0) and ((xx_out / 20.) == int(xx_out / 20.))): print(' xx_out = ', xx_out)

                for yy_out in range(0, asize_2):
                    if (verb > 1): print(' This pixel has xx, yy = ', xx_out, yy_out)
                    for ii_slice in range(0, n_ind_reads):   # over all slices
                        ran_lim = reset_data[yy_out, xx_out]
                        poiss_ran = random.gauss(-np.sqrt(ran_lim), np.sqrt(ran_lim))
                        total_noise = np.sqrt(poiss_ran * poiss_ran + self._read_noise * self._read_noise)

                        acube[ii_slice, yy_out, xx_out] += total_noise
                        sum_abs_total_noise += abs(total_noise)

                        if (verb > 1):
                            print(' total_noise =', total_noise, ' self._read_noise = ', self._read_noise)
                            print('  acube[', ii_slice, yy_out, xx_out, '] = ', acube[ii_slice, yy_out, xx_out])

        else:  # noiseless
            pass

        for ii_slice in range(0, n_ind_reads):   # for diagnostics only
            print('for ii_slice = ', ii_slice, ' acube now has min, mean, max, std = ',
                acube[ii_slice, :, :].min(), acube[ii_slice, :, :].mean(),
                acube[ii_slice, :, :].max(), acube[ii_slice, :, :].std())


        self.write_file(acube, 'acube.fits', n_ind_reads)

        print('The cube in the output file acube.fits has size =   ', acube.shape)


    def write_file(self, data, output_fname, n_ind_reads):
        """ Write data cube to file, with relevant keywords

        @param data: data cube created
        @type data: floats
        @param output_fname: name of output file
        @type output_fname: string
        @param n_ind_reads: number of individual reads
        @type n_ind_reads: int
        """
        fitsobj = fits.HDUList()
        hdu = fits.PrimaryHDU()

        prihdr = hdu.header
        prihdr['NGROUP'] = (self._ngroup, 'number of groups')
        prihdr['NREAD'] = (self._nread, 'number of reads')
        prihdr['READMODE'] = (self._mode, 'read mode')
        prihdr['INSTRU'] = (self._inst, 'instrument')
        prihdr['NSKIP'] = (self._nskip, 'number of frames skipped')
        prihdr['NFRAMES'] = (self._nframes, 'number of frames')
        prihdr['NOISE'] = (self._noisy, 'noise added?')
        prihdr['NINDREAD'] = (n_ind_reads, 'number of individual reads')
        prihdr['NAVGIMAG'] = (int(n_ind_reads / self._nread), 'number of average images')

        hdu.data = data
        fitsobj.append(hdu)
        fitsobj.writeto(output_fname, overwrite=True)
        fitsobj.close()

        print(' Wrote to output_fname = ', output_fname)

def check_pars(mode, nread, asize, sample_image, ngroup, noisy):
    """ Verify that the input values are valid.

    @param mode: mode
    @type mode: string
    @param nread: number of reads
    @type nread: int
    @param asize: size of 1 edge of input sample_image array
    @type asize: int
    @param sample_image: image to sample
    @type sample_image: string
    @param ngroup: maximum number of reads
    @type ngroup: int
    @param noisy: if noise is to be added
    @type noisy: string that is either True or False
    @return: mode, nread
    @rtype: string, int
    """

    if mode == 'INVALID':
        print('Fatal ERROR: unsupported mode: ', mode)
        sys.exit(ERROR_RETURN)

    read_range = np.arange(ngroup) + 1
    if int(nread) not in read_range:
        print(' Fatal ERROR: requested read value ', nread, ' is not an allowed value. Try again.')
        sys.exit(ERROR_RETURN)

    if (asize < 0):
        print(' Negative value of asize; asize must be positive ')
        sys.exit(ERROR_RETURN)

    if ((noisy != 'True') and (noisy != 'False')):
        print(' WARNING : noisy is neither True or False so setting noisy to False')
        noisy = 'False'

    try:
        fh_reset = fits.open(sample_image)
        fh_reset[0].data
    except Exception:
        print('Fatal ERROR: unable to access data from sample_image ', sample_image)
        sys.exit(ERROR_RETURN)

    return mode, int(nread)

if __name__ == "__main__":
    """ Get mode (and other arguments?), and call make_cube.
          @param cmdline: command-line arguments
          @type cmdline: list of strings
      """
    usage = "usage:  ./create_cube mode nread asize sample_image [noisy [verbosity]]"

    if len(sys.argv) < 5:
        print("usage:  ./create_cube mode nread asize sample_image [noisy [verbosity]]")
        sys.exit(ERROR_RETURN)

    mode = sys.argv[1]
    nread = sys.argv[2]
    asize = sys.argv[3]
    sample_image = sys.argv[4]

    if (len(sys.argv) > 5):
        noisy = sys.argv[5]
    else:
        noisy = False

    if (len(sys.argv) > 6):
        verb = sys.argv[6]
    else:
        verb = 0

    try:
        tstart0 = time.time()
        print('Starting the program create_cube : ')
        print('The current time (start) = ', time.asctime())
        aa = create_cube(mode=mode, nread=nread, asize=asize, sample_image=sample_image, noisy=noisy, verb=verb)
        status = aa.make_cube()
        print('The current time (end) = ', time.asctime())
        tstop = time.time()
        print(' elapsed time = tstop-tstart0 = ', tstop - tstart0)
    except Exception as errmess:
        print('Fatal ERROR: ', errmess)
        sys.exit(ERROR_RETURN)
