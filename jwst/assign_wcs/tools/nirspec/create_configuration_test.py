from astropy.io import fits
import numpy as np
import numpy.lib.recfunctions

def add_configuration_records(filename):
    # filename = 'jw95065006001_01_msa.fits'
    # Read in the base fits file.
    msa_conf = fits.open(filename)

    # Add in the different configuration test conditions

    rectype = np.dtype([('slitlet_id', '>i2'),
                              ('msa_metadata_id', '>i2'),
                              ('shutter_quadrant', '>i2'),
                              ('shutter_row', '>i2'),
                              ('shutter_column', '>i2'),
                              ('source_id', '>i2'),
                              ('background', 'S1'),
                              ('shutter_state', 'S6'),
                              ('estimated_source_in_shutter_x', '>f4'),
                              ('estimated_source_in_shutter_y', '>f4')])

    # The base one is
    base = np.array([(12, 2, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 2, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 2, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 64, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 64, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 64, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 65, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 65, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 65, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 66, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 66, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 66, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 67, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 67, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 67, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 68, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 68, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 68, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 69, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 69, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 69, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 70, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 70, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 70, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 71, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
               (12, 71, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               (12, 71, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan)],\
              rectype)

    # Test 1: Kinda normal
    test1 = np.array([
                (55, 12, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
                (55, 12, 4, 251, 23, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (55, 12, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
                (55, 12, 4, 251, 25, 1, 'Y', 'OPEN', np.nan, np.nan),
                (55, 12, 4, 251, 26, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               ],\
              rectype)


    # Test 2: Create slitlet_id set with no background open
    #         This should fail as there should be only one "N"
    test2 = np.array([
                (56, 13, 4, 251, 22, 1, 'N', 'OPEN', np.nan, np.nan),
                (56, 13, 4, 251, 23, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
                (56, 13, 4, 251, 24, 1, 'Y', 'OPEN', np.nan, np.nan),
                (56, 13, 4, 251, 25, 1, 'Y', 'OPEN', np.nan, np.nan),
                (56, 13, 4, 251, 26, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
               ],\
              rectype)

    # Test 3: All background
    test3 = np.array([
                (57, 14, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
                (57, 14, 4, 251, 23, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (57, 14, 4, 251, 24, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
               ],\
              rectype)

    # Test 4: Empty in between
    test4 = np.array([
                (58, 15, 4, 251, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
                (58, 15, 4, 251, 23, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (58, 15, 4, 251, 24, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
                (58, 15, 4, 251, 25, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (58, 15, 4, 251, 27, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (58, 15, 4, 251, 28, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
               ],\
              rectype)

    # Test 5: Empty in between
    test5 = np.array([
                (59, 16, 4, 256, 22, 1, 'Y', 'OPEN', np.nan, np.nan),
                (59, 16, 4, 256, 23, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (59, 16, 4, 256, 24, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
                (59, 16, 4, 256, 25, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (59, 16, 4, 256, 27, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (59, 16, 4, 256, 28, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (60, 16, 4, 258, 30, 1, 'Y', 'OPEN', np.nan, np.nan),
                (60, 16, 4, 258, 31, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (60, 16, 4, 258, 32, 1, 'N', 'OPEN', 0.18283921, 0.31907734),
                (60, 16, 4, 258, 33, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (60, 16, 4, 258, 34, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
                (60, 16, 4, 258, 35, 1, 'Y', 'OPEN', 0.18283921, 0.31907734),
               ],\
              rectype)


    # And now write out the msa_conf to a new file
    tt = numpy.lib.recfunctions.stack_arrays((base, test1, test2, test3, test4, test5), usemask=False)
    msa_conf[2].data = tt
    msa_conf.writeto('test_configuration_msa.fits')

if __name__ == '__main__':
    print('Call add_configuration_records(filename) directly.')
