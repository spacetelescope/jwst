import os
import numpy as np
from astropy.io import fits
from jwst import datamodels
from . import input_file_types
from . import main

#
# Subarrays:
#  {Instrument:{Name: (colstart, colstop, rowstart, rowstop, Name)}
#  Values are 1-indexed
#
subarrays = {'MIRI': {'FULL': (1, 1032, 1, 1024, 'FULL'),
                     'MASK1065': (1, 256, 1, 256, 'MASK1065'),
                     'MASK1140': (1, 256, 229, 484, 'MASK1140'),
                     'MASK1550': (1, 256, 452, 707, 'MASK1550'),
                     'MASKLYOT': (1, 320, 705, 1024, 'MASKLYOT'),
                     'BRIGHTSKY': (1, 864, 1, 512, 'BRIGHTSKY'),
                     'SUB256': (1, 608, 1, 256, 'SUB256'),
                     'SUB128': (1, 132, 897, 1024, 'SUB128'),
                     'SUB64': (1, 68, 897, 960, 'SUB64'),
                     'SLITLESSPRISM': (1, 68, 321, 1024, 'SLITLESSPRISM')
                     },
             'NIRCAM': {'FULL': (1, 2048, 1, 2048, 'FULL')},
             'NIRSPEC': {'FULL': (1, 2048, 1, 2048, 'FULL'),
                        'ALLSLITS': (1, 2048, 890, 1145, 'ALLSLITS'),
                        'SUBS200A1': (1, 2048, 911, 974, 'SUBS200A1'),
                        'SUBS200A2': (1, 2048, 951, 1014, 'SUBS200A2'),
                        'SUBS200B1': (1, 2048, 1061, 1124, 'SUBS200B1'),
                        'SUBS400A1': (1, 2048, 993, 1056, 'SUBS400A1'),
                        'SUB2048': (1, 2048, 1044, 1075, 'SUB2048')
                        },
             'NIRISS': {'FULL': (1, 2048, 1, 2048, 'FULL'),
                        'SUB256': (1793, 2048, 1793, 2048, 'SUB256'),
                        'SUB128': (1921, 2048, 1921, 2048, 'SUB128'),
                        'SUB64': (1985, 2048, 1985, 2048, 'SUB64'),
                        'AMITA': (1953, 2016, 1985, 2048, 'AMITA'),
                        'SOSSTA': (1202, 1265, 392, 455, 'SOSSTA'),
                        'SOSSFT1': (1793, 2048, 1, 2048, 'SOSSFT1'),
                        'SOSSFT4': (1793, 2048, 1, 2048, 'SOSSFT4'),
                        'SOSSBT1': (1969, 2048, 1, 2048, 'SOSSBT1'),
                        'SOSSBT4': (1969, 2048, 1, 2048, 'SOSSBT4')
                        },
             'TFI': {'FULL': (1, 2048, 1, 2048, 'FULL'),
                     'SUB512': (1, 512, 769, 1280, 'SUB512'),
                     'SUB256': (1, 256, 897, 1152, 'SUB256'),
                     'SUB128': (1, 128, 961, 1088, 'SUB128'),
                     'SUB64': (1, 64, 993, 1056, 'SUB64'),
                     'SUBC050': (1, 308, 1741, 2048, 'SUBC050'),
                     'SUBC075': (513, 820, 1741, 2048, 'SUBC075'),
                     'SUBC150': (1025, 1332, 1741, 2048, 'SUBC150'),
                     'SUBC200': (1537, 1844, 1741, 2048, 'SUBC200'),
                     'SUBT050': (261, 308, 1741, 1788, 'SUBT050'),
                     'SUBT075': (773, 820, 1741, 1788, 'SUBT075'),
                     'SUBT150': (1285, 1332, 1741, 1788, 'SUBT150'),
                     'SUBT200': (1797, 1844, 1741, 1788, 'SUBT200'),
                     'SUBPA': (1, 308, 1025, 1332, 'SUBPA'),
                     'SUBPB': (717, 1024, 1025, 1332, 'SUBPB'),
                     'SUBPC': (1025, 1332, 1025, 1332, 'SUBPC')
                     },
             'GUIDER': {'FULL': (1, 2048, 1, 2048, 'FULL')
                        }
              }

SCAs = {'NRCA1': 481,
        'NRCA2': 482,
        'NRCA3': 483,
        'NRCA4': 484,
        'NRCALONG': 485,
        'NRCB1': 486,
        'NRCB2': 487,
        'NRCB3': 488,
        'NRCB4': 489,
        'NRCBLONG': 490,
        'NRS1': 491,
        'NRS2': 492,
        'MIRIMAGE': 493,
        'MIRIFUSHORT': 494,
        'MIRIFULONG': 495,
        'NIRISS': 496,
        'GUIDER1': 497,
        'GUIDER2': 498}

def rotated_and_flipped(input_hdulist):
    try:
        if input_hdulist[0].header['DMSCOORD'] == 'JWST':
            return True
        else:
            return False
    except KeyError:
        return False

def getSCA_ID(input_hdulist):
    if input_file_types.is_ncont(input_hdulist):
        module = input_hdulist[0].header['MODULE'].strip()
        detector = input_hdulist[0].header['DETECTOR'].strip()
        SCA = input_hdulist[0].header['SCA']
        if detector == 'SW':
            detector_name = 'NRC' + module + str(SCA)
        else:
            detector_name = 'NRC' + module + 'LONG'
        sca_id = SCAs[detector_name]
        return sca_id
    else:
        return int(input_hdulist[0].header['SCA_ID'])

def get_subarray_name(instrument, colstart, colstop, rowstart, rowstop):
    for subarray in subarrays[instrument].values():
        if rowstart == subarray[2] and rowstop == subarray[3] and colstart == subarray[0] \
                and colstop == subarray[1]:
            return subarray[4]
    return 'GENERIC'

def sanitize(header, keyword, type=float, default=0.0):
    """Sometimes the input headers aren't what are wanted in the output,
    so they need to be cleaned up first.  We consider the following cases:

    1. Keyword absent

    In this case, we just return the value of the 'default' parameter

    2. Keyword has wrong type

    The desired type is given by the parameter 'type'.  We try a conversion
    to this type - if it succeeds, we return the converted value, if it fails,
    we return the 'default' parameter

    3. Keyword is OK, has the correct type

    In this case, we return the value of the keyword

    Parameters:

    header: fits Header object
    ------

    The input data header object, usually from the primary header

    keyword: keyword string
    -------

    The keyword that is to be sanitized

    type: desired type
    ----

    This is the required type of the keyword.  Defaults to float

    default: default value
    -------

    The default value returned if the keyword doesn't exist, or is of the wrong type
    and cannot be converted
    """
    try:
        value = header[keyword]
    except KeyError:
        return type(default)
    #
    # try converting the keyword
    try:
        converted_value = type(value)
    except ValueError:
        return type(default)
    return converted_value

def get_nircam_subarray(header):
    """Get subarray parameters from NIRCAM data header

    Parameters
    ----------

    header : fits Header object
        The header of the data whose subarray parameters are to be determined

    Returns
    -------

    detector_row_start, detector_column_start : tuple of ints
        The 1-indexed values of the starting row and column
    """

    #
    # ROWSTART and COLSTART are zero-indexed, ROWCORNR and COLCORNR
    # are 1-indexed
    # Try to get ROWCORNR from header.  If that doesn't work, try ROWSTART
    detector_row_start = None
    try:
        detector_row_start = int(header['ROWCORNR'])
    except KeyError:
        try:
            detector_row_start = int(float(header['ROWSTART'])) + 1
        except KeyError:
            pass
    if detector_row_start is None:
        print('Unable to get subarray ROWSTART, using 1')
        detector_row_start = 1

    #
    # Now try to get COLCORNR from header.  If that doesn't work, try COLSTART
    detector_column_start = None
    try:
        detector_column_start = int(header['COLCORNR'])
    except KeyError:
        try:
            detector_column_start = int(float(header['COLSTART'])) + 1
        except KeyError:
            pass
    if detector_column_start is None:
        print('Unable to get subarray COLSTART, using 1')
        detector_column_start = 1

    return detector_row_start, detector_column_start

def flip_rotate(input_hdulist):
    """Given a FITS HDUList, flip and rotate the data as appropriate.
    Decide on how to flip and rotate by the SCA_ID keyword.
    NIRSpec data always needs to be flipped and rotated: data from SCA 491
    needs to be flipped across the line x=y (achieved using the numpy
    'swapaxes' function), while data from SCA 492 needs to be flipped across
    the other diagonal (achieved using swapaxes then flipping both axes, i.e.
    a flip across x=y then a 180 degree rotation).
    NIRCAM detectors A1, A3, ALONG, B2 and B4 (481, 483, 485, 487 and 489)
    need to be flipped in the X direction.
    NIRCAM detectors A2, A4, B1, B3 and BLONG (482, 484, 486, 488 and 490)
    need to be flipped in the Y direction.
    NIRISS (SCA 496) and FGS1 (#497) need to be flipped across x=y and then
    rotated 180 degrees like NRS2.
    FGS 2 (#498) needs to be flipped across x=y and then flipped along the
    Y direction.
    Returns a FITS HDUList with the data flipped and rotated accordingly"""

    header = input_hdulist[0].header

    sca = getSCA_ID(input_hdulist)
    cube = input_hdulist[0].data

    instrument = header['INSTRUME']

    if sca in [493, 494, 495]:
        # MIRI Imaging, IFUSHORT, IFULONG
        print('Data from SCA %d does not need to be flipped and/or rotated' % sca)
        try:
            rowstart = header['ROWSTART']
            rowstop = header['ROWSTOP']
            #
            # columns are per amplifier, so there are only 258 in the full frame
            colstart = int(4 * (header['COLSTART'] - 1) + 1)
            colstop = int(4 * header['COLSTOP'])
        except KeyError:
            #
            # VM2 data doesn't have any subarray keywords at all, so for
            # now we'll make some up
            rowstart = 1
            rowstop = 1024
            colstart = 1
            colstop = 1032

        ncols = colstop - colstart + 1
        nrows = rowstop - rowstart + 1
        fastaxis = 1
        slowaxis = 2
        rcube = cube
        #
        # Set EXP_TYPE for MIRI LRS data
        if sca == 493:
            if header['EXP_TYPE'] == '':
                miri_filter = header['FWA_POS']
                if miri_filter == 'P750L':
                    subarray_name = get_subarray_name(instrument, colstart,
                                                      colstop, rowstart,
                                                      rowstop)
                    if subarray_name == 'SLITLESSPRISM':
                        header['EXP_TYPE'] = 'MIR_LRS-SLITLESS'
                    else:
                        header['EXP_TYPE'] = 'MIR_LRS-FIXEDSLIT'
        else:
            dgaapos = header['DGAA_POS'].strip()
            dgabpos = header['DGAB_POS'].strip()
            if dgaapos != dgabpos:
                header['DGAA_POS'] = dgaapos + '-' + dgabpos

    elif sca in [481, 483, 485, 487, 489]:
        #
        # NIRCAM A1, A3, ALONG, B2, B4
        # Flip horizontally
        #
        rcube = cube[:, :, ::-1]
        detector_row_start, detector_column_start = get_nircam_subarray(header)

        ncols = header['NAXIS1']
        nrows = header['NAXIS2']
        #
        # Since we're flipping these data in the X-direction only, ROWSTART and ROWSTOP
        # will be unchanged.
        # COLSTART and COLSTOP will swap and subtract from 2049
        # FASTAXIS is -1, as the detector is now read from right to left
        rowstart = detector_row_start
        rowstop = rowstart + nrows - 1
        colstop = 2049 - detector_column_start
        colstart = colstop - ncols + 1
        fastaxis = -1
        slowaxis = 2
    elif sca in [482, 484, 486, 488, 490]:
        #
        # NIRCAM A2, A4, B1, B3, BLONG
        # Flip vertically
        #
        rcube = cube[:, ::-1]
        detector_row_start, detector_column_start = get_nircam_subarray(header)

        ncols = header['NAXIS1']
        nrows = header['NAXIS2']
        #
        # Since we're flipping these data in the Y-direction only, COLSTART and COLSTOP
        # will be unchanged.
        # ROWSTART and ROWSTOP will swap and subtract from 2049
        # FASTAXIS is 1, as the detector is still read from left to right
        colstart = detector_column_start
        colstop = colstart + ncols - 1
        rowstop = 2049 - detector_row_start
        rowstart = rowstop - nrows + 1
        fastaxis = 1
        slowaxis = -2
    elif sca == 491:
        #
        # NIRSpec NRS1
        #
        rcube = np.swapaxes(cube, 1, 2)
        if input_file_types.is_nirspec_ips(input_hdulist):
            #
            # IPS data is all full-frame
            detector_row_start = 0
            detector_column_start = 0
        else:
            try:
                detector_row_start = int(float(header['A1_ROW_C']))
                print(('Detector row start = %d, from keyword A1_ROW_C' % detector_row_start))
            except KeyError:
                try:
                    detector_row_start = int(float(header['A1WINVSA']))
                    print(('Detector row start = %d, from keyword A1WINVSA' % detector_row_start))
                except KeyError:
                    print('Unable to get keyword A1_ROW_C or A1WINVSA, using 1')
                    detector_row_start = 0  # corresponds to 1 in detector/IRAF coordinates
            try:
                detector_column_start = int(float(header['A1_COL_C']))
                print(('Detector column start = %d, from keyword A1_COL_C' % detector_column_start))
            except KeyError:
                try:
                    detector_column_start = int(float(header['A1WINHSA']))
                    print(('Detector column start = %d, from keyword A1WINHSA' % detector_column_start))
                except KeyError:
                    print('Unable to get keyword A1_COL_C or A1WINHSA, using 1')
                    detector_column_start = 0  # corresponds to 1 in detector/IRAF coordinates
        ncols = int(float(header['NAXIS1']))
        nrows = int(float(header['NAXIS2']))
        if input_file_types.is_nirspec_irs2(input_hdulist):
            colstart = 1
            colstop = 2048
            rowstart = 1
            rowstop = 2048
        else:
            #
            #  Interchange X and Y coordinates for detector #491
            #  FASTAXIS is 2, as the detector is now read from bottom to top
            colstart = detector_row_start + 1
            colstop = colstart + nrows - 1
            rowstart = detector_column_start + 1
            rowstop = rowstart + ncols - 1
        ncols, nrows = nrows, ncols
        fastaxis = 2
        slowaxis = 1
        gwaxtilt = sanitize(header, 'GWA_XTIL', type=float)
        gwaytilt = sanitize(header, 'GWA_YTIL', type=float)
        gwa_tilt = sanitize(header, 'GWA_TTIL', type=float)
        header['GWA_XTIL'] = gwaxtilt
        header['GWA_YTIL'] = gwaytilt
        header['GWA_TTIL'] = gwa_tilt
    elif sca == 492:
        #
        # NIRSpec NRS2
        #
        rcube = np.swapaxes(cube, 1, 2)[:, ::-1, ::-1]
        if input_file_types.is_nirspec_ips(input_hdulist):
            #
            # IPS data is all full-frame
            detector_row_start = 0
            detector_column_start = 0
        else:
            try:
                detector_row_start = int(float(header['A2_ROW_C']))
                print(('Detector row start = %d from keyword A2_ROW_C' % detector_row_start))
            except KeyError:
                try:
                    detector_row_start = int(float(header['A1WINVSA']))
                    print(('Detector row start = %d from keyword A1WINVSA' % detector_row_start))
                except KeyError:
                    print('Unable to get keyword A2_ROW_C or A1WINVSA, using 1')
                    detector_row_start = 0  # corresponds to 1 in detector/IRAF coordinates
            try:
                detector_column_start = int(float(header['A2_COL_C']))
                print(('Detector columns start = %d from keyword A2_COL_C' % detector_column_start))
            except KeyError:
                try:
                    detector_column_start = int(float(header['A1WINHSA']))
                    print(('Detector column start = %d from keyword A1WINHSA' % detector_column_start))
                except KeyError:
                    print('Unable to get keyword A2_COL_C or A1WINHSA, using 1')
                    detector_column_start = 0
        ncols = int(float(header['NAXIS1']))
        nrows = int(float(header['NAXIS2']))
        if input_file_types.is_nirspec_irs2(input_hdulist):
            colstart = 1
            colstop = 2048
            rowstart = 1
            rowstop = 2048
        else:
            #
            #  Interchange X and Y coordinates for detector #492 and then subtract from 2049
            #  Fastaxis is now -2, as the detector is read from top to bottom
            det_xmin = detector_column_start + 1
            det_xmax = detector_column_start + ncols
            det_ymin = detector_row_start + 1
            det_ymax = detector_row_start + nrows
            colstart = 2049 - det_ymax
            colstop = 2049 - det_ymin
            rowstart = 2049 - det_xmax
            rowstop = 2049 - det_xmin
        nrows, ncols = ncols, nrows
        fastaxis = -2
        slowaxis = -1
        gwaxtilt = sanitize(header, 'GWA_XTIL', type=float)
        gwaytilt = sanitize(header, 'GWA_YTIL', type=float)
        gwa_tilt = sanitize(header, 'GWA_TTIL', type=float)
        header['GWA_XTIL'] = gwaxtilt
        header['GWA_YTIL'] = gwaytilt
        header['GWA_TTIL'] = gwa_tilt
    elif sca == 496:
        #
        # TFI/NIRISS is like NRS2: flipped across the line X=Y and then
        # rotated 180 degrees
        rcube = np.swapaxes(cube, 1, 2)[:, ::-1, ::-1]
        #
        # TFI and NIRISS data have different keywords
        if input_file_types.is_tfi(input_hdulist):
            #
            #  TFI data doesn't have any subarray keywords at all,
            #  so for now we'll make some up
            rowstart = 1
            rowstop = header['NAXIS2']
            colstart = 1
            colstop = header['NAXIS1']
            #
            # Fix up the headers so we can treat TFI data as NIRISS
#            if instrument == 'TFI':
#                instrument = 'NIRISS'
#                header.update('INSTRUME', 'NIRISS')
#                header.update('DETECTOR', 'NIRISS')
        elif input_file_types.is_niriss(input_hdulist):
            detector_rowstart = header['ROWCORNR']
            detector_rowstop = detector_rowstart + header['NAXIS2'] - 1
            detector_colstart = header['COLCORNR']
            detector_colstop = detector_colstart + header['NAXIS1'] - 1
            #
            # Flip across X=Y and rotate 180 degrees
            rowstart = 2049 - detector_colstop
            rowstop = 2049 - detector_colstart
            colstart = 2049 - detector_rowstop
            colstop = 2049 - detector_rowstart
            #
            # Choose the EXP_TYPE if it's not in the proposal
            if header['EXP_TYPE'] == '':
                if header['PUPIL'].strip() == 'GR700XD':
                    header['EXP_TYPE'] = 'NIS_NOSS'
                elif header['FILTER'].startswith('GR'):
                    header['EXP_TYPE'] = 'NIS_WFSS'
                elif header['PUPIL'] == 'NRM':
                    header['EXP_TYPE'] = 'NIS_AMI'
                else:
                    header['EXP_TYPE'] = 'NIS_IMAGING'
        else:
            print("SCA 496 must be either TFI or NIRISS")
        fastaxis = -2
        slowaxis = -1
        nrows = header['NAXIS2']
        ncols = header['NAXIS1']

    elif sca == 497:
        #
        # GUIDER1 is like NIRISS, equivalent
        # to a flip across X=Y and 180 degree rotation
        rcube = np.swapaxes(cube, 1, 2)[:, ::-1, ::-1]
        detector_rowstart = header['ROWCORNR']
        detector_rowstop = detector_rowstart + header['NAXIS2'] - 1
        detector_colstart = header['COLCORNR']
        detector_colstop = detector_colstart + header['NAXIS1'] - 1
        #
        # Flip across X=Y and rotate 180 degrees
        rowstart = 2049 - detector_colstop
        rowstop = 2049 - detector_colstart
        colstart = 2049 - detector_rowstop
        colstop = 2049 - detector_rowstart
        fastaxis = -2
        slowaxis = -1
        nrows = header['NAXIS2']
        ncols = header['NAXIS1']
#
# DMS spec says INSTRUME should be FGS, while FITSWriter puts GUIDER
        header['INSTRUME'] = 'FGS'
    elif sca == 498:
        #
        # GUIDER2 is rotated 90 degrees anticlockwise
        rcube = np.rot90(cube, 1, (2,1))
        detector_rowstart = header['ACROWCOR']
        detector_rowstop = detector_rowstart + header['NAXIS2'] - 1
        detector_colstart = header['ACCOLCOR']
        detector_colstop = detector_colstart + header['NAXIS1'] - 1
        #
        colstart = 2049 - detector_rowstop
        colstop = 2049 - detector_rowstart
        rowstart = detector_colstart
        rowstop = detector_colstop
        fastaxis = 2
        slowaxis = -1
        nrows = header['NAXIS2']
        ncols = header['NAXIS1']
#
# DMS spec says INSTRUME should be FGS, while FITSWriter puts GUIDER
        header['INSTRUME'] = 'FGS'
    del cube
    header['DMSCOORD'] = 'JWST'
    header['COLSTART'] = colstart
    header['COLSTOP'] = colstop
    header['ROWSTART'] = rowstart
    header['ROWSTOP'] = rowstop
    header['NAXIS1'] = ncols
    header['NAXIS2'] = nrows
    header['FASTAXIS'] = fastaxis
    header['SLOWAXIS'] = slowaxis
    #
    #  Update the subarray parameters
    subarray_name = get_subarray_name(instrument, colstart, colstop, rowstart, rowstop)
    header['SUBARRAY'] = subarray_name
    input_hdulist[0].data = rcube
    return

def create_subarrays(input_hdulist, subarrays):
    """Given a full-frame image, create subarrays from it.  Returns a list
    of filenames corresponding to the subarrays created"""
    subarray_files = []
    #
    #  Extract the data
    #
    for subarray in subarrays:
        outputfilename = create_single_subarray(input_hdulist, subarray)
        subarray_files.append(outputfilename)
    #
    #  All done, now we return the names of the subarray files
    #
    return subarray_files

def create_single_subarray(input_hdulist, subarray):
    full = {'MIRI': (1032, 1280),
            'NIRCAM': (2048, 2048),
            'NIRSPEC': (2048, 2048),
            'NIRISS': (2048, 2048),
            'TFI': (2048, 2048),
            'FGS': (2048, 2048)
            }
    #
    #  Make sure we are extracting from a full-frame exposure
    #
    header = input_hdulist[0].header
    cube = input_hdulist[0].data
    instrument = header['INSTRUME']
    (nx, ny) = (cube.shape[2], cube.shape[1])
    if nx != full[instrument][0] and ny != full[instrument][1]:
        print("Can only take subarrays from full-frame exposures")
        return None
    print('Subarray = %s' % subarray)
    print(subarrays[instrument][subarray])
    xstart, xstop, ystart, ystop, name = subarrays[instrument][subarray]
    subarray_data = cube[:, ystart - 1:ystop, xstart - 1:xstop]
    if (instrument == 'MIRI'):
        #
        #  Add the Reference Output if MIRI
        #
        # Start by reforming Reference output to be contiguous
        #
        refout = cube[:, 1024:, :]
        continuous_refout = refout.reshape(refout.shape[0], 1024, 258)
        print('Reference Output subarray in contiguous goes from')
        print('x:  %d to %d' % ((xstart - 1) // 4, xstop // 4))
        print('y:  %d to %d' % ((ystart - 1), ystop))
        refout_subarray = continuous_refout[:, (ystart - 1):ystop,
                                          (xstart - 1) // 4:(xstop + 1) // 4]
##        print('Reshaped array has dimensions %d by %d' % \
##              (refout_subarray.shape[2]*4, refout_subarray.shape[1]/4))
##        reshaped_refout_subarray = refout_subarray.reshape(refout_subarray.shape[0],
##                                                        refout_subarray.shape[1]/4,
##                                                        refout_subarray.shape[2]*4)
##        out_shape = (cube.shape[0], ((ystop-ystart+1)/4)*5, xstop-xstart+1)
##        out_data = np.zeros(out_shape)
##        print("Shape of output array is ", out_data.shape)
##        out_data[:, :(ystop-ystart+1),:(xstop-xstart+1)] = subarray_data
##        out_data[:, (ystop-ystart+1):, :] = reshaped_refout_subarray
##    else:
    out_data = subarray_data
    #
    #  Now write the subarray out to file
    outputfilename = (''.join((name, '.fits'))).lower()
    #
##     #  Get rid of anything left behind by the last run
##     try:
##         os.remove(outputfilename)
##     except:
##         pass
    output_hdulist = fits.HDUList()
    output_primaryhdu = fits.PrimaryHDU()
    output_primary_header = header.copy()
    #
    #  Add the keywords to the header
    #    output_primary_header.update('SUBARRAY', name)
    output_primary_header['ROWSTART'] = ystart
    output_primary_header['ROWSTOP'] = ystop
    output_primary_header['COLSTART'] = xstart
    output_primary_header['COLSTOP'] = xstop
    output_primary_header['SUBARRAY'] = name
    output_primaryhdu.header = output_primary_header
    output_primaryhdu.data = out_data
    output_hdulist.append(output_primaryhdu)
#    try:
#        output_hdulist.append(input_hdulist[1].copy())
#    except:
#        pass
    if instrument == 'MIRI':
        output_refhdu = fits.ImageHDU()
        output_refhdu.data = refout_subarray
        output_hdulist.append(output_refhdu)

    output_hdulist.writeto(outputfilename)
    output_hdulist.close()
    return outputfilename

def split_data_and_refout(hdulist):

    hdr = hdulist[0].header
    #
    # ncols and nrows refer to the science data dimensions, not the keywords
    # in the raw data - they were fixed up in flip_rotate
    ncols = hdr['COLSTOP'] - hdr['COLSTART'] + 1
    nrows = hdr['ROWSTOP'] - hdr['ROWSTART'] + 1
    #
    # make sure they're integers with a nearest integer calculation
    ncols = int(ncols + 0.5)
    nrows = int(nrows + 0.5)
    fulldata = hdulist[0].data
    detectordata = fulldata[:, :nrows]
    refoutdata = fulldata[:, nrows:]
    #
    # Reference output has 1/4 the columns of science data
    refout = refoutdata.reshape((fulldata.shape[0], nrows, ncols // 4))
    hdulist[0].data = detectordata
    try:
        del hdulist[1]
    except:
        pass
    refhdu = fits.ImageHDU()
    refhdu.data = refout
    hdulist.append(refhdu)
    return

def create_dms(base_file, level="1b", parfile=None, subarray=None, exp_type='UNKNOWN'):
    #
    # Create a Level 1b file from a full-frame base file
    # base_file is the full-format image the subarray is created from
    # subarray is a tuple describing the subarray:
    # (nxstart, nxstop, nystart, nystop, name)
    #

    base_root = base_file.split('.')[0]

    base_hdulist = fits.open(base_file)

    #
    # Add the exp_type to the input hdulist
    base_hdulist[0].header['EXP_TYPE'] = exp_type

    #
    #  If we've already rotated/flipped the file, it will have a keyword set
    #  and we skip this step
    #
    if not rotated_and_flipped(base_hdulist):
        print('Rotating and flipping data %s' % (base_file))
        flip_rotate(base_hdulist)

    instrument = base_hdulist[0].header['INSTRUME']

    if subarray is not None:
        subarray_file = create_single_subarray(base_hdulist,
                                               subarray)
        print('Parfile = %s' % parfile)
        if parfile is None:

            parfile = ''.join((base_root, '_', subarray[4].lower(),
                               '_observation_identifiers.dat'))
            print(parfile)

        hdulist = main.generate([subarray_file, parfile], level=level)
        filename = main.guess_filename(hdulist)
        hdulist[0].header['FILENAME'] = filename
        hdulist.writeto(filename, output_verify='silentfix')
        os.remove(subarray_file)
    else:
        if instrument == 'MIRI':
            split_data_and_refout(base_hdulist)
        hdulist = main.generate([base_hdulist, parfile],
                                          level=level)
        filename = main.guess_filename(hdulist)
        hdulist[0].header['FILENAME'] = filename
        hdulist.writeto(filename, output_verify='silentfix')

    base_hdulist.close()
    #
    # Try and open the file as a JWST datamodel
    a = datamodels.open(filename)
    print("%s opened as a %s" % (filename, a.__module__))
    a.close()
    return
