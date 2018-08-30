# Copyright (C) 2010-2011 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

# THIRD-PARTY
from astropy.io import fits as pyfits

# STDLIB
import glob
import os
from io import StringIO

# LOCAL
from . import input_file_types
from . import objects as O
from . import template
from . import util

def select_type(fitsfile, level):
    """
    Auto-selects the appropriate output filetype based on an input
    file.

    Parameters
    ----------
    fitsfile : string or `pyfits.HDUList` instance
        Input FITS file to use to determine output file type.

    level : string
        The level of FITS file to generate.  May be '1a' or '1b'.
    """
    if level not in ('1a', '1b', '2a'):
        raise ValueError("level must be 1a, 1b or 2a")

    if level == '1a':
        return 'level1a'
    elif level == '1b' or level == '2a':
        fitsfile, hdulist, opened = util.get_fitsfile(fitsfile)

        if input_file_types.is_fitswriter(hdulist):
            if input_file_types.is_miri(hdulist):
                if input_file_types.is_miri_ifu(hdulist):
                    si = 'miri_ifu'
                elif input_file_types.is_miri_lrs(hdulist):
                    si = 'miri_lrs'
                else:
                    si = 'miri_imaging'
            elif input_file_types.is_nircam_fm1(hdulist):
                si = 'nircam_fm1'
            elif input_file_types.is_nirspec_fm1(hdulist):
                si = 'nirspec_fm1'
            elif input_file_types.is_tfi(hdulist):
                si = 'tfi'
            elif input_file_types.is_niriss(hdulist):
                if input_file_types.is_niriss_spec_horiz(hdulist):
                    si = 'niriss_spec_horiz'
                elif input_file_types.is_niriss_spec_vert(hdulist):
                    si = 'niriss_spec_vert'
                else:
                    si = 'niriss'
            else:
                sca_to_si_mapping = {
                    0x1e1: 'nircam',
                    0x1e2: 'nircam',
                    0x1e3: 'nircam',
                    0x1e4: 'nircam',
                    0x1e5: 'nircam',
                    0x1e6: 'nircam',
                    0x1e7: 'nircam',
                    0x1e8: 'nircam',
                    0x1e9: 'nircam',
                    0x1ea: 'nircam',
                    0x1eb: 'nirspec',
                    0x1ec: 'nirspec',
                    0x1ed: 'miri',
                    0x1ee: 'miri',
                    0x1ef: 'miri',
                    0x1f0: 'niriss',
                    0x1f1: 'guider',
                    0x1f2: 'guider'
                    }

                try:
                    sca_id = int(hdulist[0].header['SCA_ID'])
                    try:
                        si = sca_to_si_mapping[sca_id]
                    except KeyError:
                        raise ValueError(
                            "Can not automatically determine filetype from file")
                except:
                    raise ValueError("No SCA_ID keyword present")
        elif input_file_types.is_ncont(hdulist):
            si = 'nircam_ncont'
        elif input_file_types.is_swts(hdulist):
            si = 'nircam_swts'
        elif input_file_types.is_nirspec_ips(hdulist):
            si = 'nirspec_ips'
        else:
            raise ValueError("Can not automatically determine filetype from file")
        if opened:
            hdulist.close()
        template = 'level%s_%s' % (level, si)
        print('Template: %s' % template)
        return template


def verify(fitsfile, filetype, error_collector=None):
    """
    Verifies that a FITS file matches a particular filetype
    definition.

    Parameters
    ----------
    fitsfile : string or `pyfits.HDUList` instance
        The FITS file to verify.  May be a path or a `pyfits.HDUList`
        object.

    filetype : string
        Type of FITS file to verify.  To get a list of supported file
        types, use `list_filetypes`.

    error_collector : callable
        A function to collect errors and warnings.  See
        `error_collector functions`_.  If none is provided, a default
        collector will write errors and warnings to stderr.

    Returns
    -------
    success : boolean
        Returns `True` if verification passes 100%, `False` otherwise.
    """
    if filetype is None:
        filetype = select_type(fitsfile, '1b')
    filetype_object = util.get_filetype(filetype)
    error_collector = util.get_error_collector(error_collector)
    fitsfile, hdulist, opened = util.get_fitsfile(fitsfile)

    try:
        state = O.ParseState()
        state.file = fitsfile
        result = filetype_object.verify(hdulist, error_collector, state)
    finally:
        if opened:
            hdulist.close()

    return result


def generate(input_files, filetype=None, level='1b', error_collector=None, verify=True):
    """
    Converts a FITS file to one that matches a particular filetype
    definition.

    Parameters
    ----------
    input_files : list
        This is a list containing the input files.  This must contain
        exactly 1 FITS file, and any number of data source files.  The
        FITS file will be accessible through the ``input`` function in
        the generator template.  The data files will be accessible
        through a function *name* where *name* is the name of the data
        file as defined in its ``<<file>>`` line.

    filetype : string
        The type of FITS file to generate.  To get a list of supported
        file types, use `list_filetypes`.  If no filetype is
        specified, the result will be a level 1B FITS file of the
        appropriate instrument type, determined automatically from the
        input FITS file.

    error_collector : callable
        A function to collect errors and warnings.  See
        `error_collector functions`_.  If none is provided, a default
        collector will write errors and warnings to stderr.

    verify : bool
        If `True` (default) verify the output FITS file.

    Returns
    -------
    hdulist : `pyfits.HDUList`
      The content of the output FITS file.
    """
    data_files = input_file_types.get_inputfiles(input_files)

    if filetype is None:
        filetype = select_type(data_files['input'], level)
    filetype_object = template.get_generator(filetype)
    error_collector = util.get_error_collector(error_collector)

    # TODO: Maybe this shouldn't be hardcoded but specified in the
    # file somehow
    # if 'target' not in data_files or 'program' not in data_files:
    #     raise ValueError(
    #         "At least a target and program data file must be specified")

    try:
        state = O.ParseState()
        output_hdulist = filetype_object.generate(
            data_files, error_collector, state)
    finally:
        for data_file in data_files.values():
            data_file.close()

    if verify:
        state.file = 'RESULT'
        filetype_object.verify(output_hdulist, error_collector, state)

    return output_hdulist


def describe(filetype, format='rst', output=None):
    """
    Generates a document describing the structure of a given FITS
    file.

    Parameters
    ----------
    filetype : string
        The type of FITS file to generate.  To get a list of supported
        file types, use `list_filetypes`.

    format : string
        Must be one of the following:

        - 'rst': reStructuredText

        - 'html': HTML (requires docutils to be installed)

    output : file-like object
        The content will be written to the given file-like object,
        otherwise the content will be returned as a string.
    """
    class Unclosable(StringIO):
        def close(self):
            pass

    filetype_object = util.get_filetype(filetype)

    state = O.ParseState()
    state.file = "N/A"

    if output is None or format != 'rst':
        stream = StringIO()
    else:
        stream = output

    filetype_object.describe(stream, state)

    if format == 'rst':
        if output is None:
            return stream.getvalue()
    else:
        try:
            import docutils
        except ImportError:
            raise ImportError("docutils is required for HTML output")

        from docutils.core import publish_file

        if output is None:
            html_stream = Unclosable()
        else:
            html_stream = output

        stream.seek(0)
        publish_file(source=stream, source_path='<tmp>',
                     destination=html_stream, destination_path='<tmp>',
                     writer_name='html')

        if output is None:
            return html_stream.getvalue()


def list_filetypes():
    """
    Return a list of all of the available filetypes.
    """
    types_dir = util.get_filetypes_dir()

    filetypes = []
    for filename in glob.glob(os.path.join(types_dir, '*.gen.txt')):
        filetype = os.path.basename(filename)[:-8]
        filetypes.append(filetype)

    filetypes.sort()

    return filetypes


def guess_filename(hdulist):
    """
    Based on the given HDUList, guess at a standards-compliant output
    file name.
    """
    # TODO: Or, we could just get this from the FILENAME keyword
    header = hdulist[0].header
    try:
        _program = header['PROGRAM']
    except KeyError:
        print("Error getting PROGRAM from header")
    try:
        _observtn = header['OBSERVTN']
    except KeyError:
        print("Error getting OBSERVTN from header")
    try:
        _visit = header['VISIT']
    except KeyError:
        print("Error getting VISIT from header")
    try:
        _visitgrp = header['VISITGRP']
    except KeyError:
        print("Error getting VISITGRP from header")
    try:
        _seq_id = header['SEQ_ID']
    except KeyError:
        print("Error getting SEQ_ID from header")
    try:
        _act_id = header['ACT_ID']
    except KeyError:
        print("Error getting ACT_ID from header")
    try:
        _exposure = header['EXPOSURE']
    except KeyError:
        print("Error getting EXPOSURE from header")
    try:
        _detector = header['DETECTOR']
    except KeyError:
        print("Error getting DETECTOR from header")
    try:
        #
        #  Latest version (Revision A of JWST-STScI-0021111) has no underscores between
        #  program and observations, nor between observation and visit
        #
        filename = ('jw%s%s%s_%s%s%s_%s_%s_uncal.fits' % (_program,
                                                          _observtn,
                                                          _visit,
                                                          _visitgrp,
                                                          _seq_id,
                                                          _act_id,
                                                          _exposure,
                                                          _detector))
    except KeyError:
        raise ValueError("Could not automatically generate output filename")
    return filename

