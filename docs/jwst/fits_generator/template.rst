.. _template:

Template file format
====================

File types are described using a simple file format that vaguely
resembles FITS headers.

Since it is necessary to create templates for several different
flavors of data (FITSWriter, NIRSpec simulations, NIRCam homebrew etc)
as well as different EXP_TYPEs that share many sections of data header
but differ in other sections, the templates are divided into sections
that are included.  So a typical template for a particular flavor of
data might look like this::

    <<file nirspec_ifu_level1b>>
    <<header primary>>
    #include "level1b.gen.inc"
    #include 'observation_identifiers.gen.inc'
    #include 'exposure_parameters.gen.inc'
    #include 'program_information.gen.inc'
    #include 'observation_information.gen.inc'
    #include 'visit_information.gen.inc'
    #include 'exposure_information.gen.inc'
    #include 'target_information.gen.inc'
    #include 'exposure_times.gen.inc'
    #include 'exposure_time_parameters.gen.inc'
    #include 'subarray_parameters.gen.inc'
    #include 'nirspec_configuration.gen.inc'
    #include 'lamp_configuration.gen.inc'
    #include 'guide_star_information.gen.inc'
    #include 'jwst_ephemeris_information.gen.inc'
    #include 'spacecraft_pointing_information.gen.inc'
    #include 'aperture_pointing_information.gen.inc'
    #include 'wcs_parameters.gen.inc'
    #include 'velocity_aberration_correction.gen.inc'
    #include 'nirspec_ifu_dither_pattern.gen.inc'
    #include 'time_related.gen.inc'
    
    <<data>>
    
    <<header science>>
    #include 'level1b_sci_extension_basic.gen.inc'
    
    <<data>>
    input[0].data.reshape((input[0].header['NINT'], \
                           input[0].header['NGROUP'], \
                           input[0].header['NAXIS2'], \
                           input[0].header['NAXIS1'])). \
                           astype('uint16')
    
    <<header error>>
    EXTNAME = 'ERR'
    
    <<data>>
    np.ones((input[0].header['NINT'], \
             input[0].header['NGROUP'], \
             input[0].header['NAXIS2'], \
             input[0].header['NAXIS1'])). \
             astype('float32')
    
    <<header data_quality>>
    EXTNAME = "DQ"
    
    <<data>>
    np.zeros((input[0].header['NINT'], \
              input[0].header['NGROUP'], \
              input[0].header['NAXIS2'], \
              input[0].header['NAXIS1']), dtype='int16')

This has some regular generator syntax, but the bulk of the
content comes from the ``#include`` directives.

By convention, templates have the extension ``gen.txt``, while
include files have the extension ``inc``.

Basic syntax
------------

Template files are in a line-based format.

Sections of the file are delimited with lines surrounded by ``<<`` and
``>>``.  For example::

    <<header primary>>

indicates the beginning of the primary header section.

Comments are lines beginning with ``#``.

Lines can be continued by putting a backslash character (``\``) at the
end of the line::

    DETECTOR  = { 0x1e1: 'NIR', \
                  0x1e2: 'NIR', \
                  0x1ee: 'MIR', \
                }[input('SCA_ID')] / Detector type

Other files can be included using the include directive::

    #include "other.file.txt"

.. _generator_template:

Generator template
------------------

The generator template follows this basic structure:

  - ``file`` line

  - Zero or more HDUs, each of which has

    - a ``header`` section defining how keywords are generated

    - an optional ``data`` section defining how the data is converted

``file`` line
'''''''''''''

The template must begin with a file line to give the file type a
name.  The name must be a valid Python identifier.  For example::

    <<file level1b>>

HDUs
''''

Each HDU is defined in two sections, the header and data.

Header
``````

The header begins with a header section line, giving the header a
name, which must be a valid Python identifier.  For example::

    <<header primary>>

Following that is a list of keyword definitions.  Each line is of the
form::

    KEYWORD = expression / comment

``KEYWORD`` is a FITS keyword, may be up to 8 characters, and must
contain only A through Z, ``_`` and ``-``.

The expression section is a Python expression that defines how the
keyword value is generated.  Within the namespace of the expression
are the following:

    - **Source functions:** Functions to retrieve keyword values from
      the input files.  ``input`` gets values from the input FITS
      file, and there are any number of additional functions which get
      values from the input data files.  For example, if the input
      data files include a file for program data, the function
      ``program`` is available to the expression that retrieves values
      from the program data file.  If the function is provided with no
      arguments, it retrieves the value with the same key as the
      output keyword.  If the function is provided with one argument,
      it is the name of the source keyword.  For example::

          OBS_ID = input()

      copies the OBS_ID value from the corresponding HDU in the source
      FITS file to the OBS_ID keyword in the output FITS file.  It is
      also possible to copy from a keyword value of a different name::

          CMPLTCND = input('CMPLTCON')

      To grab a value from the program data file, use the ``program``
      function instead::

          TARGET = program()

    - **Generator functions:** There are a number of helper functions
      in the ``generators`` module that help convert
      and generate values of different kinds.  For example::

          END_TIME = date_and_time_to_cds(input('DATE-END'), input('TIME-END'))

      creates a CDS value from an input date and time.

    - **Python expression syntax:** It's possible to do a lot of
      useful things just by using regular Python expression syntax.
      For example, to make the result a substring of a source
      keyword::

          PARASEQN = input('OBS_ID')[13:14] / Parallel Sequence ID

      or to calculate the difference of two values::

          DURATION = input('END_TIME') - input('START_TIME')

The optional comment section following a ``/`` character will be
attached to the keyword in the output FITS file.  There is an
important distinction between these comments which end up in the
output FITS file, and comments beginning with ``#`` which are included
in the template for informational purposes only and are ignored by the
template parser.

It is also possible to include comments on their own lines to create
section headings in the output FITS file.  For example::

    / MIRI-specific keywords
    FILTER    = '' / Filter element used
    FLTSUITE  = '' / Flat field element used
    WAVLNGTH  = '' / Wavelength requested in the exposure specification
    GRATING   = '' / Grating/dichroic wheel position
    LAMPON    = '' / Internal calibration lamp
    CCCSTATE  = '' / Contamination control cover state

    / Exposure parameters
    READPATT  = '' / Readout pattern
    NFRAME    = 1 / Number of frames per read group
    NSKIP     = 0 / Number of frames dropped
    FRAME0    = 0 / zero-frame read
    INTTIME   = 0 / Integration time
    EXPTIME   = 0 / Exposure time
    DURATION  = 0 / Total duration of exposure
    OBJ_TYPE  = 'FAINT' / Object type

``#include`` files will typically be just lines defining keyword definitions
as above, for example, the file ``target_information.gen.inc`` looks like this::

    / Target information
    
    TARGPROP = input('TARGNAME') / proposer's name for the target
    TARGNAME = 'NGC 104' / standard astronomical catalog name for target
    TARGTYPE = 'FIXED' / fixed target, moving target, or generic target
    TARG_RA  = 0.0 / target RA computed at time of exposure
    TARGURA  = 0.0 / target RA uncertainty
    TARG_DEC = 0.0 / target DEC computed at time of exposure
    TARRUDEC =  0.0  / target Dec uncertainty
    PROP_RA  =  0.0  / proposer specified RA for the target
    PROP_DEC =  0.0  / proposer specified Dec for the target
    PROPEPOC = 2000.0  / proposer specified epoch for RA and Dec

and is used in many of the top-level level1b templates.

Data
````
The data section consists of a single expression that returns a Numpy
array containing the output data.

The following are available in the namespace:

  - ``np``: ``import numpy as np``

  - ``input``: A fits HDUList object containing the content of the
    input FITS file.

  - ``output``: A fits HDUList object containing the content of the
    output FITS file.  Note that the output FITS file may only be
    partially contructed.  Importantly, higher-number HDUs will not
    yet exist.

A complete example
''''''''''''''''''

::

  # This file defines the structure of a MIRI level 1b file
  <<file miri_level1b>>
  <<header primary>>
  SIMPLE    = T
  BITPIX    = 32
  NAXIS     = 0
  EXTEND    = T
  ORIGIN    = 'STScI'
  TELESCOP  = 'JWST'
  FILENAME  = '' / The filename
  DATE      = now() / Date this file was generated

  #include "level1a.gen.inc"

  #include "level1b.gen.inc"

  / MIRI-specific keywords
  FILTER    = '' / Filter element used
  FLTSUITE  = '' / Flat field element used
  WAVLNGTH  = '' / Wavelength requested in the exposure specification
  GRATING   = '' / Grating/dichroic wheel position
  LAMPON    = '' / Internal calibration lamp
  CCCSTATE  = '' / Contamination control cover state

  / Exposure parameters
  READPATT  = '' / Readout pattern
  NFRAME    = 1 / Number of frames per read group
  NSKIP     = 0 / Number of frames dropped
  FRAME0    = 0 / zero-frame read
  INTTIME   = 0 / Integration time
  EXPTIME   = 0 / Exposure time
  DURATION  = 0 / Total duration of exposure
  OBJ_TYPE  = 'FAINT' / Object type

  / Subarray parameters
  SUBARRAY  = '' / Name of subarray used
  SUBXSTRT  = 0 / x-axis pixel number of subarray origin
  SUBXSIZE  = 0 / length of subarray along x-axis
  SUBTSTRT  = 0 / y-axis pixel number of subarray origin
  SUBYSIZE  = 0 / length of subarray along y-axis
  LIGHTCOL  = 0 / Number of light-sensitive columns

  <<data>>

  <<header science>>
  XTENSION  = 'IMAGE' /        FITS extension type
  BITPIX    =         /        bits per data value
  NAXIS     =         /        number of data array dimensions
  NAXIS1    =         /        length of first data axis (#columns)
  NAXIS2    =         /        length of second data axis (#rows)
  NAXIS3    =         /        length of third data axis (#groups/integration)
  NAXIS4    =         /        length of fourth data axis (#integrations)
  PCOUNT    = 0       /        number of parameter bytes following data table
  GCOUNT    = 1       /        number of groups
  EXTNAME   = 'SCI'   /        extension name
  BSCALE    = 1.0     /        scale factor for array value to physical value
  BZERO     = 32768   /        physical value for an array value of zero
  BUNIT     = 'DN'    /        physical units of the data array values

  <<data>>
  input[0].data.reshape((input[0].header['NINT'], \
                         input[0].header['NGROUP'], \
                         input[0].header['NAXIS2'], \
                         input[0].header['NAXIS1'])). \
                        astype('uint16')
