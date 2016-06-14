.. _template:

Template file format
====================

File types are described using a simple file format that vaguely
resembles FITS headers.

Defining a new filetype requires writing two templates: one for
:ref:`generation <generator_template>` and one for
:ref:`verification <verifier_template>`.  The generation template is 
the "driver", defining the order and existence of keywords.  The 
verification template keywords are considered unordered, and may contain
definitions of keywords not ultimately used in the file format.  This
allows each specific format to be defined in a generator template, and
avoids duplication of effort by including verification definitions
from the same set of files.

A third form of template also exists that is used as :ref:`input data
<data_files>`.

By convention, generator templates have the extension ``gen.txt``,
verifier templates have the extension ``val.txt`` and data files have
the extension ``dat``.

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

      If an APT file is provided as input, the ``apt()`` function is
      available.  It takes as input an `XPath expression
      <http://effbot.org/zone/element-xpath.htm>`_ to look up an
      element in the APT XML file.  For example::

          PI = apt("./ProposalInformation/PrincipalInvestigator/InvestigatorAddress/LastName") + ", " + \
               apt("./ProposalInformation/PrincipalInvestigator/InvestigatorAddress/FirstName")

    - **Generator functions:** There are a number of helper functions
      in the :ref:`generators <generators>` module that help convert
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

Data
````
The data section consists of a single expression that returns a Numpy
array containing the output data.

The following are available in the namespace:

  - ``np``: ``import numpy as np``

  - ``input``: A pyfits HDUList object containing the content of the
    input FITS file.

  - ``output``: A pyfits HDUList object containing the content of the
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

.. _verifier_template:

Verifier template
-----------------

The verifier template follows this basic structure:

  - ``file`` line

  - Any number of optional ``inherit`` lines

  - Zero or more HDUs, each of which has

    - a ``header`` section defining how keywords are generated

    - an optional ``data`` section defining how the data is converted

``file`` line
'''''''''''''

The template must begin with a file line to give the file type a
name.  The name must be a valid Python identifier.  For example::

    <<file level1b>>

``inherit`` lines
'''''''''''''''''

The ``inherit`` lines include another file containing keyword
definitions that this file will inherit from.  For example, the basic
FITS keywords may be defined in a centrally located file that all of
the file types inherit from.  The specific file types may override any
of the inherited definitions.

.. note::
   ``inherit`` works differently from ``#include``.  ``inherit`` loads
   in an entire set of keyword definitions that are used as fallbacks
   for keywords not defined in a given file.  ``#include`` merely
   includes another file verbatim.

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

The expression section is a Python expression that is evaluated to see
whether the value is valid.  The expression should be a boolean
expression returning `True` if the value is valid, otherwise `False`.
Within the namespace of the expression are the following:

    - **Value variable x:** The variable ``x`` is available in the
      namespace and stores the value of the keyword.  You can easily
      test if the value is a particular constant using::

          NAXIS = x == 4

      Or that a value is in a particular range::

          NAXIS = 1 <= x <= 4

      Or ensure a value is a member of a particular set::

          BITPIX = x in (8, 16, 32, 64, -32, -64)

    - **Accessing other keywords:** Keywords can be compared to other
      keywords by using the ``output`` function.  For example, to
      ensure that a keyword has the same value as another::

          NINT = x == output('NAXIS4')

    - **Verifier functions:** There are a number of helper functions
      in the :ref:`verifiers <verifiers>` module that test certain
      properties of the value.  For example::

          DATE = is_date(x)

Optional comments may be added to the definition line, but they are
ignored.  Only the comments in the generator template are written to
the output FITS file.

Data
````

TODO: The functionality here has not been fleshed out.

.. _data_files:

Data files
----------

A data file follows this basic structure:

  - ``file`` line

  - Zero or more ``header`` sections containing keyword values.  Note
    data files do not contain ``data`` sections.

``file`` line
'''''''''''''

The template must begin with a file line to give the file type a
name.  The name must be a valid Python identifier.

This name is used to indicate what type of data file this is.  For
example, if the data file contains the line::

    <<file program>>

then a function ``program`` is available in the generator template to
pull values from this data file.

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

The expression section is a Python literal expression.  It is
evaluated at file load time.  Here are examples for all of the basic
FITS datatypes::

    TARGNAME  = "R2-D2"
    RA_TARG   = 32.19
    COUNT     = 42
    LAMPON    = T
