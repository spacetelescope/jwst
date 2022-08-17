"""SIAF Database Access

Provide a common interface to different versions of the SIAF.

Under operations, the SIAF is found in a sqlite database.
Otherwise, use the standard interface defined by the `pysiaf` package
"""
from collections import namedtuple
from datetime import date
import logging
import os
from pathlib import Path

from .basic_utils import LoggingContext

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Map instrument three character mnemonic to full name
INSTRUMENT_MAP = {
    'fgs': 'fgs',
    'mir': 'miri',
    'nis': 'niriss',
    'nrc': 'nircam',
    'nrs': 'nirspec'
}

# SIAF container
# The names should correspond to the names in the ``wcsinfo`` schema.
# It is populated by the SIAF values in the PRD database based
# on APERNAME and UseAfterDate and used to populate the keywords
# in Level1bModel data models.
SIAF = namedtuple("SIAF", ["v2_ref", "v3_ref", "v3yangle", "vparity",
                           "crpix1", "crpix2", "cdelt1", "cdelt2",
                           "vertices_idl"])
# Set default values for the SIAF.
# Values which are needed by the pipeline are set to None which
# triggers a ValueError if missing in the SIAF database.
# Quantities not used by the pipeline get a default value -
# FITS keywords and aperture vertices.
SIAF.__new__.__defaults__ = (None, None, None, None, 0, 0, 3600, 3600,
                             (0, 1, 1, 0, 0, 0, 1, 1))

SIAF_REQUIRED = ['V2Ref', 'V3Ref', 'V3IdlYAngle', 'VIdlParity']
SIAF_OPTIONAL = ['XSciRef', 'YSciRef', 'XSciScale', 'YSciScale']
SIAF_VERTICIES = ['XIdlVert1', 'XIdlVert2', 'XIdlVert3', 'XIdlVert4',
                  'YIdlVert1', 'YIdlVert2', 'YIdlVert3', 'YIdlVert4']
SIAF_MAP = {'V2Ref': 'v2_ref', 'V3Ref': 'v3_ref', 'V3IdlYAngle': 'v3yangle', 'VIdlParity': 'vparity',
            'XSciRef': 'crpix1', 'YSciRef': 'crpix2', 'XSciScale': 'cdelt1', 'YSciScale': 'cdelt2'}

class SiafDb:
    """Use pysiaf as the source of siaf information

    Parameters
    ----------
    source : None, str, or a file-like object
        If None, then the latest PRD version in `pysiaf` is used.
        Otherwise, it should be a string or Path-like object pointing to a folder containing the
        SIAF XML files.

    prd : None or str
        The PRD version to use from the pysiaf application. If `source` has also been
        specified, `source` will be used instead.

    Notes
    -----
    The interpretation of `source` is as follows:

    """
    def __init__(self, source=None, prd=None):
        logger_pysiaf = logging.getLogger('pysiaf')
        log_level = logger_pysiaf.getEffectiveLevel()
        if not source:
            log_level = logging.ERROR
        try:
            with LoggingContext(logger_pysiaf, level=log_level):
                import pysiaf
        except ImportError:
            raise ValueError('Package "pysiaf" is not installed. Cannot use the pysiaf api')
        self.pysiaf = pysiaf

        self.xml_path = self.get_xml_path(source, prd)

    def get_aperture(self, aperture, useafter=None):
        """Get the pysiaf.Aperture for an aperture

        Parameters
        ----------
        aperture : str
            The name of the aperture to retrieve.
        useafter : str
            The date of observation (``model.meta.date``)

        Returns
        -------
        aperture : pysiaf.Aperture
            The aperture specification.
        """
        if not useafter:
            useafter = date.today().strftime('%Y-%m-%d')

        instrument = INSTRUMENT_MAP[aperture[:3].lower()]
        siaf = self.pysiaf.Siaf(instrument, basepath=self.xml_path)
        aperture = siaf[aperture.upper()]
        return aperture

    def get_wcs(self, aperture, useafter=None):
        """
        Query the SIAF database file and get WCS values.

        Given an ``APERTURE_NAME`` and a ``USEAFTER`` date query the SIAF database
        and extract the following keywords:
        ``V2Ref``, ``V3Ref``, ``V3IdlYAngle``, ``VIdlParity``,
        ``XSciRef``, ``YSciRef``, ``XSciScale``, ``YSciScale``,
        ``XIdlVert1``, ``XIdlVert2``, ``XIdlVert3``, ``XIdlVert4``,
        ``YIdlVert1``, ``YIdlVert2``, ``YIdlVert3``, ``YIdlVert4``

        Parameters
        ----------
        aperture : str
            The name of the aperture to retrieve.
        useafter : str
            The date of observation (``model.meta.date``)

        Returns
        -------
        siaf : namedtuple
            The SIAF namedtuple with values from the PRD database.
        """
        aperture = self.get_aperture(aperture, useafter=useafter)

        # Build the SIAF entry. Missing required values is an error.
        # Otherwise, use defaults.
        default_siaf = SIAF()
        values = {SIAF_MAP[key]: getattr(aperture, key) for key in SIAF_REQUIRED}
        if not all(values):
            raise RuntimeError(f'Required SIAF entries for {aperture} are not all defined: {values}')
        for key in SIAF_OPTIONAL:
            value = getattr(aperture, key)
            value = value if value else getattr(default_siaf, SIAF_MAP[key])
            values[SIAF_MAP[key]] = value
        vertices = list()
        for key in SIAF_VERTICIES:
            value = getattr(aperture, key)
            value = value if value else getattr(default_siaf, SIAF_MAP[key])
            vertices.append(value)
        vertices = tuple(vertices)

        # Fill out the Siaf
        siaf = SIAF(**values, vertices_idl=vertices)

        return siaf

    def get_xml_path(self, source, prd):
        """Determine the XML source to use

        Parameters
        ----------
        source : None, str, or a file-like object
            If None, then the latest PRD version in `pysiaf` is used.
            Otherwise, it should be a string or Path-like object pointing to a folder containing the
            SIAF XML files.

        prd : None or str
            The PRD version to use from the pysiaf application. If `source` has also been
            specified, `source` will be used instead.

        Returns
        -------
        xml_path : Path
            Either the Path to the XML files.

        Raises
        ------
        ValueError
            If `source` does not resolve to a folder or `prd` is not a valid PRD version.
        """

        # If `source` is defined and valid, use that.
        xml_path = None
        if source is not None:
            xml_path = Path(source)
            if not xml_path.is_dir():
                raise ValueError('Source %s: Needs to be a folder for use with pysiaf', xml_path)

        # If a PRD version is defined, attempt to use that.
        if not xml_path and prd:
            prd = prd.upper()
            prd_root = Path(self.pysiaf.JWST_PRD_DATA_ROOT).parent.parent.parent
            if not (prd_root / prd).is_dir():
                raise ValueError('PRD specification %s does not exist', prd)
            xml_path = prd_root / prd / 'SIAFXML' / 'SIAFXML'

        # If nothing has been specified, see if XML_DATA says what to do.
        if not xml_path:
            xml_path = os.environ.get('XML_DATA', None)
            if xml_path:
                xml_path = Path(xml_path) / 'SIAFXML'

        # If nothing else, use the `pysiaf` default.
        if not xml_path:
            xml_path = Path(self.pysiaf.JWST_PRD_DATA_ROOT)
            logger.info('pysiaf: Using latest installed PRD %s', self.pysiaf.JWST_PRD_VERSION)
        else:
            logger.info('pysiaf: Using SIAF XML folder %s', xml_path)

        return xml_path
