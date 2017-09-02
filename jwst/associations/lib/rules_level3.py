"""Association Definitions: DMS Level3 product associations
"""
import logging

from jwst.associations.lib.rules_level3_base import *

__all__ = [
    'Asn_Image',
    'Asn_MIRI_IFU',
    'Asn_MIRI_LRS_FIXEDSLIT',
    'Asn_MIRI_LRS_SLITLESS',
    'Asn_NRS_FIXEDSLIT',
    'Asn_NRS_IFU',
    'Asn_NRS_MSA',
    'Asn_NIR_SO_SLITLESS',
    'Asn_WFSCMB',

]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# --------------------------------
# Start of the User-level rules
# --------------------------------


# ----------------------------------
# Image associations
class Asn_Image(
        AsnMixin_Image,
        AsnMixin_OpticalPath,
        AsnMixin_Target,
        AsnMixin_Base
):
    """Non-Association Candidate Dither Associations"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'wfsvisit': {
                'inputs': ['wfsvisit'],
                'force_undefined': True,
            },
        })

        # Now check and continue initialization.
        super(Asn_Image, self).__init__(*args, **kwargs)


class Asn_WFSCMB(
        AsnMixin_Image,
        AsnMixin_Target,
        AsnMixin_Base
):
    """Wavefront Sensing association

    Notes
    -----
    Defined by `TRAC issue #269 <https://aeon.stsci.edu/ssb/trac/jwst/ticket/269>`_
    """
    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'wfsvisit': {
                'value': '(?!null).+',
                'inputs': ['wfsvisit'],
            },
            'asn_candidate_wfs': {
                'value': '.+mosaic.+',
                'inputs': ['asn_candidate'],
                'force_unique': True,
                'is_acid': True,
                'evaluate': True,
            },
            'activity_id': {
                'value': None,
                'inputs': ['act_id']
            }
        })

        # Now check and continue initialization.
        super(Asn_WFSCMB, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'wfs'
        super(Asn_WFSCMB, self)._init_hook(item)


# Spectrographic Associations
class Asn_MIRI_LRS_FIXEDSLIT(
        AsnMixin_Spectrum,
        AsnMixin_MIRI,
        AsnMixin_Target,
        AsnMixin_Base
):
    """MIRI LRS Fixed slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'mir_lrs-fixedslit|mir_tacq',
                'inputs': ['exp_type']
            },
            'opt_elem': {
                'value': 'p750l',
                'inputs': ['filter']
            },
            'subarray': {
                'value': 'full',
                'inputs': ['subarray']
            }
        })

        # Check and continue initialization.
        super(Asn_MIRI_LRS_FIXEDSLIT, self).__init__(*args, **kwargs)


class Asn_MIRI_LRS_SLITLESS(
        AsnMixin_Spectrum,
        AsnMixin_MIRI,
        AsnMixin_Target,
        AsnMixin_Base
):
    """MIRI LRS Slitless"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'mir_lrs-slitless|mir_tacq',
                'inputs': ['exp_type']
            },
            'opt_elem': {
                'value': 'p750l',
                'inputs': ['filter']
            },
            'subarray': {
                'value': 'subprism',
                'inputs': ['subarray']
            }
        })

        # Check and continue initialization.
        super(Asn_MIRI_LRS_SLITLESS, self).__init__(*args, **kwargs)


class Asn_NIR_SO_SLITLESS(
        AsnMixin_Spectrum,
        AsnMixin_NIRISS,
        AsnMixin_Target,
        AsnMixin_Base
):
    """NIRISS Single-Object Slitless"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'detector': {
                'value': 'nis',
                'inputs': ['detector']
            },
            'exp_type': {
                'value': 'nis_soss|nis_tacq|nis_tacnfrm',
                'inputs': ['exp_type']
            },
            'opt_elem': {
                'value': 'gr700xd',
                'inputs': ['pupil']
            },
            'subarray': {
                'value': 'full|substrip256|substrip80',
                'inputs': ['subarray'],
                'force_unique': True
            }
        })

        # Check and continue initialization.
        super(Asn_NIR_SO_SLITLESS, self).__init__(*args, **kwargs)


class Asn_NRS_FIXEDSLIT(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_OpticalPath,
        AsnMixin_Target,
        AsnMixin_Base
):
    """NIRSPEC Fixed Slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': (
                    'nrs_fixedslit'
                    '|nrs_autowave'
                    '|nrs_confirm'
                    '|nrs_tacq'
                    '|nrs_taconfirm'
                    '|nrs_taslit'
                ),
                'inputs': ['exp_type']
            },
            'fixed_slit': {
                'value': None,
                'inputs': ['fxd_slit']
            },
            'subarray': {
                'value': None,
                'inputs': ['subarray']
            },
        })

        # Check and continue initialization.
        super(Asn_NRS_FIXEDSLIT, self).__init__(*args, **kwargs)


class Asn_NRS_MSA(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_OpticalPath,
        AsnMixin_Target,
        AsnMixin_Base
):
    """NIRSPEC MSA"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': (
                    'nrs_msaspec'
                    '|nrs_autoflat'
                    '|nrs_autowave'
                    '|nrs_confirm'
                    '|nrs_taslit'
                    '|nrs_tacq'
                    '|nrs_taconfirm'
                ),
                'inputs': ['exp_type']
            },
        })

        # Check and continue initialization.
        super(Asn_NRS_MSA, self).__init__(*args, **kwargs)


class Asn_MIRI_IFU(
        AsnMixin_Spectrum,
        AsnMixin_MIRI,
        AsnMixin_Target,
        AsnMixin_Base
):
    """MIRI MRS (IFU)"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': (
                    'mir_mrs'
                    '|mir_flatmrs'
                    '|mir_tacq'
                ),
                'inputs': ['exp_type'],
                'force_unique': False,
            },
        })

        # Check and continue initialization.
        super(Asn_MIRI_IFU, self).__init__(*args, **kwargs)

    def dms_product_name(self):
        """Define product name."""
        target = self._get_target()

        instrument = self._get_instrument()

        product_name = 'jw{}-{}_{}_{}'.format(
            self.data['program'],
            self.acid.id,
            target,
            instrument
        )

        return product_name.lower()


class Asn_NRS_IFU(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_Target,
        AsnMixin_Base
):
    """NIRSPEC IFU"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': (
                    'nrs_ifu'
                    '|nrs_autowave'
                    '|nrs_confirm'
                    '|nrs_taslit'
                    '|nrs_tacq'
                    '|nrs_taconfirm'
                ),
                'inputs': ['exp_type']
            },
        })

        # Check and continue initialization.
        super(Asn_NRS_IFU, self).__init__(*args, **kwargs)


class Asn_Coron(
        AsnMixin_OpticalPath,
        AsnMixin_Base
):
    """Coronography

    Notes
    -----

    Coronography is nearly completely defined by the association candidates
    produced by APT.

    Tracking Issues:

    - `github #311 <https://github.com/STScI-JWST/jwst/issues/311>`
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': (
                    'nrc_coron'
                    '|mir_lyot'
                    '|mir_4qpm'
                ),
                'inputs': ['exp_type'],
                'force_unique': True,
            },
            'target': {
                'value': None,
                'inputs': ['targetid'],
                'onlyif': lambda item: self.get_exposure_type(item) == 'science'
            }
        })

        # PSF is required
        self.validity.update({
            'has_psf': {
                'validated': False,
                'check': lambda entry: entry['exptype'] == 'psf'
            }
        })

        # Check and continue initialization.
        super(Asn_Coron, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'coron3'
        super(Asn_Coron, self)._init_hook(item)
