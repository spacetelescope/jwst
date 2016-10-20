"""Association Definitions: DMS Level3 product associations
"""
import logging

from jwst.associations.lib.rules_level3_base import *

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
                'inputs': ['WFSVISIT'],
                'is_invalid': True,
            },
        })

        # Now check and continue initialization.
        super(Asn_Image, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image'
        super(Asn_Image, self)._init_hook(member)


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
                'value': '(?!NULL).+',
                'inputs': ['WFSVISIT'],
            },
            'asn_candidate_wfs': {
                'value': '.+MOSAIC.+',
                'inputs': ['ASN_CANDIDATE'],
                'force_unique': True,
                'is_acid': True,
            },
            'activity_id': {
                'value': None,
                'inputs': ['ACT_ID']
            }
        })

        # Now check and continue initialization.
        super(Asn_WFSCMB, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'wfs'
        super(Asn_WFSCMB, self)._init_hook(member)


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
                'value': 'MIR_LRS-FIXEDSLIT|MIR_TACQ',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': 'P750L',
                'inputs': ['FILTER']
            },
            'subarray': {
                'value': 'FULL',
                'inputs': ['SUBARRAY']
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
                'value': 'MIR_LRS-SLITLESS|MIR_TACQ',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': 'P750L',
                'inputs': ['FILTER']
            },
            'subarray': {
                'value': 'SUBPRISM',
                'inputs': ['SUBARRAY']
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
                'value': 'NIS',
                'inputs': ['DETECTOR']
            },
            'exp_type': {
                'value': 'NIS_SOSS|NIS_TACQ|NIS_TACNFRM',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': 'GR700XD',
                'inputs': ['PUPIL']
            },
            'subarray': {
                'value': 'FULL|SUBSTRIP256|SUBSTRIP80',
                'inputs': ['SUBARRAY'],
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
                    'NRS_FIXEDSLIT'
                    '|NRS_AUTOWAVE'
                    '|NRS_CONFIRM'
                    '|NRS_TACQ'
                    '|NRS_TACONFIRM'
                    '|NRS_TASLIT'
                ),
                'inputs': ['EXP_TYPE']
            },
            'fixed_slit': {
                'value': None,
                'inputs': ['FXD_SLIT']
            },
            'subarray': {
                'value': None,
                'inputs': ['SUBARRAY']
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
                    'NRS_MSASPEC'
                    '|NRS_AUTOFLAT'
                    '|NRS_AUTOWAVE'
                    '|NRS_CONFIRM'
                    '|NRS_TASLIT'
                    '|NRS_TACQ'
                    '|NRS_TACONFIRM'
                ),
                'inputs': ['EXP_TYPE']
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
                    'MIR_MRS'
                    '|MIR_FLATMRS'
                    '|MIR_TACQ'
                ),
                'inputs': ['EXP_TYPE'],
                'force_unique': False,
            },
        })

        # Check and continue initialization.
        super(Asn_MIRI_IFU, self).__init__(*args, **kwargs)

    def product_name(self):
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

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        super(Asn_MIRI_IFU, self)._init_hook(member)
        self.data['asn_type'] = 'mirifu'


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
                    'NRS_IFU'
                    '|NRS_AUTOWAVE'
                    '|NRS_CONFIRM'
                    '|NRS_TASLIT'
                    '|NRS_TACQ'
                    '|NRS_TACONFIRM'
                ),
                'inputs': ['EXP_TYPE']
            },
        })

        # Check and continue initialization.
        super(Asn_NRS_IFU, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        super(Asn_NRS_IFU, self)._init_hook(member)
        self.data['asn_type'] = 'nrsifu'
