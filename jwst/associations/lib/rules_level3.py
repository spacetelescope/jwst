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
# Associations defined by Candidates
class Asn_Mosaic(
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """Association Candidate type of Mosaic"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'asn_candidate_id': {
                'value': None,
                'inputs': ['ASN_CANDIDATE_ID'],
            },
            'asn_candidate_type': {
                'value': 'MOSAIC',
                'inputs': ['ASN_CANDIDATE_TYPE'],
            },
            'wfsvisit': {
                'value': 'NULL',
                'inputs': ['WFSVISIT']
            }
        })

        # Now check and continue initialization.
        super(Asn_Mosaic, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image'
        super(Asn_Mosaic, self)._init_hook(member)


# Image associations
class Asn_Dither(
        AsnMixin_Image,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """Non-Association Candidate Dither Associations"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'pointing_type': {
                'value': 'SCIENCE',
                'inputs': ['PNTGTYPE']
            },
            'wfsvisit': {
                'value': 'NULL',
                'inputs': ['WFSVISIT'],
            },
            'asn_candidate_type': {
                'value': 'OBSERVATION',
                'inputs': ['ASN_CANDIDATE_TYPE'],
                'required': False,
            },
        })

        # Now check and continue initialization.
        super(Asn_Dither, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image'
        super(Asn_Dither, self)._init_hook(member)


class Asn_WFSCMB(
        AsnMixin_Image,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """Wavefront Sensing association"""
    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'wfsvisit': {
                'value': '(?!NULL).+',
                'inputs': ['WFSVISIT'],
            },
            'asn_candidate_type': {
                'value': '(?!OBSERVATION).+',
                'inputs': ['ASN_CANDIDATE_TYPE']
            },
            'asn_candidate_id': {
                'value': None,
                'inputs': ['ASN_CANDIDATE_ID']
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
        AsnMixin_Unique_Config
):
    """MIRI LRS Fixed slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'patttype': {
                'value': None,
                'inputs': ['PATTTYPE'],
                'force_unique': True
            },
            'exp_type': {
                'value': 'MIR_LRS-FIXEDSLIT',
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
        AsnMixin_Unique_Config
):
    """MIRI LRS Slitless"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'MIR_LRS-SLITLESS',
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
        AsnMixin_Unique_Config
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
                'value': 'NIS_SOSS',
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


class Asn_NIR_FIXEDSLIT(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """NIRSPEC Fixed Slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'NRS_FIXEDSLIT',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER']
            },
            'grating': {
                'value': None,
                'inputs': ['GRATING']
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
        super(Asn_NIR_FIXEDSLIT, self).__init__(*args, **kwargs)


class Asn_NIR_MSA(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """NIRSPEC MSA"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'pointing_type': {
                'value': 'SCIENCE',
                'inputs': ['PNTGTYPE']
            },
            'exp_type': {
                'value': 'NRS_MSASPEC',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER']
            },
            'grating': {
                'value': None,
                'inputs': ['GRATING']
            },
        })

        # Check and continue initialization.
        super(Asn_NIR_MSA, self).__init__(*args, **kwargs)


class Asn_MIRI_MRS(
        AsnMixin_Spectrum,
        AsnMixin_MIRI,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """MIRI MRS (IFU)"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'MIR_MRS',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['BAND']
            },
        })

        # Check and continue initialization.
        super(Asn_MIRI_MRS, self).__init__(*args, **kwargs)


class Asn_NIR_IFU(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """NIRSPEC IFU"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'NRS_IFU',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER']
            },
            'grating': {
                'value': None,
                'inputs': ['GRATING']
            }
        })

        # Check and continue initialization.
        super(Asn_NIR_IFU, self).__init__(*args, **kwargs)
