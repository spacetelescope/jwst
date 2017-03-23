"""Association Definitions: DMS Level2b product associations
"""
import logging

from jwst.associations.lib.rules_level2_base import *

__all__ = ['Asn_Lv2_Image']

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# --------------------------------
# Start of the User-level rules
# --------------------------------


class Asn_Lv2Image(DMSLevel2bBase):
    """Level 2 Image association"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': 'NRC_IMAGE|MIR_IMAGE|NIS_IMAGE|FGS_IMAGE',
                'inputs': ['EXP_TYPE'],
                'force_unique': False,
            }
        })

        super(Asn_Lv2Image, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2Image, self)._init_hook(member)
        self.data['asn_type'] = 'image2'


class Asn_Lv2Spec(DMSLevel2bBase):
    """Level 2 Spectral association"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': (
                    'NRC_GRISM'
                    '|NRC_TSGRISM'
                    '|MIR_LRS-FIXEDSLIT'
                    '|MIR_LRS-SLITLESS'
                    '|NRS_FIXEDSLIT'
                    '|NRS_IFU'
                    '|NRS_MSASPEC'
                    '|NRS_BRIGHTOBJ'
                    '|NIS_WFSS'
                    '|NIS_SOSS'
                ),
                'inputs': ['EXP_TYPE'],
                'force_unique': False,
            }
        })

        super(Asn_Lv2Spec, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2Spec, self)._init_hook(member)
        self.data['asn_type'] = 'spec2'
