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


class Asn_Lv2_Image(DMS_Level2b_Base):
    """Level 2 Image association"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': 'NRC_IMAGE|MIR_IMAGE|NIS_IMAGE|FGS_IMAGE',
                'inputs': ['EXP_TYPE'],
                'force_unique': True,
            }
        })

        super(Asn_Lv2_Image, self).__init__(*args, **kwargs)
