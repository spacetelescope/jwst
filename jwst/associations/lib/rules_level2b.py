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


class Asn_Lv2Image(AsnMixin_Lv2Mode, AsnMixin_Lv2Image):
    """Level 2 Image association"""


class Asn_Lv2Spec(AsnMixin_Lv2Mode, AsnMixin_Lv2Spec):
    """Level 2 Spectral association"""
