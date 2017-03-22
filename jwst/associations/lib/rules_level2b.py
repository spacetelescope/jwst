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
