"""Association Definitions: DMS Level2b product associations
"""
import logging

from jwst.associations.lib.dms_base import format_list
from jwst.associations.lib.rules_level2_base import *
from jwst.associations.lib.rules_level3_base import (DMS_Level3_Base, _EMPTY)

__all__ = [
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# --------------------------------
# Start of the User-level rules
# --------------------------------
