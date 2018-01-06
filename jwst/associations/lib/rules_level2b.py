"""Association Definitions: DMS Level2b product associations
"""
import logging

from jwst.associations.lib.constraint import Constraint
from jwst.associations.lib.rules_level2_base import *

__all__ = [
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------
# Start of the User-level rules
# --------------------------------
class Asn_Lv2Image(
        AsnMixin_Lv2Singleton,
        AsnMixin_Lv2Image,
        DMSLevel2bBase
):
    """Level2b Image"""

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            CONSTRAINT_BASE,
            CONSTRAINT_MODE,
            CONSTRAINT_IMAGE_SCIENCE,
        ])

        # Now check and continue initialization.
        super(Asn_Lv2Image, self).__init__(*args, **kwargs)
