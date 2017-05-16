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


class Asn_Lv2Image(
        AsnMixin_Lv2Singleton,
        AsnMixin_Lv2Image,
        AsnMixin_Lv2Mode
):
    """Level2b Image"""


class Asn_Lv2Spec(
        AsnMixin_Lv2Singleton,
        AsnMixin_Lv2Spec,
        AsnMixin_Lv2Mode
):
    """Level2b Spectra"""


class Asn_Lv2SpecBkg(
        AsnMixin_Lv2Spec,
        AsnMixin_Lv2Bkg,
        AsnMixin_Lv2Mode
):
    """Level2b Spectra with backgrounds"""

    def __init__(self, *args, **kwargs):
        self.validity.update({
            'has_background': {
                'validated': False,
                'check': lambda entry: entry['exptype'] == 'BACKGROUND'
            }
        })

        super(Asn_Lv2SpecBkg, self).__init__(*args, **kwargs)
