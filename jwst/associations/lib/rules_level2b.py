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


class Asn_Lv2Image(AsnMixin_Lv2Image):
    """Level 2 Image association"""


class Asn_Lv2Spec(AsnMixin_Lv2Spec):
    """Level 2 Spectral association"""


class Asn_Lv2SpecBkg(AsnMixin_Lv2Spec, AsnMixin_Lv2Mode):
    """Level2b Image with backgrounds"""

    def __init__(self, *args, **kwargs):
        self.add_constraints({
            'background': {
                'inputs': ['ASN_CANDIDATE'],
                'value': '.+BACKGROUND.+',
                'force_unique': True,
                'is_acid': False,
            }
        })

        # Now, lets see if member belongs to us.
        super(Asn_Lv2SpecBkg, self).__init__(*args, **kwargs)

    def _add(self, member, check_extra_flags=None):
        if not check_extra_flags:
            check_extra_flags = []
        check_extra_flags.append('BACKGROUND')
        super(Asn_Lv2SpecBkg, self)._add(member, check_extra_flags)
