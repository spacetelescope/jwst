"""Association Definitions: DMS Level2b product associations
"""
import logging

from jwst.associations.lib.rules_level2_base import *

__all__ = [
    'Asn_Lv2_Image',
    'Asn_Lv2ImageNonScience',
    'Asn_Lv2ImageSpecial',
    'Asn_Lv2Spec',
    'Asn_Lv2SpecNonScience',
    'Asn_Lv2SpecSpecial'
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
        AsnMixin_Lv2Mode
):
    """Level2b Image"""


class Asn_Lv2ImageNonScience(
        AsnMixin_Lv2Singleton,
        AsnMixin_Lv2ImageNonScience,
        AsnMixin_Lv2Mode
):
    """Level2b Image"""


class Asn_Lv2ImageSpecial(
        AsnMixin_Lv2Special,
        AsnMixin_Lv2Singleton,
        AsnMixin_Lv2Image,
        AsnMixin_Lv2Mode,
):
    """Level2b Image that are marked special

    Image exposures that are marked as backgrounds, imprints, etc.,
    still get 2b processing just as normal science. However, no other
    exposures should get included into the association.

    """


class Asn_Lv2Spec(
        AsnMixin_Lv2Singleton,
        AsnMixin_Lv2Spec,
        AsnMixin_Lv2Mode,
):
    """Level2b Spectra"""


class Asn_Lv2SpecNonScience(
        AsnMixin_Lv2Singleton,
        AsnMixin_Lv2SpecNonScience,
        AsnMixin_Lv2Mode,
):
    """Level2b Spectra"""


class Asn_Lv2SpecSpecial(
        AsnMixin_Lv2Special,
        AsnMixin_Lv2Singleton,
        AsnMixin_Lv2Spec,
        AsnMixin_Lv2Mode,
):
    """Level2b Spectra that are marked special

    Spectral exposures that are marked as backgrounds, imprints, etc.,
    still get 2b processing just as normal science. However, no other
    exposures should get included into the association.

    """
