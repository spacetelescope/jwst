"""Association Definitions: DMS Level3 product associations."""  # noqa: INP001

from jwst.associations import Association
from jwst.associations.registry import RegistryMarker

# The schema that these associations must adhere to.
_ASN_SCHEMA_LEVEL3 = "asn_schema_jw_level3.json"
_DMS_POOLNAME_REGEX = r"jw(\d{5})_(\d{8}[Tt]\d{6})_pool"


class DMSLevel3BaseSet1(Association):
    """Basic class for DMS Level3 associations."""


@RegistryMarker.rule
class AsnDitherSet1(DMSLevel3BaseSet1):
    """Non-Association Candidate Dither Associations."""


@RegistryMarker.rule
class AsnWFSSet1(DMSLevel3BaseSet1):
    """Wavefront Sensing association."""
