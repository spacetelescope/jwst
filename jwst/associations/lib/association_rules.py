"""Association Definitions: DMS-specific

Notes
-----
These associations are specifically defined for use in DMS.
"""
from jwst.associations import Association
from jwst.associations.registry import RegistryMarker

from jwst.associations.lib import (rules_level2b, rules_level3)
RegistryMarker.mark(rules_level2b)
RegistryMarker.mark(rules_level3)
