# To avoid relative imports and mimic actual usage
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '../../../../..'))

from jwst.associations import Association
from jwst.associations.association import finalize as general_asn_finalize
from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.constraint import ConstraintTrue
from jwst.associations.lib.dms_base import DMSBaseMixin


@RegistryMarker.rule
class Asn_DMS_Base(DMSBaseMixin, Association):
    """Basic DMS rule"""

    def __init__(self, version_id=None):
        self.constraints = ConstraintTrue()
        super(Asn_DMS_Base, self).__init__(version_id=version_id)
        self.data['members'] = list()

    def make_member(self, item):
        return item

    def _add(self, item):
        self.data['members'].append(item)

    def finalize(self):
        """Peform finalization steps"""
        return [self]


# Use the generic finalization
RegistryMarker.callback('finalize')(general_asn_finalize)


class Utility:
    """Should not be part of the utilities"""

    @staticmethod
    def not_valid_function():
        """Should not be part of the utilities"""


@RegistryMarker.utility
class ValidUtility:
    """Yea, valid!"""

    @staticmethod
    def valid_function():
        "yes, i'm good"
