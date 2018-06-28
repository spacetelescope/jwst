# To avoid relative imports and mimic actual usage
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '../../../..'))

from associations import Association
from associations.lib.dms_base import DMSBaseMixin


class Asn_DMS_Base(DMSBaseMixin, Association):
    """Basic DMS rule"""

    def __init__(self, version_id=None):
        super(Asn_DMS_Base, self).__init__(version_id=version_id)
        self.data['members'] = list()

    def make_member(self, item):
        return item

    def _add(self, item):
        self.data['members'].append(item)
