"""Association Definitions: DMS Level3 product associations
"""
import logging

from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.rules_level3_base import *

__all__ = [
    'Asn_Coron',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------
# Start of the User-level rules
# --------------------------------

@RegistryMarker.rule
class Asn_Coron(
        AsnMixin_OpticalPath,
        AsnMixin_Base
):
    """Coronography

    Notes
    -----

    Coronography is nearly completely defined by the association candidates
    produced by APT.

    Tracking Issues:

    - `github #311 <https://github.com/STScI-JWST/jwst/issues/311>`
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': (
                    'nrc_coron'
                    '|mir_lyot'
                    '|mir_4qpm'
                ),
                'inputs': ['exp_type'],
                'force_unique': True,
            },
            'target': {
                'value': None,
                'inputs': ['targetid'],
                'onlyif': lambda item: self.get_exposure_type(item) == 'science'
            }
        })

        # PSF is required
        self.validity.update({
            'has_psf': {
                'validated': False,
                'check': lambda entry: entry['exptype'] == 'psf'
            }
        })

        # Check and continue initialization.
        super(Asn_Coron, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'coron3'
        super(Asn_Coron, self)._init_hook(item)
