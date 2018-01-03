"""Association Definitions: DMS Level3 product associations
"""
import logging

from jwst.associations.lib.rules_level3_base import *

__all__ = [
    'Asn_Image',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# --------------------------------
# Start of the User-level rules
# --------------------------------


# ----------------------------------
# Image associations
class Asn_Image(DMS_Level3_Base):
    """Non-Association Candidate Dither Associations"""

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            CONSTRAINT_BASE,
            CONSTRAINT_OPTICAL_PATH,
            CONSTRAINT_TARGET,
            CONSTRAINT_IMAGE,
            AttrConstraint(
                name='wfsvisit',
                sources=['visitype'],
                value='((?!wfsc).)*'
            ),
        ])

        # Now check and continue initialization.
        super(Asn_Image, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image3'
        super(Asn_Image, self)._init_hook(item)


class Asn_WFSCMB(DMS_Level3_Base):
    """Wavefront Sensing association

    Notes
    -----
    Defined by `TRAC issue #269 <https://aeon.stsci.edu/ssb/trac/jwst/ticket/269>`_
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            CONSTRAINT_BASE,
            CONSTRAINT_OPTICAL_PATH,
            CONSTRAINT_TARGET,
            CONSTRAINT_IMAGE,
            AttrConstraint(
                name='wfsvisit',
                sources=['visitype'],
                value='.+wfsc.+',
            ),
            AttrConstraint(
                name='asn_candidate_wfs',
                sources=['asn_candidate'],
                value='.+mosaic.+',
                force_unique=True,
                is_acid=True,
                evaluate=True,
            ),
            AttrConstraint(
                name='activity_id',
                sources=['act_id']
            )
        ])

        super(Asn_WFSCMB, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'wfs'
        super(Asn_WFSCMB, self)._init_hook(item)
