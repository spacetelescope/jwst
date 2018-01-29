"""Association Definitions: DMS Level3 product associations
"""
import logging

from jwst.associations.lib.dms_base import ACQ_EXP_TYPES
from jwst.associations.lib.rules_level3_base import *

__all__ = [
    'Asn_AMI',
    'Asn_Coron',
    'Asn_Image',
    'Asn_MIRI_LRS_FIXEDSLIT',
    'Asn_MIRI_LRS_SLITLESS',
    'Asn_MIRI_IFU',
    'Asn_NRS_FIXEDSLIT',
    'Asn_NRS_IFU',
    'Asn_NRS_MSA',
    'Asn_NIS_SO_SLITLESS',
    'Asn_TSO_EXPTYPE',
    'Asn_TSO_Flag',
    'Asn_WFSCMB',
    'Asn_WFSS',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------
# Start of the User-level rules
# --------------------------------
class Asn_Image(AsnMixin_Science):
    """Non-Association Candidate Dither Associations"""

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Optical_Path(),
            Constraint_Target(),
            Constraint_Image(),
            DMSAttrConstraint(
                name='wfsvisit',
                sources=['visitype'],
                value='((?!wfsc).)*',
                required=False
            ),
        ])

        # Now check and continue initialization.
        super(Asn_Image, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image3'
        super(Asn_Image, self)._init_hook(item)


class Asn_WFSCMB(AsnMixin_Science):
    """Wavefront Sensing association

    Notes
    -----
    Defined by `TRAC issue #269 <https://aeon.stsci.edu/ssb/trac/jwst/ticket/269>`_
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Optical_Path(),
            Constraint_Target(),
            Constraint_Image(),
            DMSAttrConstraint(
                name='wfsvisit',
                sources=['visitype'],
                value='.+wfsc.+',
            ),
            DMSAttrConstraint(
                name='asn_candidate_wfs',
                sources=['asn_candidate'],
                value='.+mosaic.+',
                force_unique=True,
                is_acid=True,
                evaluate=True,
            ),
            DMSAttrConstraint(
                name='activity_id',
                sources=['act_id']
            )
        ])

        super(Asn_WFSCMB, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'wfs'
        super(Asn_WFSCMB, self)._init_hook(item)


class Asn_MIRI_LRS_FIXEDSLIT(AsnMixin_Spectrum):
    """MIRI LRS Fixed slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_NotTSO(),
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_lrs-fixedslit'
                ),
                force_unique=False,
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['filter'],
                value='p750l',
            ),
            DMSAttrConstraint(
                name='subarray',
                sources=['subarray'],
                value='full',
            )
        ])

        # Check and continue initialization.
        super(Asn_MIRI_LRS_FIXEDSLIT, self).__init__(*args, **kwargs)


class Asn_MIRI_LRS_SLITLESS(AsnMixin_Spectrum):
    """MIRI LRS Slitless"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_NotTSO(),
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_lrs-slitless'
                ),
                force_unique=False,
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['filter'],
                value='p750l',
            ),
            DMSAttrConstraint(
                name='subarray',
                sources=['subarray'],
                value='subprism',
            )
        ])

        # Check and continue initialization.
        super(Asn_MIRI_LRS_SLITLESS, self).__init__(*args, **kwargs)


class Asn_NIS_SO_SLITLESS(AsnMixin_Spectrum):
    """NIRISS Single-Object Slitless"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nis_soss'
                ),
                force_unique=False,
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['pupil'],
                value='gr700xd',
            ),
            DMSAttrConstraint(
                name='subarray',
                sources=['subarray'],
                value=(
                    'full'
                    '|substrip256'
                    '|substrip80'
                ),
            )
        ])

        # Check and continue initialization.
        super(Asn_NIS_SO_SLITLESS, self).__init__(*args, **kwargs)


class Asn_NRS_FIXEDSLIT(AsnMixin_Spectrum):
    """NIRSPEC Fixed Slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Optical_Path(),
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrs_fixedslit'
                    '|nrs_autowave'
                ),
                force_unique=False,
            ),
            DMSAttrConstraint(
                name='fixed_slit',
                sources=['fxd_slit']
            ),
            DMSAttrConstraint(
                name='subarray',
                sources=['subarray']
            ),
        ])

        # Check and continue initialization.
        super(Asn_NRS_FIXEDSLIT, self).__init__(*args, **kwargs)


class Asn_NRS_MSA(AsnMixin_Spectrum):
    """NIRSPEC MSA"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Optical_Path(),
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrs_msaspec'
                    '|nrs_autoflat'
                    '|nrs_autowave'
                ),
                force_unique=False,
            ),
        ])

        # Check and continue initialization.
        super(Asn_NRS_MSA, self).__init__(*args, **kwargs)


class Asn_MIRI_IFU(AsnMixin_Spectrum):
    """MIRI MRS (IFU)"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_mrs'
                    '|mir_flatmrs'
                ),
                force_unique=False,
            ),
        ])

        # Check and continue initialization.
        super(Asn_MIRI_IFU, self).__init__(*args, **kwargs)

    def dms_product_name(self):
        """Define product name."""
        target = self._get_target()

        instrument = self._get_instrument()

        product_name = 'jw{}-{}_{}_{}'.format(
            self.data['program'],
            self.acid.id,
            target,
            instrument
        )

        return product_name.lower()


class Asn_NRS_IFU(AsnMixin_Spectrum):
    """NIRSPEC IFU"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrs_ifu'
                    '|nrs_autowave'
                ),
                force_unique=False,
            ),
        ])

        # Check and continue initialization.
        super(Asn_NRS_IFU, self).__init__(*args, **kwargs)


class Asn_Coron(AsnMixin_Science):
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
        self.constraints = Constraint(
            [
                Constraint_Optical_Path(),
                DMSAttrConstraint(
                    name='exp_type',
                    sources=['exp_type'],
                    value=(
                        'nrc_coron'
                        '|mir_lyot'
                        '|mir_4qpm'
                    ),
                ),
                DMSAttrConstraint(
                    name='target',
                    sources=['targetid'],
                    onlyif=lambda item: self.get_exposure_type(item) == 'science',
                    force_reprocess=ProcessList.EXISTING,
                    only_on_match=True,
                ),
            ],
            name='asn_coron'
        )

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


class Asn_AMI(AsnMixin_Science):
    """Aperture Mask Interferometry
    Notes
    -----
    AMI is nearly completely defined by the association candidates
    produced by APT.
    Tracking Issues:
    - `github #310 <https://github.com/STScI-JWST/jwst/issues/310>`
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Optical_Path(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nis_ami'
                ),
            ),
            DMSAttrConstraint(
                name='target',
                sources=['targetid'],
                onlyif=lambda item: self.get_exposure_type(item) == 'science',
                force_reprocess=ProcessList.EXISTING,
                only_on_match=True,
            ),
        ])

        # Check and continue initialization.
        super(Asn_AMI, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'ami3'
        super(Asn_AMI, self)._init_hook(item)


class Asn_WFSS(AsnMixin_Spectrum):
    """WFSS/Grism modes"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='nis_wfss',
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['filter'],
            ),
            DMSAttrConstraint(
                name='opt_elem2',
                sources=['grating'],
                value='gr150r|gr150c',
                force_unique=False,
            ),
        ])

        # Check and continue initialization.
        super(Asn_WFSS, self).__init__(*args, **kwargs)


class Asn_TSO_Flag(AsnMixin_Science):
    """Time-Series observations"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(),
            Constraint_Optical_Path(),
            DMSAttrConstraint(
                name='is_tso',
                sources=['tsovisit'],
                value='t',
            ),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type']
            )
        ])

        super(Asn_TSO_Flag, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'tso3'
        super(Asn_TSO_Flag, self)._init_hook(item)


class Asn_TSO_EXPTYPE(AsnMixin_Science):
    """Time-Series observations"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(),
            Constraint_Optical_Path(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_lrs-slitless'
                    '|nis_soss'
                    '|nrc_tsimage'
                    '|nrc_tsgrism'
                    '|nrs_bota'
                    '|nrs_brightobj'
                ),
            ),
            DMSAttrConstraint(
                name='no_tso_flag',
                sources=['tsovisit'],
                required=False,
                force_undefined=True
            )
        ])

        super(Asn_TSO_EXPTYPE, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'tso3'
        super(Asn_TSO_EXPTYPE, self)._init_hook(item)


class Asn_ACQ_Reprocess(DMS_Level3_Base):
    """For first loop, simply send acquisitions and confirms back"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
                DMSAttrConstraint(
                    sources=['exp_type'],
                    value='|'.join(ACQ_EXP_TYPES),
                    force_unique=False
                ),
                SimpleConstraint(
                    name='force_fail',
                    test=lambda x, y: False,
                    value='anything but None',
                    force_reprocess=ProcessList.NONSCIENCE
                )
            ])

        super(Asn_ACQ_Reprocess, self).__init__(*args, **kwargs)
