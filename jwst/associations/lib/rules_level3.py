"""Association Definitions: DMS Level3 product associations
"""
import logging

from jwst.associations.lib.rules_level3_base import *

__all__ = [
    'Asn_Image',
    'Asn_MIRI_LRS_FIXEDSLIT',
    'Asn_MIRI_LRS_SLITLESS',
    'Asn_NIS_SO_SLITLESS',
    'Asn_WFSCMB',
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
            LV3AttrConstraint(
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
            LV3AttrConstraint(
                name='wfsvisit',
                sources=['visitype'],
                value='.+wfsc.+',
            ),
            LV3AttrConstraint(
                name='asn_candidate_wfs',
                sources=['asn_candidate'],
                value='.+mosaic.+',
                force_unique=True,
                is_acid=True,
                evaluate=True,
            ),
            LV3AttrConstraint(
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
            CONSTRAINT_BASE,
            CONSTRAINT_NOTTSO,
            CONSTRAINT_TARGET,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_lrs-fixedslit'
                    '|mir_tacq'
                ),
                force_unique=False,
            ),
            LV3AttrConstraint(
                name='opt_elem',
                sources=['filter'],
                value='p750l',
            ),
            LV3AttrConstraint(
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
            CONSTRAINT_BASE,
            CONSTRAINT_NOTTSO,
            CONSTRAINT_TARGET,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_lrs-slitless'
                    '|mir_tacq'
                ),
                force_unique=False,
            ),
            LV3AttrConstraint(
                name='opt_elem',
                sources=['filter'],
                value='p750l',
            ),
            LV3AttrConstraint(
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
            CONSTRAINT_BASE,
            CONSTRAINT_TARGET,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nis_soss'
                    '|nis_tacq'
                    '|nis_tacnfrm'
                ),
                force_unique=False,
            ),
            LV3AttrConstraint(
                name='opt_elem',
                sources=['pupil'],
                value='gr700xd',
            ),
            LV3AttrConstraint(
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
            CONSTRAINT_BASE,
            CONSTRAINT_OPTICAL_PATH,
            CONSTRAINT_TARGET,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrs_fixedslit'
                    '|nrs_autowave'
                    '|nrs_confirm'
                    '|nrs_taconfirm'
                    '|nrs_tacq'
                    '|nrs_taslit'
                ),
                force_unique=False,
            ),
            LV3AttrConstraint(
                name='fixed_slit',
                sources=['fxd_slit']
            ),
            LV3AttrConstraint(
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
            CONSTRAINT_BASE,
            CONSTRAINT_OPTICAL_PATH,
            CONSTRAINT_TARGET,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrs_msaspec'
                    '|nrs_autoflat'
                    '|nrs_autowave'
                    '|nrs_confirm'
                    '|nrs_taslit'
                    '|nrs_taconfirm'
                    '|nrs_tacq'
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
            CONSTRAINT_BASE,
            CONSTRAINT_TARGET,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_mrs'
                    '|mir_flatmrs'
                    '|mir_tacq'
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
            CONSTRAINT_BASE,
            CONSTRAINT_TARGET,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrs_ifu'
                    '|nrs_autowave'
                    '|nrs_confirm'
                    '|nrs_taconfirm'
                    '|nrs_tacq'
                    '|nrs_taslit'
                ),
                force_unique=False,
            ),
        ])

        # Check and continue initialization.
        super(Asn_NRS_IFU, self).__init__(*args, **kwargs)


class Asn_Coron(DMS_Level3_Base):
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
        self.constraints = Constraint([
            CONSTRAINT_BASE,
            CONSTRAINT_OPTICAL_PATH,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrc_coron'
                    '|mir_lyot'
                    '|mir_4qpm'
                ),
            ),
            LV3AttrConstraint(
                name='target',
                sources=['targetid'],
                onlyif=lambda item: self.get_exposure_type(item) == 'science',
                force_reprocess=ProcessList.EXISTING,
                only_on_match=True,
            ),
        ])

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


class Asn_AMI(DMS_Level3_Base):
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
            CONSTRAINT_BASE,
            CONSTRAINT_OPTICAL_PATH,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nis_ami'
                    '|nis_taconfirm'
                    '|nis_tacq'
                ),
            ),
            LV3AttrConstraint(
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
            CONSTRAINT_BASE,
            CONSTRAINT_TARGET,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='nis_wfss',
            ),
            LV3AttrConstraint(
                name='opt_elem',
                sources=['filter'],
            ),
            LV3AttrConstraint(
                name='opt_elem2',
                sources=['grating'],
                value='gr150r|gr150c',
                force_unique=False,
            ),
        ])

        # Check and continue initialization.
        super(Asn_WFSS, self).__init__(*args, **kwargs)



class Asn_TSO_Flag(DMS_Level3_Base):
    """Time-Series observations"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            CONSTRAINT_BASE,
            CONSTRAINT_TARGET,
            CONSTRAINT_OPTICAL_PATH,
            LV3AttrConstraint(
                name='is_tso',
                sources=['tsovisit'],
                value='t',
            ),
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type']
            )
        ])

        super(Asn_TSO_Flag, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'tso3'
        super(Asn_TSO_Flag, self)._init_hook(item)


class Asn_TSO_EXPTYPE(DMS_Level3_Base):
    """Time-Series observations"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            CONSTRAINT_BASE,
            CONSTRAINT_TARGET,
            CONSTRAINT_OPTICAL_PATH,
            LV3AttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_lrs-slitless'
                    '|nis_soss'
                    '|nis_taconfirm'
                    '|nis_tacq'
                    '|nrc_tsimage'
                    '|nrc_tsgrism'
                    '|nrs_bota'
                    '|nrs_brightobj'
                    '|nrs_taconfirm'
                ),
            ),
            LV3AttrConstraint(
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
