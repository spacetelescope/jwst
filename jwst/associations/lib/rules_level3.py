"""Association Definitions: DMS Level3 product associations
"""
import logging

from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.dms_base import (Constraint_TargetAcq, Constraint_TSO, nrsfss_valid_detector, nrsifu_valid_detector)
from jwst.associations.lib.dms_base import (nrccoron_valid_detector)
from jwst.associations.lib.process_list import ListCategory
from jwst.associations.lib.rules_level3_base import *
from jwst.associations.lib.rules_level3_base import (
    dms_product_name_sources, dms_product_name_noopt, dms_product_name_coronimage,
    format_product
)

__all__ = [
    'Asn_Lv3ACQ_Reprocess',
    'Asn_Lv3AMI',
    'Asn_Lv3Image',
    'Asn_Lv3ImageBackground',
    'Asn_Lv3MIRCoron',
    'Asn_Lv3MIRMRS',
    'Asn_Lv3MIRMRSBackground',
    'Asn_Lv3NRCCoron',
    'Asn_Lv3NRCCoronImage',
    'Asn_Lv3NRSFSS',
    'Asn_Lv3NRSIFU',
    'Asn_Lv3NRSIFUBackground',
    'Asn_Lv3SlitlessSpectral',
    'Asn_Lv3SpecAux',
    'Asn_Lv3SpectralSource',
    'Asn_Lv3SpectralTarget',
    'Asn_Lv3TSO',
    'Asn_Lv3WFSCMB',
    'Asn_Lv3WFSSNIS',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------
# Start of the User-level rules
# --------------------------------
@RegistryMarker.rule
class Asn_Lv3ACQ_Reprocess(DMS_Level3_Base):
    """Level 3 Gather Target Acquisitions

    Characteristics:
        - Association type: Not applicable
        - Pipeline: Not applicable
        - Used to populate other related associations

    Notes
    -----
    For first loop, simply send acquisitions and confirms back.
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_TargetAcq(),
            SimpleConstraint(
                name='force_fail',
                test=lambda x, y: False,
                value='anything but None',
                reprocess_on_fail=True,
                work_over=ListCategory.NONSCIENCE,
                reprocess_rules=[]
            )
        ])

        super(Asn_Lv3ACQ_Reprocess, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv3AMI(AsnMixin_Science):
    """Level 3 Aperture Mask Interferometry Association

    Characteristics:
        - Association type: ``ami3``
        - Pipeline: ``calwebb_ami3``
        - Gather science and related PSF exposures

    Notes
    -----
    AMI is nearly completely defined by the association candidates
    produced by APT.
    Tracking Issues:

        - `github #310 <https://github.com/STScI-JWST/jwst/issues/310>`_
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
                force_reprocess=ListCategory.EXISTING,
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
        super(Asn_Lv3AMI, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'ami3'
        super(Asn_Lv3AMI, self)._init_hook(item)


@RegistryMarker.rule
class Asn_Lv3Image(AsnMixin_Science):
    """Level 3 Science Image Association

    Characteristics:
        - Association type: ``image3``
        - Pipeline: ``calwebb_image3``
        - Non-TSO
        - Non-WFS&C
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Optical_Path(),
            Constraint_Target(association=self),
            Constraint_Image(),
            DMSAttrConstraint(
                name='wfsvisit',
                sources=['visitype'],
                value='((?!wfsc).)*',
                required=False
            ),
            Constraint(
                [
                    DMSAttrConstraint(
                        name='bkgdtarg',
                        sources=['bkgdtarg'],
                    ),
                    Constraint_TSO()
                ],
                reduce=Constraint.notany
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv3Image, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image3'
        super(Asn_Lv3Image, self)._init_hook(item)


@RegistryMarker.rule
class Asn_Lv3ImageBackground(AsnMixin_AuxData, AsnMixin_Science):
    """Level 3 Background Image Association

    Characteristics:
        - Association type: ``image3``
        - Pipeline: ``calwebb_image3``
        - Non-TSO
        - Non-WFS&C
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Optical_Path(),
            Constraint_Target(association=self),
            Constraint_Image(),
            DMSAttrConstraint(
                name='wfsvisit',
                sources=['visitype'],
                value='((?!wfsc).)*',
                required=False
            ),
            DMSAttrConstraint(
                name='bkgdtarg',
                sources=['bkgdtarg'],
            ),
            Constraint(
                [Constraint_TSO()],
                reduce=Constraint.notany
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv3ImageBackground, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image3'
        super(Asn_Lv3ImageBackground, self)._init_hook(item)


@RegistryMarker.rule
class Asn_Lv3MIRCoron(AsnMixin_Coronagraphy, AsnMixin_Science):
    """Level 3 Coronagraphy Association

    Characteristics:
        - Association type: ``coron3``
        - Pipeline: ``calwebb_coron3``
        - MIRI Coronagraphy
        - Gather science and related PSF exposures

    Notes
    -----
    Coronagraphy is nearly completely defined by the association candidates
    produced by APT.
    Tracking Issues:

        - `github #311 <https://github.com/STScI-JWST/jwst/issues/311>`_
        - `JP-3219 <https://jira.stsci.edu/browse/JP-3219>`_
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint(
            [
                Constraint_Optical_Path(),
                DMSAttrConstraint(
                    name='exp_type',
                    sources=['exp_type'],
                    value='mir_lyot|mir_4qpm',
                ),
                DMSAttrConstraint(
                    name='target',
                    sources=['targetid'],
                    onlyif=lambda item: self.get_exposure_type(item) == 'science',
                    force_reprocess=ListCategory.EXISTING,
                    only_on_match=True,
                ),
                Constraint(
                    [DMSAttrConstraint(
                        name='bkgdtarg',
                        sources=['bkgdtarg'],
                        force_unique=False,
                    )],
                    reduce=Constraint.notany
                ),
            ],
            name='asn_coron'
        )

        # Check and continue initialization.
        super(Asn_Lv3MIRCoron, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv3MIRMRS(AsnMixin_Spectrum):
    """Level 3 MIRI MRS Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - Just MIRI MRS
        - optical path determined by calibration
        - Cannot be TSO
        - Must have pattern type defined
    """

    def __init__(self, *args, **kwargs):
        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(association=self),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='mir_mrs',
            ),
            Constraint(
                [
                    Constraint_TSO(),
                ],
                reduce=Constraint.notany
            ),
        ])

        # Check and continue initialization.
        super(Asn_Lv3MIRMRS, self).__init__(*args, **kwargs)

    @property
    def dms_product_name(self):
        return dms_product_name_noopt(self)


@RegistryMarker.rule
class Asn_Lv3MIRMRSBackground(AsnMixin_AuxData, AsnMixin_Spectrum):
    """Level 3 MIRI MRS Association Auxiliary data

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - Just MIRI MRS
        - optical path determined by calibration
        - Cannot be TSO
        - Must have pattern type defined
    """

    def __init__(self, *args, **kwargs):
        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='mir_mrs',
            ),
            Constraint(
                [
                    Constraint_TSO(),
                ],
                reduce=Constraint.notany
            ),
            DMSAttrConstraint(
                name='bkgdtarg',
                sources=['bkgdtarg'],
                value='T',
            ),
        ])

        # Check and continue initialization.
        super(Asn_Lv3MIRMRSBackground, self).__init__(*args, **kwargs)

    @property
    def dms_product_name(self):
        return dms_product_name_noopt(self)


@RegistryMarker.rule
class Asn_Lv3NRCCoron(AsnMixin_Coronagraphy, AsnMixin_Science):
    """Level 3 Coronagraphy Association

    Characteristics:
        - Association type: ``coron3``
        - Pipeline: ``calwebb_coron3``
        - Gather science and related PSF exposures
        - Exclude "extra" NIRCam detectors that don't have target on them

    Notes
    -----
    Coronagraphy is nearly completely defined by the association candidates
    produced by APT.
    Tracking Issues:

        - `github #311 <https://github.com/STScI-JWST/jwst/issues/311>`_
        - `JP-3219 <https://jira.stsci.edu/browse/JP-3219>`_
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint(
            [
                Constraint_Optical_Path(),
                DMSAttrConstraint(
                    name='exp_type',
                    sources=['exp_type'],
                    value=('nrc_coron'),
                ),
                DMSAttrConstraint(
                    name='target',
                    sources=['targetid'],
                    onlyif=lambda item: self.get_exposure_type(item) == 'science',
                    force_reprocess=ListCategory.EXISTING,
                    only_on_match=True,
                ),
                Constraint(
                    [DMSAttrConstraint(
                        name='bkgdtarg',
                        sources=['bkgdtarg'],
                        force_unique=False,
                    )],
                    reduce=Constraint.notany
                ),
                SimpleConstraint(
                    value=True,
                    test=lambda value, item: nrccoron_valid_detector(item),
                    force_unique=False
                ),
            ],
            name='asn_coron'
        )

        # Check and continue initialization.
        super(Asn_Lv3NRCCoron, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv3NRCCoronImage(AsnMixin_Science):
    """Level 3 Coronagraphy Association handled as regular imaging

    Characteristics:
        - Association type: ``image3``
        - Pipeline: ``calwebb_image3``
        - Gather science exposures only, no psf exposures
        - Only include NRC SW images taken in full-frame

    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint(
            [
                Constraint_Optical_Path(),
                DMSAttrConstraint(
                    name='exp_type',
                    sources=['exp_type'],
                    value=('nrc_coron'),
                ),
                DMSAttrConstraint(
                    name='target',
                    sources=['targetid'],
                    onlyif=lambda item: self.get_exposure_type(item) == 'science',
                    force_reprocess=ListCategory.EXISTING,
                    only_on_match=True,
                ),
                Constraint(
                    [DMSAttrConstraint(
                        name='bkgdtarg',
                        sources=['bkgdtarg'],
                        force_unique=False,
                    ),
                    DMSAttrConstraint(
                        name='is_psf',
                        sources=['is_psf'],
                        value = ('T')
                    )],
                    reduce=Constraint.notany
                ),  
                DMSAttrConstraint(
                    name='channel',
                    sources=['channel'],
                    value=('short'),
                ),
                DMSAttrConstraint(
                    name='subarray',
                    sources=['subarray'],
                    value=('full'),
                ),
            ],
        )

        # Check and continue initialization.
        super(Asn_Lv3NRCCoronImage, self).__init__(*args, **kwargs)

    @property
    def dms_product_name(self):
        return dms_product_name_coronimage(self)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image3'
        super(Asn_Lv3NRCCoronImage, self)._init_hook(item)

    def is_item_coron(self, item):
        """Override to ignore coronographic designation

        Coronagraphic data is to be processed both as coronagraphic
        (by default), but also as just plain imaging. Coronagraphic
        data is processed using the Asn_Lv3Coron rule. This rule
        will handle the creation of the image version. It causes
        the input members to be of type "cal", instead of "calints".
        """
        return False


@RegistryMarker.rule
class Asn_Lv3NRSFSS(AsnMixin_Spectrum):
    """Level 3 NIRSpec Fixed-slit Science

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - NIRSpec Fixed-slit Science
        - Non-TSO
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint(
                [Constraint_TSO()],
                reduce=Constraint.notany
            ),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrs_autoflat'
                    '|nrs_autowave'
                    '|nrs_fixedslit'
                ),
                force_unique=False
            ),
            SimpleConstraint(
                value=True,
                test=lambda value, item: nrsfss_valid_detector(item),
                force_unique=False
            ),
            Constraint_Optical_Path(),
            Constraint_Target(association=self),
        ])

        # Check and continue initialization.
        super(Asn_Lv3NRSFSS, self).__init__(*args, **kwargs)

    @property
    def dms_product_name(self):
        return dms_product_name_sources(self)


@RegistryMarker.rule
class Asn_Lv3NRSIFU(AsnMixin_Spectrum):
    """Level 3 IFU gratings Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - optical path determined by calibration
    """

    def __init__(self, *args, **kwargs):
        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(association=self),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nrs_autowave'
                    '|nrs_ifu'
                ),
                force_unique=False
            ),
            SimpleConstraint(
                value=True,
                test=lambda value, item: nrsifu_valid_detector(item),
                force_unique=False
            ),
            DMSAttrConstraint(
                name='patttype',
                sources=['patttype'],
                required=True
            ),
            Constraint(
                [
                    Constraint_TSO(),
                ],
                reduce=Constraint.notany
            ),
            DMSAttrConstraint(
                    name='opt_elem',
                    sources=['grating'],
            )
        ])

        # Check and continue initialization.
        super(Asn_Lv3NRSIFU, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv3NRSIFUBackground(AsnMixin_AuxData, AsnMixin_Spectrum):

    """Level 3 Spectral Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
    """
    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(association=self),
            Constraint(
                [
                    Constraint_TSO(),
                ],
                reduce=Constraint.notany
            ),
            DMSAttrConstraint(
                name='bkgdtarg',
                sources=['bkgdtarg'],
                value='T',
            ),
            DMSAttrConstraint(
                name='allowed_bkgdtarg',
                sources=['exp_type'],
                value='nrs_ifu',
            ),
            SimpleConstraint(
                value=True,
                test=lambda value, item: nrsifu_valid_detector(item),
                force_unique=False
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['grating'],
                force_unique=True,
            ),
        ])

        # Check and continue initialization.
        super(Asn_Lv3NRSIFUBackground, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv3SlitlessSpectral(AsnMixin_Spectrum):
    """Level 3 slitless, target-based or single-object spectrographic Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - Single target
        - Non-TSO
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint(
                [Constraint_TSO()],
                reduce=Constraint.notany
            ),
            Constraint_Optical_Path(),
            Constraint_Target(association=self),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'nis_soss'
                ),
                force_unique=False
            ),
            Constraint(
                [
                    DMSAttrConstraint(
                        name='patttype_spectarg',
                        sources=['patttype'],
                    ),
                ],
                reduce=Constraint.notany
            ),
            # Constraint to prevent calibration data from level 3 processing
            Constraint(
                [
                    DMSAttrConstraint(
                        name='restricted_slitless',
                        sources=['exp_type'],
                        value = ('mir_lrs-slitless')
                    ),
                    DMSAttrConstraint(
                        name='tso_obs',
                        sources=['tso_visit'],
                        value = ('T')
                    ),
                ],
                reduce=Constraint.notany
            )
        ])

        # Check and continue initialization.
        super(Asn_Lv3SlitlessSpectral, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv3SpecAux(AsnMixin_AuxData, AsnMixin_Spectrum):

    """Level 3 Spectral Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
    """
    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(association=self),
            Constraint(
                [
                    Constraint_TSO(),
                ],
                reduce=Constraint.notany
            ),
            DMSAttrConstraint(
                name='bkgdtarg',
                sources=['bkgdtarg'],
                value='T',
            ),
            DMSAttrConstraint(
                name='allowed_bkgdtarg',
                sources=['exp_type'],
                value='mir_lrs-fixedslit|nrs_fixedslit',
            ),
            Constraint_Optical_Path(),
        ])

        # Check and continue initialization.
        super(Asn_Lv3SpecAux, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv3SpectralSource(AsnMixin_Spectrum):
    """Level 3 slit-like, multi-object spectrographic Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - Multi-object
        - Non-TSO
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint(
                [Constraint_TSO()],
                reduce=Constraint.notany
            ),
            Constraint_Optical_Path(),
            Constraint_Target(association=self),
            Constraint(
                [
                    DMSAttrConstraint(
                        name='exp_type',
                        sources=['exp_type'],
                        value=(
                            'nrc_wfss'
                            '|nrs_autoflat'
                            '|nrs_autowave'
                        ),
                        force_unique=False
                    ),
                    Constraint_MSA()
                ],
                reduce=Constraint.any
            )
        ])

        # Check and continue initialization.
        super(Asn_Lv3SpectralSource, self).__init__(*args, **kwargs)

    @property
    def dms_product_name(self):
        return dms_product_name_sources(self)


@RegistryMarker.rule
class Asn_Lv3SpectralTarget(AsnMixin_Spectrum):
    """Level 3 slit-like, target-based or single-object spectrographic Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - Single target
        - Non-TSO
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint(
                [Constraint_TSO()],
                reduce=Constraint.notany
            ),
            Constraint_Optical_Path(),
            Constraint_Target(association=self),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'mir_lrs-fixedslit'
                    '|nis_soss'
                ),
                force_unique=False
            ),
            Constraint(
                [
                    DMSAttrConstraint(
                        name='patttype_spectarg',
                        sources=['patttype'],
                        value='2-point-nod|4-point-nod|along-slit-nod',
                    ),
                ],
                reduce=Constraint.any
            )
        ])

        # Check and continue initialization.
        super(Asn_Lv3SpectralTarget, self).__init__(*args, **kwargs)

    def finalize(self):
        """Finalize association

        For NRS Fixed-slit, finalization means creating new members for the
        background nods.

        Returns
        -------
        associations: [association[, ...]] or None
            List of fully-qualified associations that this association
            represents.
            `None` if a complete association cannot be produced.

        """
        if self.is_valid:
            return self.make_fixedslit_bkg()

        return None


@RegistryMarker.rule
class Asn_Lv3TSO(AsnMixin_Science):
    """Level 3 Time-Series Association

    Characteristics:
        - Association type: ``tso3``
        - Pipeline: ``calwebb_tso3``
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(association=self),
            Constraint_Optical_Path(),
            Constraint_TSO(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
            ),
            # Don't allow IFU exposures in tso3
            Constraint(
                [
                    Constraint_IFU(),
                ],
                reduce=Constraint.notany
            ),
            # Don't allow NIRCam engineering mode
            # with PUPIL='CLEAR' in tso3
            Constraint(
                [
                    Constraint([
                        DMSAttrConstraint(
                            name='restricted_grism',
                            sources=['exp_type'],
                            value = ('nrc_tsgrism')
                        ),
                        DMSAttrConstraint(
                            name='grism_clear',
                            sources=['pupil'],
                            value = ('clear')
                        ),
                    ]),
                    Constraint([
                        DMSAttrConstraint(
                            name='restricted_ts',
                            sources=['exp_type'],
                            value = 'nrc_tsgrism'
                        ),
                        DMSAttrConstraint(
                            name='module',
                            sources=['detector'],
                            value='nrcblong'
                        ),
                    ]),
                ],
                reduce=Constraint.notany
            )
        ])

        super(Asn_Lv3TSO, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'tso3'
        super(Asn_Lv3TSO, self)._init_hook(item)


@RegistryMarker.rule
class Asn_Lv3WFSCMB(AsnMixin_Science):
    """Level 3 Wavefront Control & Sensing Association

    For coarse and fine phasing, dither pairs need to
    be associated to be combined.  The optical path
    is assumed to be equivalent within an activity.

    Characteristics:
        - Association type: ``wfs-image3``
        - Pipeline: ``calwebb_wfs-image3``
        - Coarse and fine phasing dithers
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Optical_Path(),
            Constraint_Target(association=self),
            Constraint_Image(),
            DMSAttrConstraint(
                name='patttype',
                sources=['patttype'],
                value='wfsc'
            ),
            DMSAttrConstraint(
                name='detector',
                sources=['detector']
            ),
            DMSAttrConstraint(
                name='obs_id',
                sources=['obs_id']
            ),
            DMSAttrConstraint(
                name='act_id',
                sources=['act_id']
            ),
            Constraint(
                [
                    DMSAttrConstraint(
                        name='dms_note',
                        sources=['dms_note'],
                        value='wfsc_los_jitter'
                    ),
                ],
                reduce=Constraint.notany,
            ),
        ])

        # Only valid if two members exist and candidate is not a GROUP.
        self.validity.update({
            'has_pair': {
                'validated': False,
                'check': self._has_pair
            },
            'is_not_group': {
                'validated': False,
                'check': self._validate_candidates
            }
        })

        super(Asn_Lv3WFSCMB, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'wfs-image3'
        super(Asn_Lv3WFSCMB, self)._init_hook(item)

    @property
    def dms_product_name(self):
        """Define product name

        Modification is to append the `expspcin` value
        after the calibration suffix.
        """
        product_name_format = '{existing}-{detector}_{suffix}-{expspcin}'

        existing = super().dms_product_name

        product_name = format_product(
            product_name_format,
            existing=existing,
            detector=self.constraints['detector'].value,
            expspcin=self.constraints['act_id'].value
        )

        return product_name.lower()

    def _has_pair(self, entry=None):
        """Check if current product has two members

        If `entry` is given, it is counted as one of the
        members. If not, the existing member list is only
        accounted for.

        Parameters
        ----------
        entry : dict or None
            New entry to add to member list.

        Returns
        -------
        bool
            True if there are two members.
        """
        if entry is None:
            count = 2
        else:
            count = 1

        return len(self.current_product['members']) == count

    def _validate_candidates(self, member):
        """Disallow GROUP candidates

        Parameters
        ----------
        member : Member
            Member being added. Ignored.

        Returns
        -------
        False if candidate is GROUP.
        True otherwise.
        """

        # If a group candidate, reject.
        if self.acid.type.lower() == 'group':
            return False

        return True

@RegistryMarker.rule
class Asn_Lv3WFSSNIS(AsnMixin_Spectrum):
    """Level 3 WFSS/Grism Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - Gather all grism exposures
    """

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(association=self),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='nis_wfss',
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['filter'],
                value='gr150r|gr150c',
                force_unique=True,
            ),
            DMSAttrConstraint(
                name='opt_elem2',
                sources=['pupil'],
            ),
        ])

        # Check and continue initialization.
        super(Asn_Lv3WFSSNIS, self).__init__(*args, **kwargs)

    @property
    def dms_product_name(self):
        return dms_product_name_sources(self)
