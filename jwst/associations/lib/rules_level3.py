"""Association Definitions: DMS Level3 product associations
"""
import logging

from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.dms_base import (ACQ_EXP_TYPES, Constraint_TSO)
from jwst.associations.lib.rules_level3_base import *
from jwst.associations.lib.rules_level3_base import (
    dms_product_name_sources,
    format_product
)

__all__ = [
    'Asn_AMI',
    'Asn_ACQ_Reprocess',
    'Asn_Coron',
    'Asn_IFU',
    'Asn_Lv3SpecAux',
    'Asn_Image',
    'Asn_SpectralSource',
    'Asn_SpectralTarget',
    'Asn_TSO',
    'Asn_WFSCMB',
    'Asn_WFSS_NIS',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------
# Start of the User-level rules
# --------------------------------
@RegistryMarker.rule
class Asn_Image(AsnMixin_Science):
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
                [Constraint_TSO()],
                reduce=Constraint.notany
            )
        ])

        # Now check and continue initialization.
        super(Asn_Image, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image3'
        super(Asn_Image, self)._init_hook(item)


@RegistryMarker.rule
class Asn_WFSCMB(AsnMixin_Science):
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
            )
        ])

        # Only valid if two members exist.
        self.validity.update({
            'has_pair': {
                'validated': False,
                'check': self._has_pair
            }
        })

        super(Asn_WFSCMB, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'wfs-image3'
        super(Asn_WFSCMB, self)._init_hook(item)

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


@RegistryMarker.rule
class Asn_SpectralTarget(AsnMixin_Spectrum):
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
                    '|mir_lrs-slitless'
                    '|nis_soss'
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
            )
        ])

        # Check and continue initialization.
        super(Asn_SpectralTarget, self).__init__(*args, **kwargs)

    def finalize(self):
        """Finalize assocation

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
        else:
            return None

@RegistryMarker.rule
class Asn_SpectralSource(AsnMixin_Spectrum):
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
                            '|nrs_fixedslit'
                        ),
                        force_unique=False
                    ),
                    Constraint_MSA()
                ],
                reduce=Constraint.any
            )
        ])

        # Check and continue initialization.
        super(Asn_SpectralSource, self).__init__(*args, **kwargs)

    @property
    def dms_product_name(self):
        return dms_product_name_sources(self)


@RegistryMarker.rule
class Asn_IFU(AsnMixin_Spectrum):
    """Level 3 IFU Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
        - optical path determined by calibration
    """

    def __init__(self, *args, **kwargs):
        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(association=self),
            Constraint_IFU(),
            Constraint(
                [
                    Constraint_TSO(),
                    DMSAttrConstraint(
                        name='patttype',
                        sources=['patttype'],
                        value=['2_point|4_point_nod|along_slit_nod'],
                    )
                ],
                reduce=Constraint.notany
            )        ])
        # Check and continue initialization.
        super(Asn_IFU, self).__init__(*args, **kwargs)

    @property
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

@RegistryMarker.rule
class Asn_Lv3SpecAux(AsnMixin_AuxData, AsnMixin_BkgScience):

    """Level 3 Spectral Association

    Characteristics:
        - Association type: ``spec3``
        - Pipeline: ``calwebb_spec3``
    """
    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.constraints = Constraint([
            Constraint_Target(association=self),
            Constraint_IFU(),
            Constraint(
                [
                    Constraint_TSO(),
                ],
                reduce=Constraint.notany
            ),
            Constraint(
                [
                    DMSAttrConstraint(
                        name='bkgdtarg',
                        sources=['bkgdtarg'],
                        value=['T'],)
                ],
                reduce=Constraint.any
                                    ),
            Constraint(
                [
                    DMSAttrConstraint(
                        name='allowed_bkgdtarg',
                        sources=['exp_type'],
                        value=['mir_mrs','nrs_ifu'],)
                ],
            reduce=Constraint.any
                    ),
                ])

        # Check and continue initialization.
        super(Asn_Lv3SpecAux, self).__init__(*args, **kwargs)

@RegistryMarker.rule
class Asn_Coron(AsnMixin_Science):
    """Level 3 Coronography Association

    Characteristics:
        - Association type: ``coron3``
        - Pipeline: ``calwebb_coron3``
        - Gather science and related PSF exposures

    Notes
    -----
    Coronography is nearly completely defined by the association candidates
    produced by APT.
    Tracking Issues:

        - `github #311 <https://github.com/STScI-JWST/jwst/issues/311>`_
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


@RegistryMarker.rule
class Asn_AMI(AsnMixin_Science):
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
        super(Asn_AMI, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'ami3'
        super(Asn_AMI, self)._init_hook(item)


@RegistryMarker.rule
class Asn_WFSS_NIS(AsnMixin_Spectrum):
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
                force_unique=False,
            ),
            DMSAttrConstraint(
                name='opt_elem2',
                sources=['pupil'],
            ),
        ])

        # Check and continue initialization.
        super(Asn_WFSS_NIS, self).__init__(*args, **kwargs)

    @property
    def dms_product_name(self):
        return dms_product_name_sources(self)


@RegistryMarker.rule
class Asn_TSO(AsnMixin_Science):
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
        ])

        super(Asn_TSO, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'tso3'
        super(Asn_TSO, self)._init_hook(item)


@RegistryMarker.rule
class Asn_ACQ_Reprocess(DMS_Level3_Base):
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
            DMSAttrConstraint(
                sources=['exp_type'],
                value='|'.join(ACQ_EXP_TYPES),
                force_unique=False
            ),
            SimpleConstraint(
                name='force_fail',
                test=lambda x, y: False,
                value='anything but None',
                reprocess_on_fail=True,
                work_over=ProcessList.NONSCIENCE,
                reprocess_rules=[]
            )
        ])

        super(Asn_ACQ_Reprocess, self).__init__(*args, **kwargs)
