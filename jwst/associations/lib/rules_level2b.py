"""Association Definitions: DMS Level2b product associations
"""
from collections import deque
import logging

from jwst.associations.exceptions import AssociationNotValidError
from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.constraint import (Constraint, SimpleConstraint)
from jwst.associations.lib.dms_base import (
    Constraint_TSO,
    Constraint_WFSC,
    format_list,
    item_getattr,
    nrsfss_valid_detector,
    nrsifu_valid_detector,
)
from jwst.associations.lib.member import Member
from jwst.associations.lib.process_list import ProcessList
from jwst.associations.lib.utilities import (getattr_from_list, getattr_from_list_nofail)
from jwst.associations.lib.rules_level2_base import *
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

__all__ = [
    'Asn_Lv2FGS',
    'Asn_Lv2Image',
    'Asn_Lv2ImageNonScience',
    'Asn_Lv2ImageSpecial',
    'Asn_Lv2ImageTSO',
    'Asn_Lv2MIRLRSFixedSlitNod',
    'Asn_Lv2NRSFSS',
    'Asn_Lv2NRSIFUNod',
    'Asn_Lv2NRSLAMPSpectral',
    'Asn_Lv2NRSMSA',
    'Asn_Lv2Spec',
    'Asn_Lv2SpecSpecial',
    'Asn_Lv2SpecTSO',
    'Asn_Lv2WFSS',
    'Asn_Lv2WFSC',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------
# Start of the User-level rules
# --------------------------------
@RegistryMarker.rule
class Asn_Lv2Image(
        AsnMixin_Lv2Image,
        DMSLevel2bBase
):
    """Level2b Non-TSO Science Image Association

    Characteristics:
        - Association type: ``image2``
        - Pipeline: ``calwebb_image2``
        - Image-based science exposures
        - Single science exposure
        - Non-TSO
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint_Image_Science(),
            Constraint_Single_Science(self.has_science),
            Constraint(
                [Constraint_TSO()],
                reduce=Constraint.notany
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv2Image, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2ImageNonScience(
        AsnMixin_Lv2Special,
        AsnMixin_Lv2Image,
        DMSLevel2bBase
):
    """Level2b Non-science Image Association

    Characteristics:
        - Association type: ``image2``
        - Pipeline: ``calwebb_image2``
        - Image-based non-science exposures, such as target acquisitions
        - Single science exposure
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Image_Nonscience(),
            Constraint_Single_Science(self.has_science),
        ])

        # Now check and continue initialization.
        super(Asn_Lv2ImageNonScience, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2ImageSpecial(
        AsnMixin_Lv2Special,
        AsnMixin_Lv2Image,
        DMSLevel2bBase
):
    """Level2b Auxiliary Science Image Association

    Characteristics:
        - Association type: ``image2``
        - Pipeline: ``calwebb_image2``
        - Image-based science exposures that are to be used as background or PSF exposures
        - Single science exposure
        - No other exposure can be part of the association
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint_Image_Science(),
            Constraint_Single_Science(self.has_science),
            Constraint_Special(),
        ])

        # Now check and continue initialization.
        super(Asn_Lv2ImageSpecial, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2ImageTSO(
        AsnMixin_Lv2Image,
        DMSLevel2bBase
):
    """Level2b Time Series Science Image Association

    Characteristics:
        - Association type: ``tso-image2``
        - Pipeline: ``calwebb_tso-image2``
        - Image-based Time Series exposures
        - Single science exposure
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint_Image_Science(),
            Constraint_Single_Science(self.has_science),
            Constraint_TSO(),
        ])

        # Now check and continue initialization.
        super(Asn_Lv2ImageTSO, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2ImageTSO, self)._init_hook(item)
        self.data['asn_type'] = 'tso-image2'


@RegistryMarker.rule
class Asn_Lv2FGS(
        AsnMixin_Lv2Image,
        DMSLevel2bBase
):
    """Level2b FGS Association

    Characteristics:
        - Association type: ``image2``
        - Pipeline: ``calwebb_image2``
        - Image-based FGS science exposures
        - Single science exposure
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Single_Science(self.has_science),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value=(
                    'fgs_image'
                    '|fgs_focus'
                ),
            ),
            Constraint(
                [Constraint_WFSC()],
                reduce=Constraint.notany
            )
        ])

        super(Asn_Lv2FGS, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2Spec(
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b Science Spectral Association

    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - Spectral-based single target science exposures
        - Single science exposure
        - Non-TSO
        - Not part of a background dither observation
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint_Spectral_Science(
                exclude_exp_types=['nis_wfss', 'nrc_wfss', 'nrs_fixedslit', 'nrs_msaspec']
            ),
            Constraint(
                [
                    Constraint_Single_Science(self.has_science),
                    SimpleConstraint(
                        value='science',
                        test=lambda value, item: self.get_exposure_type(item) != value,
                        force_unique=False,
                    )
                ],
                reduce=Constraint.any
            ),
            Constraint(
                [
                    Constraint_TSO(),
                    DMSAttrConstraint(
                        name='patttype',
                        sources=['patttype'],
                        value=['2-point-nod|4-point-nod'],
                    )
                ],
                reduce=Constraint.notany
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv2Spec, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2SpecSpecial(
        AsnMixin_Lv2Special,
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b Auxiliary Science Spectral Association

    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - Spectral-based single target science exposures that are background exposures
        - Single science exposure
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint_Spectral_Science(),
            Constraint_Single_Science(self.has_science),
            Constraint_Special(),
        ])

        # Now check and continue initialization.
        super(Asn_Lv2SpecSpecial, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2SpecTSO(
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b Time Series Science Spectral Association

    Characteristics:
        - Association type: ``tso-spec2``
        - Pipeline: ``calwebb_tso-spec2``
        - Spectral-based single target time series exposures
        - Single science exposure
        - No other exposure can be part of the association
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint_Spectral_Science(
                exclude_exp_types=['nrs_msaspec', 'nrs_fixedslit']
            ),
            Constraint_Single_Science(self.has_science),
            Constraint_TSO(),
        ])

        # Now check and continue initialization.
        super(Asn_Lv2SpecTSO, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2SpecTSO, self)._init_hook(item)
        self.data['asn_type'] = 'tso-spec2'


@RegistryMarker.rule
class Asn_Lv2MIRLRSFixedSlitNod(
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b MIRI LRS Fixed Slit background nods Association

    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - MIRI LRS Fixed slit
        - Single science exposure
        - Include slit nods as backgrounds
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='mir_lrs-fixedslit'
            ),
            DMSAttrConstraint(
                name='patttype',
                sources=['patttype'],
                value=['along-slit-nod'],
            ),
            SimpleConstraint(
                value=True,
                test=lambda value, item: self.acid.type != 'background',
                force_unique=False
            ),
            Constraint(
                [
                    Constraint(
                        [
                            DMSAttrConstraint(
                                name='patt_num',
                                sources=['patt_num'],
                            ),
                            Constraint_Single_Science(
                                self.has_science,
                                reprocess_on_match=True,
                                work_over=ProcessList.EXISTING
                            )
                        ]
                    ),
                    Constraint(
                        [
                            DMSAttrConstraint(
                                name='is_current_patt_num',
                                sources=['patt_num'],
                                value=lambda: '((?!{}).)*'.format(self.constraints['patt_num'].value),
                            ),
                            SimpleConstraint(
                                name='force_match',
                                value=None,
                                sources=lambda item: False,
                                test=lambda constraint, obj: True,
                                force_unique=True,
                            )
                        ]
                    )
                ],
                reduce=Constraint.any
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv2MIRLRSFixedSlitNod, self).__init__(*args, **kwargs)

    def get_exposure_type(self, item, default='science'):
        """Modify exposure type depending on dither pointing index

        Behaves as the superclass method. However, if the constraint
        `is_current_patt_num` is True, mark the exposure type as
        `background`.
        """
        exp_type = super(Asn_Lv2MIRLRSFixedSlitNod, self).get_exposure_type(
            item, default
        )
        if exp_type == 'science' and self.constraints['is_current_patt_num'].matched:
            exp_type = 'background'

        return exp_type


@RegistryMarker.rule
class Asn_Lv2NRSLAMPSpectral(
        AsnMixin_Lv2Special,
        DMSLevel2bBase
):
    """Level2b NIRSpec spectral Lamp Calibrations Association

    Characteristics:
        - Association type: ``nrslamp-spec2``
        - Pipeline: ``calwebb_nrslamp-spec2``
        - Spectral-based calibration exposures
        - Single science exposure
    """

    def __init__(self, *args, **kwargs):

        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Single_Science(self.has_science),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='nrs_autoflat|nrs_autowave|nrs_lamp'
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['filter'],
                value='opaque'
            ),
            DMSAttrConstraint(
                name='opmode',
                sources=['opmode'],
            ),
            DMSAttrConstraint(
                name='lamp',
                sources=['lamp'],
            ),
            Constraint(
                [
                    DMSAttrConstraint(
                        sources=['grating'],
                        value='mirror',
                        force_unique=False,
                    ),
                    DMSAttrConstraint(
                        sources=['opmode'],
                        value='grating-only',
                        force_unique=False,
                    ),
                    DMSAttrConstraint(
                        sources=['lamp'],
                        value='nolamp',
                        force_unique=False,
                    ),
                ],
                reduce=Constraint.notany
            ),
        ])

        super(Asn_Lv2NRSLAMPSpectral, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2NRSLAMPSpectral, self)._init_hook(item)
        self.data['asn_type'] = 'nrslamp-spec2'


@RegistryMarker.rule
class Asn_Lv2WFSS(
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b WFSS/GRISM Association

    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - Mutli-object science exposures
        - Single science exposure
        - Require a source catalog from processing of the corresponding direct imagery.
    """

    def __init__(self, *args, **kwargs):

        self.constraints = Constraint([

            # Basic constraints
            Constraint_Base(),
            Constraint_Target(),

            # Allow WFSS exposures but account for the direct imaging.
            Constraint([

                # Constrain on the WFSS exposure
                Constraint([
                    DMSAttrConstraint(
                        name='exp_type',
                        sources=['exp_type'],
                        value='nis_wfss|nrc_wfss',
                    ),
                    Constraint_Mode(),
                    Constraint_Single_Science(self.has_science),
                ]),

                # Or select related imaging exposures.
                DMSAttrConstraint(
                    name='image_exp_type',
                    sources=['exp_type'],
                    value='nis_image|nrc_image',
                    force_reprocess=ProcessList.EXISTING,
                    only_on_match=True,
                ),
            ], reduce=Constraint.any)
        ])

        super(Asn_Lv2WFSS, self).__init__(*args, **kwargs)

    def add_catalog_members(self):
        """Add catalog and direct image member based on direct image members"""
        directs = self.members_by_type('direct_image')
        if not directs:
            raise AssociationNotValidError(
                '{} has no required direct image exposures'.format(
                    self.__class__.__name__
                )
            )

        sciences = self.members_by_type('science')
        if not sciences:
            raise AssociationNotValidError(
                '{} has no required science exposure'.format(
                    self.__class__.__name__
                )
            )
        science = sciences[0]

        # Get the exposure sequence for the science. Then, find
        # the direct image greater than but closest to this value.
        closest = directs[0]  # If the search fails, just use the first.
        try:
            expspcin = int(getattr_from_list(science.item, ['expspcin'], _EMPTY)[1])
        except KeyError:
            # If exposure sequence cannot be determined, just fall through.
            logger.debug('Science exposure %s has no EXPSPCIN defined.', science)
        else:
            min_diff = -1         # Initialize to an invalid value.
            for direct in directs:
                try:
                    direct_expspcin = int(getattr_from_list(
                        direct.item, ['expspcin'], _EMPTY
                    )[1])
                except KeyError:
                    # Try the next one.
                    logger.debug('Direct image %s has no EXPSPCIN defined.', direct)
                    continue
                diff = direct_expspcin - expspcin
                if diff > min_diff:
                    min_diff = diff
                    closest = direct

        # Note the selected direct image. Used in `Asn_Lv2WFSS._get_opt_element`
        self.direct_image = closest

        # Remove all direct images from the association.
        members = self.current_product['members']
        direct_idxs = [
            idx
            for idx, member in enumerate(members)
            if member['exptype'] == 'direct_image'
        ]
        deque((
            list.pop(members, idx)
            for idx in sorted(direct_idxs, reverse=True)
        ))

        # Add the Level3 catalog and direct image members
        lv3_direct_image_root = DMS_Level3_Base._dms_product_name(self)
        members.append(
            Member({
                'expname': lv3_direct_image_root + '_i2d.fits',
                'exptype': 'direct_image'
            })
        )
        members.append(
            Member({
                'expname': lv3_direct_image_root + '_cat.ecsv',
                'exptype': 'sourcecat'
            })
        )

    def finalize(self):
        """Finalize the association

        For WFSS, this involves taking all the direct image exposures,
        determine which one is first after last science exposure,
        and creating the catalog name from that image.
        """
        try:
            self.add_catalog_members()
        except AssociationNotValidError as err:
            logger.debug(
                '%s: %s',
                self.__class__.__name__, str(err)
            )
            return None

        return super(Asn_Lv2WFSS, self).finalize()

    def get_exposure_type(self, item, default='science'):
        """Modify exposure type depending on dither pointing index

        If an imaging exposure as been found, treat is as a direct image.
        """
        exp_type = super(Asn_Lv2WFSS, self).get_exposure_type(
            item, default
        )
        if exp_type == 'science' and item['exp_type'] in ['nis_image', 'nrc_image']:
            exp_type = 'direct_image'

        return exp_type

    def _get_opt_element(self):
        """Get string representation of the optical elements

        Returns
        -------
        opt_elem: str
            The Level3 Product name representation
            of the optical elements.
        Notes
        -----
        This is an override for the method in `DMSBaseMixin`.
        The optical element is retieved from the chosen direct image
        found in `self.direct_image`, determined in the `self.finalize`
        method.
        """
        item = self.direct_image.item
        opt_elems = []
        for keys in [['filter', 'band'], ['pupil', 'grating']]:
            opt_elem = getattr_from_list_nofail(
                item, keys, _EMPTY
            )[1]
            if opt_elem:
                opt_elems.append(opt_elem)
        opt_elems.sort(key=str.lower)
        full_opt_elem = '-'.join(opt_elems)
        if full_opt_elem == '':
            full_opt_elem = 'clear'

        return full_opt_elem


@RegistryMarker.rule
class Asn_Lv2NRSMSA(
        AsnMixin_Lv2Nod,
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b NIRSpec MSA Association

    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - Spectral-based NIRSpec MSA multi-object science exposures
        - Single science exposure
        - Handle slitlet nodding for background subtraction
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint(
                [
                    DMSAttrConstraint(
                        name='exp_type',
                        sources=['exp_type'],
                        value='nrs_msaspec'
                    ),
                    DMSAttrConstraint(
                        sources=['msametfl']
                    ),
                    DMSAttrConstraint(
                        name='expspcin',
                        sources=['expspcin'],
                    )
                ]
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv2NRSMSA, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2NRSFSS(
        AsnMixin_Lv2Nod,
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b NIRSpec Fixed-slit Association

    Notes
    -----
    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - Spectral-based NIRSpec fixed-slit single target science exposures
        - Single science exposure
        - Handle along-the-slit background nodding

    Association includes both the background and science exposures of the nodding.
    The identified science exposure is fixed by the nod, pattern, and exposure number
    to prevent other science exposures being included.
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='nrs_fixedslit'
            ),
            SimpleConstraint(
                value=True,
                test=lambda value, item: nrsfss_valid_detector(item),
                force_unique=False
            ),
            Constraint(
                [
                    SimpleConstraint(
                        value='science',
                        test=lambda value, item: self.get_exposure_type(item) != value,
                        force_unique=False
                    ),
                    Constraint(
                        [
                            DMSAttrConstraint(
                                name='expspcin',
                                sources=['expspcin'],
                            ),
                            DMSAttrConstraint(
                                name='nods',
                                sources=['numdthpt'],
                            ),
                            DMSAttrConstraint(
                                name='subpxpts',
                                sources=['subpxpns', 'subpxpts'],
                            ),
                            SimpleConstraint(
                                value='science',
                                test=lambda value, item: self.get_exposure_type(item) == value,
                                force_unique=False
                            )
                        ]
                    ),
                ],
                reduce=Constraint.any
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv2NRSFSS, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2NRSIFUNod(
        AsnMixin_Lv2Nod,
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b NIRSpec IFU Association

    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - Spectral-based NIRSpec IFU multi-object science exposures
        - Single science exposure
        - Handle 2 and 4 point background nodding
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='nrs_ifu'
            ),
            SimpleConstraint(
                value=True,
                test=lambda value, item: nrsifu_valid_detector(item),
                force_unique=False
            ),
            DMSAttrConstraint(
                name='expspcin',
                sources=['expspcin'],
            ),
            DMSAttrConstraint(
                name='patttype',
                sources=['patttype'],
                value=['2-point-nod|4-point-nod'],
                force_unique=True
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv2NRSIFUNod, self).__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2WFSC(
        DMSLevel2bBase
):
    """Level2b Wavefront Sensing & Control Association

    Characteristics:
        - Association type: ``wfs-image2``
        - Pipeline: ``calwebb_wfs-image2``
        - WFS and WFS&C observations
        - Single science exposure
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Image_Science(),
            Constraint_Single_Science(self.has_science),
            Constraint_ExtCal(),
            Constraint_WFSC(),
        ])

        # Now check and continue initialization.
        super(Asn_Lv2WFSC, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2WFSC, self)._init_hook(item)
        self.data['asn_type'] = 'wfs-image2'


@RegistryMarker.rule
class Asn_Force_Reprocess(DMSLevel2bBase):
    """Force all backgrounds to reprocess"""

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            SimpleConstraint(
                value='background',
                sources=self.get_exposure_type,
                force_unique=False,
            ),
            SimpleConstraint(
                name='force_fail',
                test=lambda x, y: False,
                value='anything but None',
                reprocess_on_fail=True,
                work_over=ProcessList.EXISTING,
                reprocess_rules=[]
            )
        ])

        super(Asn_Force_Reprocess, self).__init__(*args, **kwargs)
