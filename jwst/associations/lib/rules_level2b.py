"""Association Definitions: DMS Level2b product associations
"""
import logging

from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.constraint import (Constraint, SimpleConstraint)
from jwst.associations.lib.dms_base import (
    Constraint_TSO,
    format_list
)
from jwst.associations.lib.process_list import ProcessList
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
    'Asn_Lv2NRSLAMPImage',
    'Asn_Lv2NRSLAMPSpectral',
    'Asn_Lv2NRSMSA',
    'Asn_Lv2Spec',
    'Asn_Lv2SpecSpecial',
    'Asn_Lv2SpecTSO',
    'Asn_Lv2WFSS_NIS',
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
            Constraint_Mode(),
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
                exclude_exp_types=['nrs_msaspec', 'nrs_fixedslit', 'nis_wfss']
            ),
            Constraint(
                [
                    Constraint_Single_Science(self.has_science),
                    SimpleConstraint(
                        value='science',
                        test=lambda value, item: self.get_exposure_type(item) != value,
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
                        value=['2-point-nod|4-point-nod|along-slit-nod'],
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
class Asn_Lv2NRSLAMPImage(
        AsnMixin_Lv2Special,
        DMSLevel2bBase
):
    """Level2b NIRSpec Image Lamp calibrations Association

    Characteristics:
        - Association type: ``nrslamp-image2``
        - Pipeline: ``calwebb_nrslamp-image2``
        - Image-based NRS calibrations
        - Single science exposure
    """

    def __init__(self, *args, **kwargs):

        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Single_Science(self.has_science),
            DMSAttrConstraint(
                name='opt_elem2',
                sources=['grating'],
                value='mirror'
            ),
            DMSAttrConstraint(
                name='instrument',
                sources=['instrume'],
                value='nirspec'
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['filter'],
                value='opaque'
            ),
        ])

        super(Asn_Lv2NRSLAMPImage, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2NRSLAMPImage, self)._init_hook(item)
        self.data['asn_type'] = 'nrslamp-image2'


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
            Constraint(
                [
                    DMSAttrConstraint(
                        name='opt_elem2',
                        sources=['grating'],
                        value='mirror'
                    )
                ],
                reduce=Constraint.notany
            ),
            DMSAttrConstraint(
                name='instrument',
                sources=['instrume'],
                value='nirspec'
            ),
            DMSAttrConstraint(
                name='opt_elem',
                sources=['filter'],
                value='opaque'
            ),
        ])

        super(Asn_Lv2NRSLAMPSpectral, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2NRSLAMPSpectral, self)._init_hook(item)
        self.data['asn_type'] = 'nrslamp-spec2'


@RegistryMarker.rule
class Asn_Lv2WFSS_NIS(
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b NIRISS WFSS/GRISM Association

    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - Spectral-based NIRISS mutli-object science exposures
        - Single science exposure
        - Require a source catalog from processing of the corresponding direct imagery.
    """

    def __init__(self, *args, **kwargs):

        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint_Target(),
            Constraint_Single_Science(self.has_science),
            DMSAttrConstraint(
                name='exp_type',
                sources=['exp_type'],
                value='nis_wfss',
            )
        ])

        super(Asn_Lv2WFSS_NIS, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""
        super(Asn_Lv2WFSS_NIS, self)._init_hook(item)

        # Get the Level3 product name of this association.
        # Except for the grism component, it should be what
        # the Level3 direct image name is.
        lv3_direct_image_catalog = DMS_Level3_Base._dms_product_name(self) + '_cat.ecsv'

        # Insert the needed catalog member
        member = {
            'expname': lv3_direct_image_catalog,
            'exptype': 'sourcecat'
        }
        members = self.current_product['members']
        members.append(member)

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
        The second optical element, the grism, would never be part
        of the direct image Level3 name.
        """
        opt_elem = ''
        try:
            value = format_list(self.constraints['opt_elem2'].found_values)
        except KeyError:
            pass
        else:
            if value not in _EMPTY and value != 'clear':
                opt_elem = value
        if opt_elem == '':
            opt_elem = 'clear'
        return opt_elem


@RegistryMarker.rule
class Asn_Lv2NRSMSA(
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

    def finalize(self):
        """Finalize assocation

        For NRS MSA, finalization means creating new associations for
        background nods.

        Returns
        -------
        associations: [association[, ...]] or None
            List of fully-qualified associations that this association
            represents.
            `None` if a complete association cannot be produced.

        """
        if self.is_valid:
            return self.make_nod_asns()
        else:
            return None


@RegistryMarker.rule
class Asn_Lv2NRSFSS(
        AsnMixin_Lv2Spectral,
        DMSLevel2bBase
):
    """Level2b NIRSpec Fixed-slit Association

    Characteristics:
        - Association type: ``spec2``
        - Pipeline: ``calwebb_spec2``
        - Spectral-based NIRSpec fixed-slit single target science exposures
        - Single science exposure
        - Handle along-the-slit background nodding
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Mode(),
            Constraint(
                [
                    Constraint(
                        [
                            DMSAttrConstraint(
                                name='exp_type',
                                sources=['exp_type'],
                                value='nrs_fixedslit'
                            ),
                            SimpleConstraint(
                                value='science',
                                test=lambda value, item: self.get_exposure_type(item) != value,
                                force_unique=False
                            )
                        ]
                    ),
                    Constraint(
                        [
                            DMSAttrConstraint(
                                name='exp_type',
                                sources=['exp_type'],
                                value='nrs_fixedslit'
                            ),
                            DMSAttrConstraint(
                                name='expspcin',
                                sources=['expspcin'],
                            ),
                            DMSAttrConstraint(
                                name='nods',
                                sources=['numdthpt'],
                            ),
                            DMSAttrConstraint(
                                name='subpxpns',
                                sources=['subpxpns'],
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

    def finalize(self):
        """Finalize assocation

        For NRS Fixed-slit, finalization means creating new associations for
        background nods.

        Returns
        -------
        associations: [association[, ...]] or None
            List of fully-qualified associations that this association
            represents.
            `None` if a complete association cannot be produced.

        """
        return self.make_nod_asns()


@RegistryMarker.rule
class Asn_Lv2NRSIFUNod(
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
            Constraint(
                [
                    DMSAttrConstraint(
                        name='exp_type',
                        sources=['exp_type'],
                        value='nrs_ifu'
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
                ]
            ),
        ])

        # Now check and continue initialization.
        super(Asn_Lv2NRSIFUNod, self).__init__(*args, **kwargs)

    def finalize(self):
        """Finalize assocation

        Finalization means creating new associations for
        background nods.

        Returns
        -------
        associations: [association[, ...]] or None
            List of fully-qualified associations that this association
            represents.
            `None` if a complete association cannot be produced.

        """
        return self.make_nod_asns()


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
            Constraint_Single_Science(self.has_science),
            DMSAttrConstraint(
                name='wfsc',
                sources=['visitype'],
                value='.+wfsc.+',
                force_unique=True
            )
        ])

        # Now check and continue initialization.
        super(Asn_Lv2WFSC, self).__init__(*args, **kwargs)

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(Asn_Lv2WFSC, self)._init_hook(item)
        self.data['asn_type'] = 'wfs-image2'
