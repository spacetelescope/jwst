"""Base classes which define the Level3 Associations"""
from collections import defaultdict
import logging
from os.path import (
    basename,
    split,
    splitext
    )
import re

from jwst.associations import (
    Association,
    ProcessList,
    libpath
)
from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.utilities import (
    evaluate,
    is_iterable
)
from jwst.associations.exceptions import (
    AssociationNotValidError,
)
from jwst.associations.lib.acid import ACID
from jwst.associations.lib.constraint import (
    Constraint,
    SimpleConstraint,
)
from jwst.associations.lib.counter import Counter
from jwst.associations.lib.dms_base import (
    _EMPTY,
    ACQ_EXP_TYPES,
    Constraint_TargetAcq,
    CORON_EXP_TYPES,
    DMSAttrConstraint,
    DMSBaseMixin,
    IMAGE2_SCIENCE_EXP_TYPES,
    IMAGE2_NONSCIENCE_EXP_TYPES,
    SPEC2_SCIENCE_EXP_TYPES,
)
from jwst.associations.lib.format_template import FormatTemplate
from jwst.associations.lib.member import Member
from jwst.associations.lib.product_utils import prune_duplicate_associations

__all__ = [
    'ASN_SCHEMA',
    'AsnMixin_AuxData',
    'AsnMixin_Science',
    'AsnMixin_Spectrum',
    'Constraint',
    'Constraint_Base',
    'Constraint_IFU',
    'Constraint_Image',
    'Constraint_MSA',
    'Constraint_Obsnum',
    'Constraint_Optical_Path',
    'Constraint_Spectral',
    'Constraint_Target',
    'DMS_Level3_Base',
    'DMSAttrConstraint',
    'ProcessList',
    'SimpleConstraint',
    'Utility',
]
from jwst.lib.suffix import remove_suffix
# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# The schema that these associations must adhere to.
ASN_SCHEMA = RegistryMarker.schema(libpath('asn_schema_jw_level3.json'))

# DMS file name templates
_LEVEL1B_REGEX = r'(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'
_DMS_POOLNAME_REGEX = r'jw(\d{5})_(\d{8}[Tt]\d{6})_pool'

# Product name regex's
_REGEX_ACID_VALUE = r'(o\d{3}|(c|a)\d{4})'

# Exposures that should have received Level2b processing
LEVEL2B_EXPTYPES = []
LEVEL2B_EXPTYPES.extend(IMAGE2_SCIENCE_EXP_TYPES)
LEVEL2B_EXPTYPES.extend(IMAGE2_NONSCIENCE_EXP_TYPES)
LEVEL2B_EXPTYPES.extend(SPEC2_SCIENCE_EXP_TYPES)

# Association Candidates that should never make Level3 associations
INVALID_AC_TYPES = ['background']


class DMS_Level3_Base(DMSBaseMixin, Association):
    """Basic class for DMS Level3 associations."""

    # Set the validation schema
    schema_file = ASN_SCHEMA.schema

    # Attribute values that are indicate the
    # attribute is not specified.
    INVALID_VALUES = _EMPTY

    # Make sequences type-dependent
    _sequences = defaultdict(Counter)

    def __init__(self, *args, **kwargs):

        super(DMS_Level3_Base, self).__init__(*args, **kwargs)

        # Initialize validity checks
        self.validity.update({
            'has_science': {
                'validated': False,
                'check': lambda member: member['exptype'] == 'science'
            },
        })

        # Other presumptions on the association
        if 'constraints' not in self.data:
            self.data['constraints'] = 'No constraints'
        if 'asn_type' not in self.data:
            self.data['asn_type'] = 'user_built'
        if 'asn_id' not in self.data:
            self.data['asn_id'] = 'a3001'
        if 'target' not in self.data:
            self.data['target'] = 'none'
        if 'asn_pool' not in self.data:
            self.data['asn_pool'] = 'none'

    @property
    def current_product(self):
        return self.data['products'][-1]

    def __eq__(self, other):
        """Compare equality of two assocaitions"""
        if isinstance(other, DMS_Level3_Base):
            result = self.data['asn_type'] == other.data['asn_type']
            result = result and (self.member_ids == other.member_ids)
            return result

        return NotImplemented

    def __ne__(self, other):
        """Compare inequality of two associations"""
        result = self.__eq__(other)
        if result is not NotImplemented:
            result = not result
        return result

    @property
    def dms_product_name(self):
        """Define product name.

        Returns
        -------
        product_name : str
            The product name
        """
        return self._dms_product_name(self)

    @staticmethod
    def _dms_product_name(association):
        """Define product name.

        Parameters
        ----------
        association : `Association`
            Association to get the name from.

        Returns
        -------
        product_name : str
            The product name
        """
        target = association._get_target()

        instrument = association._get_instrument()

        opt_elem = association._get_opt_element()

        exposure = association._get_exposure()
        if exposure:
            exposure = '-' + exposure

        subarray = association._get_subarray()
        if subarray:
            subarray = '-' + subarray

        product_name = (
            'jw{program}-{acid}'
            '_{target}'
            '_{instrument}'
            '_{opt_elem}{subarray}'
        )
        product_name = product_name.format(
            program=association.data['program'],
            acid=association.acid.id,
            target=target,
            instrument=instrument,
            opt_elem=opt_elem,
            subarray=subarray,
            exposure=exposure
        )

        return product_name.lower()

    def update_asn(self, item=None, member=None):
        """Update association meta information

        Parameters
        ----------
        item : dict or None
            Item to use as a source. If not given, item-specific
            information will be left unchanged.

        member : Member or None
            An association member to use as source.
            If not given, member-specific information will be update
            from current association/product membership.

        Notes
        -----
        If both `item` and `member` are given,
        information in `member` will take precedence.
        """
        super(DMS_Level3_Base, self).update_asn(item=item, member=member)

        # Constraints
        self.data['constraints'] = str(self.constraints)

        # ID
        self.data['asn_id'] = self.acid.id

        # Target
        self.data['target'] = self._get_target()

        # Item-based information
        if item is not None:

            # Program
            if self.data['program'] == 'noprogram':
                self.data['program'] = '{:0>5s}'.format(item['program'])

            # Pool
            if self.data['asn_pool'] == 'none':
                self.data['asn_pool'] = basename(
                    item.meta['pool_file']
                ).split('.')[0]
                parsed_name = re.search(
                    _DMS_POOLNAME_REGEX, self.data['asn_pool']
                )
                if parsed_name is not None:
                    pool_meta = {
                        'program_id': parsed_name.group(1),
                        'version': parsed_name.group(2)
                    }
                    self.meta['pool_meta'] = pool_meta

        # Product-based updates
        product = self.current_product
        product['name'] = self.dms_product_name

    def make_member(self, item):
        """Create a member from the item

        Parameters
        ----------
        item : dict
            The item to create member from.

        Returns
        -------
        member : Member
            The member
        """
        try:
            exposerr = item['exposerr']
        except KeyError:
            exposerr = None

        # Get exposure type
        exptype = self.get_exposure_type(item)

        # Determine expected member name
        expname = Utility.rename_to_level2(
            item['filename'], exp_type=item['exp_type'],
            is_tso=self.is_item_tso(item, other_exp_types=CORON_EXP_TYPES),
            member_exptype=exptype
        )

        member = Member(
            {
                'expname': expname,
                'exptype': exptype,
                'exposerr': exposerr,
                'asn_candidate': item['asn_candidate']
            },
            item=item
        )
        return member

    def make_fixedslit_bkg(self):
        """Add a background to a MIR_lrs-fixedslit observation"""

        # check to see if these are nodded backgrounds, if they are setup
        # the background members, otherwise return the original association
        # to test for the string 'nod' we need to copy and pop the value out of the set
        if 'nod' not in self.constraints['patttype_spectarg'].found_values.copy().pop():
            results = []
            results.append(self)
            return results

        for product in self['products']:
            members = product['members']
            # Split out the science exposures
            science_exps = [
                member
                for member in members
                if member['exptype'] == 'science'
            ]
            # if there is only one science observation it cannot be the background
            # return with original association.
            if len(science_exps) < 2:
                return results

            # Create new members for each science exposure in the association,
            # using the the base name + _x1d as background.
            results = []
            # Loop over all science exposures in the association
            for science_exp in science_exps:
                sci_name = science_exp['expname']
                science_exp['expname'] = sci_name
                # Construct the name for the background file
                bkg_name = remove_suffix(
                    splitext(split(science_exp['expname'])[1])[0])[0]
                bkg_name = bkg_name+'_x1d.fits'
                now_background = Member(science_exp)
                now_background['expname'] = bkg_name
                now_background['exptype'] = 'background'
                # Add the background file to the association table
                members.append(now_background)

            if self.is_valid:
                results.append(self)

            return results


    def _init_hook(self, item):
        """Post-check and pre-add initialization"""
        super(DMS_Level3_Base, self)._init_hook(item)

        # Set which sequence counter should be used.
        self._sequence = self._sequences[self.data['asn_type']]

        # Create the product.
        self.new_product()

        # Update meta data
        self.update_asn(item=item)

    def _add(self, item):
        """Add item to this association."""
        member = self.make_member(item)
        if self.is_member(member):
            # logger.debug(
            #     'Member is already part of the association:'
            #     '\n\tassociation: {}'
            #     '\n]tmember: {}'.format(self, member)
            # )
            return

        self.update_validity(member)
        members = self.current_product['members']
        members.append(member)
        if member['exposerr'] not in _EMPTY:
            logger.warning('Member {} has exposure error "{}"'.format(
                item['filename'],
                member['exposerr']
            ))

        # Update meta info
        self.update_asn(item=item, member=member)

    def _add_items(self,
                   items,
                   product_name=None,
                   with_exptype=False,
                   **kwargs):
        """ Force adding items to the association

        Parameters
        ----------
        items : [object[, ...]]
            A list of items to make members of the association.

        product_name : str or None
            The name of the product to add the items to.
            If the product does not already exist, it will be created.
            If None, the default DMS Level3 naming
            conventions will be attempted.

        with_exptype : bool
            If True, each item is expected to be a 2-tuple with
            the first element being the item to add as `expname`
            and the second items is the `exptype`

        kwargs : dict
            Allows other keyword arguments used by other subclasses.

        Notes
        -----
        This is a low-level shortcut into adding members, such as file names,
        to an association. All defined shortcuts and other initializations are
        by-passed, resulting in a potentially unusable association.
        """
        if product_name is None:
            raise AssociationNotValidError(
                'Product name needs to be specified'
            )

        self.new_product(product_name)
        members = self.current_product['members']
        for item in items:
            exptype = 'science'
            if with_exptype:
                item, exptype = item
            member = Member(
                {
                    'expname': item,
                    'exptype': exptype
                },
                item=item
            )
            self.update_validity(member)
            members.append(member)
        self.sequence = next(self._sequence)

    def __repr__(self):
        try:
            file_name, json_repr = self.ioregistry['json'].dump(self)
        except:
            return str(self.__class__)
        return json_repr

    def __str__(self):
        result_list = []
        result_list.append(
            '{} with {} products'.format(
                self.asn_name,
                len(self.data['products'])
            )
        )
        result_list.append(
            'Rule={}'.format(self.data['asn_rule'])
        )
        result_list.append(self.data['constraints'])
        result_list.append('Products:')
        for product in self.data['products']:
            result_list.append(
                '\t{} with {} members'.format(
                    product['name'],
                    len(product['members'])
                )
            )
        result = '\n'.join(result_list)
        return result

    def ok_candidate(self, member=None):
        """Validation test for acceptable candidates

        Parameters
        ----------
        member : Member
            Member being added causing check.
            Not used

        Returns
        -------
        is_valid : bool
        """
        return self.acid.type.lower() not in INVALID_AC_TYPES


@RegistryMarker.utility
class Utility():
    """Utility functions that understand DMS Level 3 associations"""

    @staticmethod
    def resequence(associations):
        """Resequence the numbering for the Level3 association types"""
        counters = defaultdict(lambda: defaultdict(Counter))
        for asn in associations:
            asn.sequence = next(
                counters[asn.data['asn_id']][asn.data['asn_type']]
            )

    @staticmethod
    def rename_to_level2(level1b_name, exp_type=None, is_tso=False, member_exptype='science'):
        """Rename a Level 1b Exposure to a Level2 name.

        The basic transform is changing the suffix `uncal` to
        `cal`, `calints`, or `rate`.

        Parameters
        ----------
        level1b_name : str
            The Level 1b exposure name.

        exp_type:
            JWST exposure type. If not specified,
            it will be presumed that the name
            should get a Level2b name

        is_tso : boolean
            Use 'calints' instead of 'cal' as
            the suffix.

        member_exptype: str
            The assocition member exposure type, such as "science".

        Returns
        -------
        str
            The Level 2b name
        """
        match = re.match(_LEVEL1B_REGEX, level1b_name)
        if match is None or match.group('type') != '_uncal':
            logger.warning((
                'Item FILENAME="{}" is not a Level 1b name. '
                'Cannot transform to Level 2b.'
            ).format(
                level1b_name
            ))
            return level1b_name

        if member_exptype == 'background':
            suffix = 'x1d'
        else:
            if exp_type in LEVEL2B_EXPTYPES:
                suffix = 'cal'
            else:
                suffix = 'rate'

        if is_tso:
            suffix += 'ints'

        level2_name = ''.join([
            match.group('path'),
            '_',
            suffix,
            match.group('extension')
        ])
        return level2_name

    @staticmethod
    def get_candidate_list(value):
        """Parse the candidate list from a item value

        Parameters
        ----------
        value : str
            The value from the item to parse. Usually
            item['ASN_CANDIDATE']

        Returns
        -------
        [ACID, ...]
            The list of parsed candidates.
        """
        result = []
        evaled = evaluate(value)
        if is_iterable(evaled):
            result = [
                ACID(v)
                for v in evaled
            ]
        return result

    @staticmethod
    @RegistryMarker.callback('finalize')
    def finalize(associations):
        """Check validity and duplications in an association list

        Parameters
        ----------
        associations:[association[, ...]]
            List of associations

        Returns
        -------
        finalized_associations : [association[, ...]]
            The validated list of associations
        """
        finalized_asns = []
        lv3_asns = []
        for asn in associations:
            if isinstance(asn, DMS_Level3_Base):
                finalized = asn.finalize()
                if finalized is not None:
                    lv3_asns.extend(finalized)
            else:
                finalized_asns.append(asn)

        lv3_asns = prune_duplicate_associations(lv3_asns)

        # Ensure sequencing is correct.
        Utility.resequence(lv3_asns)

        # Merge lists and return
        return finalized_asns + lv3_asns


# ---------
# Utilities
# ---------
# Define default product name filling
format_product = FormatTemplate(
    key_formats={
        'source_id': ['s{:05d}', 's{:s}'],
        'expspcin': ['{:0>2s}']

    }
)


def dms_product_name_sources(asn):
    """Produce source-based product names

    Parameters
    ---------
    asn : Association
        The association for which the product
        name is to be created.

    Returns
    -------
    product_name : str
        The product name
    """
    instrument = asn._get_instrument()

    opt_elem = asn._get_opt_element()

    subarray = asn._get_subarray()
    if subarray:
        subarray = '-' + subarray

    product_name_format = (
        'jw{program}-{acid}'
        '_{source_id}'
        '_{instrument}'
        '_{opt_elem}{subarray}'
    )
    product_name = format_product(
        product_name_format,
        program=asn.data['program'],
        acid=asn.acid.id,
        instrument=instrument,
        opt_elem=opt_elem,
        subarray=subarray,
    )

    return product_name.lower()


# -----------------
# Basic constraints
# -----------------
class Constraint_Base(Constraint):
    """Select on program and instrument"""
    def __init__(self):
        super(Constraint_Base, self).__init__(
            [
                DMSAttrConstraint(
                    name='program',
                    sources=['program'],
                ),
                DMSAttrConstraint(
                    name='instrument',
                    sources=['instrume'],
                ),
            ],
            name='base'
        )


class Constraint_IFU(DMSAttrConstraint):
    """Constrain on IFU exposures"""
    def __init__(self):
        super(Constraint_IFU, self).__init__(
            name='exp_type',
            sources=['exp_type'],
            value=(
                'mir_mrs'
                '|mir_flatmrs'
                '|nrs_autowave'
                '|nrs_ifu'
            ),
            force_unique=False
        )


class Constraint_Image(DMSAttrConstraint):
    """Select on exposure type"""
    def __init__(self):
        super(Constraint_Image, self).__init__(
            name='exp_type',
            sources=['exp_type'],
            value=(
                'nrc_image'
                '|mir_image'
                '|nis_image'
                '|fgs_image'
                '|nrs_mimf'
            ),
        )


class Constraint_MSA(Constraint):
    """Constrain on NIRSpec MSA exposures that are spectral"""
    def __init__(self):
        super(Constraint_MSA, self).__init__(
            [
                DMSAttrConstraint(
                    name='exp_type',
                    sources=['exp_type'],
                    value=('nrs_msaspec'),
                    force_unique=False
                ),
                DMSAttrConstraint(
                    name='is_msa',
                    sources=['msametfl'],
                    force_unique=False,
                )
            ],
            name='msa_spectral'
        )


class Constraint_Obsnum(DMSAttrConstraint):
    """Select on OBSNUM"""
    def __init__(self):
        super(Constraint_Obsnum, self).__init__(
            name='obs_num',
            sources=['obs_num'],
            force_unique=False,
            required=False,
        )


class Constraint_Optical_Path(Constraint):
    """Select on optical path"""
    def __init__(self):
        super(Constraint_Optical_Path, self).__init__([
            DMSAttrConstraint(
                name='opt_elem',
                sources=['filter'],
                required=False,
            ),
            DMSAttrConstraint(
                name='opt_elem2',
                sources=['pupil', 'grating'],
                required=False,
            ),
            DMSAttrConstraint(
                name='opt_elem3',
                sources=['fxd_slit'],
                required=False,
            ),
            DMSAttrConstraint(
                name='subarray',
                sources=['subarray']
            )
        ])


class Constraint_Spectral(DMSAttrConstraint):
    """Constrain on spectral exposure types"""
    def __init__(self):
        super(Constraint_Spectral, self).__init__(
            name='exp_type',
            sources=['exp_type'],
            value=(
                'mir_lrs-fixedslit'
                '|nrc_grism'
                '|nrc_wfss'
                '|nrs_autoflat'
                '|nrs_autowave'
                '|nrs_fixedslit'
            ),
            force_unique=False
        )


class Constraint_Target(DMSAttrConstraint):
    """Select on target

    Parameters
    ----------
    association: Association
        If specified, use the `get_exposure_type` method
        to as part of the target selection.
    """
    def __init__(self, association=None):
        if association is None:
            super(Constraint_Target, self).__init__(
                name='target',
                sources=['targetid'],
            )
        else:
            super(Constraint_Target, self).__init__(
                name='target',
                sources=['targetid'],
                onlyif=lambda item: association.get_exposure_type(item) != 'background',
                force_reprocess=ProcessList.EXISTING,
                only_on_match=True,
            )


# -----------
# Base Mixins
# -----------
class AsnMixin_AuxData:
    """Process special and non-science exposures as science.
    """
    def get_exposure_type(self, item, default='science'):
        """Override to force exposure type to always be science
        Parameters
        ----------
        item : dict
            The pool entry for which the exposure type is determined
        default : str or None
            The default exposure type.
            If None, routine will raise LookupError
        Returns
        -------
        exposure_type : 'science'
            Returns as science for most Exposures
        exposure_type : 'target_acquisition'
            Returns target_acquisition for mir_tacq
        """
        NEVER_CHANGE = ['target_acquisition']
        exp_type = super().get_exposure_type(item, default=default)
        if exp_type in NEVER_CHANGE:
            return exp_type
        return 'science'


class AsnMixin_Science(DMS_Level3_Base):
    """Basic science constraints"""

    def __init__(self, *args, **kwargs):

        # Setup target acquisition inclusion
        constraint_acqs = Constraint(
            [
                Constraint_TargetAcq(),
                DMSAttrConstraint(
                    name='acq_obsnum',
                    sources=['obs_num'],
                    value=lambda: '('
                    + '|'.join(self.constraints['obs_num'].found_values)
                    + ')',
                    force_unique=False,
                )
            ],
            name='acq_constraint',
            work_over=ProcessList.EXISTING
        )

        # Put all constraints together.
        self.constraints = Constraint(
            [
                Constraint_Base(),
                DMSAttrConstraint(
                    sources=['is_imprt'],
                    force_undefined=True
                ),
                Constraint(
                    [
                        Constraint(
                            [
                                self.constraints,
                                Constraint_Obsnum()
                            ],
                            name='rule'
                        ),
                        constraint_acqs
                    ],
                    name='acq_check',
                    reduce=Constraint.any
                ),
            ],
            name='dmsbase_top'
        )

        super(AsnMixin_Science, self).__init__(*args, **kwargs)


class AsnMixin_Spectrum(AsnMixin_Science):
    """All things that are spectrum"""

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'spec3'
        super(AsnMixin_Spectrum, self)._init_hook(item)
