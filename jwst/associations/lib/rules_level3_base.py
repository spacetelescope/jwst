"""Base classes which define the Level3 Associations"""
from collections import defaultdict
import logging
from os.path import basename
import re

from jwst.associations import (
    Association,
    AssociationRegistry,
    ProcessList,
    libpath
)
from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.utilities import (
    evaluate,
    is_iterable
)
from jwst.associations.exceptions import (
    AssociationNotAConstraint,
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
    DMSAttrConstraint,
    DMSBaseMixin,
    IMAGE2_SCIENCE_EXP_TYPES,
    IMAGE2_NONSCIENCE_EXP_TYPES,
    SPEC2_SCIENCE_EXP_TYPES,
    TSO_EXP_TYPES,
)
from jwst.associations.lib.format_template import FormatTemplate

__all__ = [
    'ASN_SCHEMA',
    'AsnMixin_Science',
    'AsnMixin_Spectrum',
    'Constraint_Base',
    'Constraint_IFU',
    'Constraint_Image',
    'Constraint_MSA',
    'Constraint_Optical_Path',
    'Constraint_Spectral',
    'Constraint_Target',
    'Constraint',
    'DMS_Level3_Base',
    'DMSAttrConstraint',
    'ProcessList',
    'SimpleConstraint',
    'Utility',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# The schema that these associations must adhere to.
ASN_SCHEMA = RegistryMarker.schema(libpath('asn_schema_jw_level3.json'))

# DMS file name templates
_LEVEL1B_REGEX = '(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'
_DMS_POOLNAME_REGEX = 'jw(\d{5})_(\d{8}[Tt]\d{6})_pool'

# Product name regex's
_REGEX_ACID_VALUE = '(o\d{3}|(c|a)\d{4})'

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
            'ok_candidate': {
                'validated': False,
                'check': self.ok_candidate
            }
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
        else:
            return NotImplemented

    def __ne__(self, other):
        """Compare inequality of two associations"""
        result = self.__eq__(other)
        if result is not NotImplemented:
            result = not result
        return result

    def dms_product_name(self):
        """Define product name.

        Returns
        -------
        product_name: str
            The product name
        """
        return self._dms_product_name(self)

    @staticmethod
    def _dms_product_name(association):
        """Define product name.

        Parameters
        ----------
        association: `Association`
            Association to get the name from.

        Returns
        -------
        product_name: str
            The product name
        """
        target = association._get_target()

        instrument = association._get_instrument()

        opt_elem = association._get_opt_element()

        exposure = association._get_exposure()
        if len(exposure):
            exposure = '-' + exposure

        subarray = association._get_subarray()
        if len(subarray):
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
        item: dict or None
            Item to use as a source. If not given, item-specific
            information will be left unchanged.

        member: dict or None
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
        product['name'] = self.dms_product_name()

    def make_member(self, item):
        """Create a member from the item

        Parameters
        ----------
        item: dict
            The item to create member from.

        Returns
        -------
        member: dict
            The member
        """
        try:
            exposerr = item['exposerr']
        except KeyError:
            exposerr = None

        # Get exposure type
        try:
            is_tso = self.constraints['is_tso'].matched
        except KeyError:
            is_tso = item['exp_type'] in TSO_EXP_TYPES

        exptype = self.get_exposure_type(item)

        # Determine expected member name
        expname = Utility.rename_to_level2(
                item['filename'], exp_type=item['exp_type'], is_tso=is_tso
            )

        member = {
            'expname': expname,
            'exptype': exptype,
            'exposerr': exposerr,
            'asn_candidate': item['asn_candidate']
        }
        return member

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
            logger.debug(
                'Member is already part of the association:'
                '\n\tassociation: {}'
                '\n]tmember: {}'.format(self, member)
            )
            return

        self.update_validity(member)
        members = self.current_product['members']
        members.append(member)
        if member['exposerr'] not in _EMPTY:
            logger.warn('Member {} has exposure error "{}"'.format(
                item['filename'],
                member['exposerr']
            ))
        self.from_items.append(item)

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
        items: [object[, ...]]
            A list of items to make members of the association.

        product_name: str or None
            The name of the product to add the items to.
            If the product does not already exist, it will be created.
            If None, the default DMS Level3 naming
            conventions will be attempted.

        with_exptype: bool
            If True, each item is expected to be a 2-tuple with
            the first element being the item to add as `expname`
            and the second items is the `exptype`

        kwargs: dict
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
            member = {
                'expname': item,
                'exptype': exptype
            }
            self.update_validity(member)
            members.append(member)
            self.from_items.append(item)
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
        member: dict
            Member being added causing check.
            Not used

        Returns
        -------
        is_valid: bool
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
    def rename_to_level2(level1b_name, exp_type=None, is_tso=False):
        """Rename a Level 1b Exposure to a Level2 name.

        The basic transform is changing the suffix `uncal` to
        `cal`, `calints`, or `rate`.

        Parameters
        ----------
        level1b_name: str
            The Level 1b exposure name.

        exp_type:
            JWST exposure type. If not specified,
            it will be presumed that the name
            should get a Level2b name

        is_tso: boolean
            Use 'calints' instead of 'cal' as
            the suffix.

        Returns
        -------
        str
            The Level 2b name
        """
        match = re.match(_LEVEL1B_REGEX, level1b_name)
        if match is None or match.group('type') != '_uncal':
            logger.warn((
                'Item FILENAME="{}" is not a Level 1b name. '
                'Cannot transform to Level 2b.'
            ).format(
                level1b_name
            ))
            return level1b_name

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
        value: str
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
        finalized_associations: [association[, ...]]
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
        'source_id': 's{:05d}'
    }
)


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
            ),
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
            ),
            DMSAttrConstraint(
                name='opt_elem2',
                sources=['pupil', 'grating'],
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


class Constraint_MSA(Constraint):
    """Constrain on NIRSpec MSA exposures that are spectral"""
    def __init__(self):
        super(Constraint_MSA, self).__init__(
            [
                DMSAttrConstraint(
                    name='exp_type',
                    sources=['exp_type'],
                    value=(
                        'nrs_autoflat'
                        '|nrs_autowave'
                        '|nrs_msaspec'
                    ),
                    force_unique=False
                ),
                DMSAttrConstraint(
                    name='is_msa',
                    sources=['msametfl'],
                )
            ],
            name='msa_spectral'
        )


class Constraint_Target(DMSAttrConstraint):
    """Select on target"""
    def __init__(self):
        super(Constraint_Target, self).__init__(
            name='target',
            sources=['targetid'],
        )


# -----------
# Base Mixins
# -----------
class AsnMixin_Science(DMS_Level3_Base):
    """Basic science constraints"""

    def __init__(self, *args, **kwargs):

        # Setup target acquisition inclusion
        constraint_acqs = Constraint(
            [
                DMSAttrConstraint(
                    name='acq_exp',
                    sources=['exp_type'],
                    value='|'.join(ACQ_EXP_TYPES),
                    force_unique=False
                ),
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
                    sources=['is_imprt', 'bkgdtarg'],
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
