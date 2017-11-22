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
from jwst.associations.association import (
    evaluate,
    is_iterable
)
from jwst.associations.exceptions import (
    AssociationNotAConstraint,
    AssociationNotValidError,
)
from jwst.associations.lib.acid import ACID
from jwst.associations.lib.counter import Counter
from jwst.associations.lib.dms_base import (DMSBaseMixin, _EMPTY)
from jwst.associations.lib.format_template import FormatTemplate

__all__ = [
    'AsnMixin_Base',
    'AsnMixin_CrossCandidate',
    'AsnMixin_Image',
    'AsnMixin_MIRI',
    'AsnMixin_NIRCAM',
    'AsnMixin_NIRISS',
    'AsnMixin_NIRSPEC',
    'AsnMixin_OpticalPath',
    'AsnMixin_Spectrum',
    'AsnMixin_Target',
    'ASN_SCHEMA',
    'DMS_Level3_Base',
    'ProcessList',
    'Utility',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# The schema that these associations must adhere to.
ASN_SCHEMA = libpath('asn_schema_jw_level3.json')

# Degraded status information
_DEGRADED_STATUS_OK = (
    'No known degraded exposures in association.'
)
_DEGRADED_STATUS_NOTOK = (
    'One or more members have an error associated with them.'
    '\nDetails can be found in the members.exposerr property.'
)

# DMS file name templates
_LEVEL1B_REGEX = '(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'
_DMS_POOLNAME_REGEX = 'jw(\d{5})_(\d{8}[Tt]\d{6})_pool'

# Product name regex's
_REGEX_ACID_VALUE = '(o\d{3}|(c|a)\d{4})'

# Key that uniquely identfies members.
KEY = 'expname'

# Exposures that are always TSO
TSO_EXP_TYPES = (
    'mir_lrs-slitless',
    'nis_soss',
    'nrc_tsimage',
    'nrc_tsgrism',
    'nrs_brightobj'
)


class DMS_Level3_Base(DMSBaseMixin, Association):
    """Basic class for DMS Level3 associations."""

    # Set the validation schema
    schema_file = ASN_SCHEMA

    # Attribute values that are indicate the
    # attribute is not specified.
    INVALID_VALUES = _EMPTY

    # Make sequences type-dependent
    _sequences = defaultdict(Counter)

    def __init__(self, *args, **kwargs):

        super(DMS_Level3_Base, self).__init__(*args, **kwargs)

        # Keep the set of members included in this association
        self.members = set()

        # Initialize validity checks
        self.validity.update({
            'has_science': {
                'validated': False,
                'check': lambda member: member['exptype'] == 'science'
            }
        })

        # Other presumptions on the association
        if 'degraded_status' not in self.data:
            self.data['degraded_status'] = _DEGRADED_STATUS_OK
        if 'program' not in self.data:
            self.data['program'] = 'noprogram'
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
            result = result and (self.members == other.members)
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

        try:
            exposure = association._get_exposure()
        except AssociationNotAConstraint:
            exposure = ''
        else:
            exposure = '-' + exposure

        product_name = 'jw{}-{}_{}_{}_{}'.format(
            association.data['program'],
            association.acid.id,
            target,
            instrument,
            opt_elem,
            exposure
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
        # Constraints
        self.data['constraints'] = '\n'.join(
            [cc for cc in self.constraints_to_text()]
        )

        # ID
        self.data['asn_id'] = self.acid.id

        # Target
        self.data['target'] = self._get_target()

        # Item-based information
        if item is not None:

            # Program
            if self.data['program'] != 'noprogram':
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

            # Degrade exposure
            if self.data['degraded_status'] == _DEGRADED_STATUS_OK:
                try:
                    exposerr = item['exposerr']
                except KeyError:
                    pass
                else:
                    if exposerr not in _EMPTY:
                        self.data['degraded_status'] = _DEGRADED_STATUS_NOTOK

        # Member-based info
        if member is not None:

            # Degraded exposure
            if self.data['degraded_status'] == _DEGRADED_STATUS_OK:
                try:
                    exposerr = member['exposerr']
                except KeyError:
                    pass
                else:
                    if exposerr not in _EMPTY:
                        self.data['degraded_status'] = _DEGRADED_STATUS_NOTOK

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

        try:
            is_tso = self.constraints['is_tso']['value'] == 't'
        except KeyError:
            is_tso = item['exp_type'] in TSO_EXP_TYPES

        member = {
            'expname': Utility.rename_to_level2b(
                item['filename'], is_tso=is_tso
            ),
            'exptype': self.get_exposure_type(item),
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

        # Add member to the short list
        self.members.add(member[KEY])

        # Update meta info
        self.update_asn(item=item, member=member)

    def _add_items(self, items, product_name=None, with_exptype=False):
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
    def rename_to_level2b(level1b_name, is_tso=False):
        """Rename a Level 1b Exposure to another level

        Parameters
        ----------
        level1b_name: str
            The Level 1b exposure name.

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

        suffix = 'cal'
        if is_tso:
            suffix = 'calints'
        level2b_name = ''.join([
            match.group('path'),
            '_',
            suffix,
            match.group('extension')
        ])
        return level2b_name

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
    @AssociationRegistry.callback('finalize')
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
        finalized = []
        lv3_asns = []
        for asn in associations:
            if isinstance(asn, DMS_Level3_Base):

                # Check validity
                if asn.is_valid:
                    lv3_asns.append(asn)

            else:
                finalized.append(asn)

        # Ensure sequencing is correct.
        Utility.resequence(lv3_asns)

        # Merge lists and return
        return finalized + lv3_asns


# ---------
# Utilities
# ---------
# Define default product name filling
format_product = FormatTemplate(
    key_formats={
        'source_id': 's{:05d}'
    }
)


# ---------------------------------------------
# Mixins to define the broad category of rules.
# ---------------------------------------------


class AsnMixin_Base(DMS_Level3_Base):
    """Restrict to Program and Instrument"""

    def __init__(self, *args, **kwargs):

        # I am defined by the following constraints
        self.add_constraints({
            'program': {
                'value': None,
                'inputs': ['program'],
            },
            'instrument': {
                'value': None,
                'inputs': ['instrume']
            },
        })

        super(AsnMixin_Base, self).__init__(*args, **kwargs)


class AsnMixin_OpticalPath(DMS_Level3_Base):
    """Ensure unique optical path"""

    def __init__(self, *args, **kwargs):
        # I am defined by the following constraints
        self.add_constraints({
            'opt_elem': {
                'value': None,
                'inputs': ['filter']
            },
            'opt_elem2': {
                'value': None,
                'inputs': ['pupil', 'grating'],
                'required': False,
            },
        })

        super(AsnMixin_OpticalPath, self).__init__(*args, **kwargs)


class AsnMixin_Target(DMS_Level3_Base):
    """Constrain by Target"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'target': {
                'value': None,
                'inputs': ['targetid'],
            },
        })

        # Check and continue initialization.
        super(AsnMixin_Target, self).__init__(*args, **kwargs)


class AsnMixin_NotTSO(DMS_Level3_Base):
    """Ensure exposure is not a TSO"""
    def __init__(self, *args, **kwargs):
        self.add_constraints({
            'is_not_tso': {
                'value': '[^t]',
                'inputs': ['tsovisit'],
                'required': False
            }
        })

        super(AsnMixin_NotTSO, self).__init__(*args, **kwargs)


class AsnMixin_MIRI(DMS_Level3_Base):
    """All things that belong to MIRI"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'instrument': {
                'value': 'miri',
                'inputs': ['instrume']
            }
        })

        # Check and continue initialization.
        super(AsnMixin_MIRI, self).__init__(*args, **kwargs)


class AsnMixin_NIRSPEC(DMS_Level3_Base):
    """All things that belong to NIRSPEC"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'instrument': {
                'value': 'nirspec',
                'inputs': ['instrume']
            }
        })

        # Check and continue initialization.
        super(AsnMixin_NIRSPEC, self).__init__(*args, **kwargs)


class AsnMixin_NIRISS(DMS_Level3_Base):
    """All things that belong to NIRISS"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'instrument': {
                'value': 'niriss',
                'inputs': ['instrume']
            },
        })

        # Check and continue initialization.
        super(AsnMixin_NIRISS, self).__init__(*args, **kwargs)


class AsnMixin_NIRCAM(DMS_Level3_Base):
    """All things that belong to NIRCAM"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'instrument': {
                'value': 'nircam',
                'inputs': ['instrume']
            },
        })

        # Check and continue initialization.
        super(AsnMixin_NIRCAM, self).__init__(*args, **kwargs)


class AsnMixin_Image(AsnMixin_NotTSO):
    """All things that are in imaging mode"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': 'nrc_image|mir_image|nis_image|fgs_image',
                'inputs': ['exp_type'],
                'force_unique': True,
            }
        })

        super(AsnMixin_Image, self).__init__(*args, **kwargs)



class AsnMixin_Spectrum(AsnMixin_NotTSO):
    """All things that are spectrum"""

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'spec3'
        super(AsnMixin_Spectrum, self)._init_hook(item)


class AsnMixin_CrossCandidate(DMS_Level3_Base):
    """Basic constraints for Cross-Candidate associations"""

    @classmethod
    def validate(cls, asn):
        super(AsnMixin_CrossCandidate, cls).validate(asn)

        if isinstance(asn, AsnMixin_CrossCandidate):
            try:
                candidates = set(
                    member['asn_candidate_id']
                    for product in asn.data['products']
                    for member in product['members']
                )
            except (AttributeError, KeyError) as err:
                raise AssociationNotValidError('Validation failed')
            if not len(candidates) > 1:
                raise AssociationNotValidError(
                    'Validation failed: No candidates found.'
                )
        return True
