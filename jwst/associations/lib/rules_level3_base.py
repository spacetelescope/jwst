"""Base classes which define the Level3 Associations"""
from collections import defaultdict
import logging
from os.path import basename
import re

from jwst.associations import (
    Association,
    AssociationRegistry,
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
from jwst.associations.lib.dms_base import DMSBaseMixin
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
    'Utility',
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Non-specified values found in DMS Association Pools
_EMPTY = (None, 'NULL', 'Null', 'null', '--')

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
                'check': lambda entry: entry['exptype'] == 'science'
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
        """Define product name."""
        target = self._get_target()

        instrument = self._get_instrument()

        opt_elem = self._get_opt_element()

        try:
            exposure = self._get_exposure()
        except AssociationNotAConstraint:
            exposure = ''
        else:
            exposure = '-' + exposure

        product_name = 'jw{}-{}_{}_{}_{}'.format(
            self.data['program'],
            self.acid.id,
            target,
            instrument,
            opt_elem,
            exposure
        )

        return product_name.lower()

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""
        super(DMS_Level3_Base, self)._init_hook(member)

        # Set which sequence counter should be used.
        self._sequence = self._sequences[self.data['asn_type']]

        self.data['target'] = member['targetid']
        self.data['program'] = str(member['program'])
        self.data['asn_pool'] = basename(
            member.meta['pool_file']
        ).split('.')[0]
        self.data['constraints'] = '\n'.join(
            [cc for cc in self.constraints_to_text()]
        )
        self.data['asn_id'] = self.acid.id
        self.new_product(product_name=self.dms_product_name())

        # Parse out information from the pool file name.
        # Necessary to carry information to the Level3 output.
        parsed_name = re.search(_DMS_POOLNAME_REGEX, self.data['asn_pool'])
        if parsed_name is not None:
            pool_meta = {
                'program_id': parsed_name.group(1),
                'version': parsed_name.group(2)
            }
            self.meta['pool_meta'] = pool_meta

    def _add(self, member):
        """Add member to this association."""
        try:
            exposerr = member['exposerr']
        except KeyError:
            exposerr = None
        entry = {
            'expname': Utility.rename_to_level2b(member['filename']),
            'exptype': self.get_exposure_type(member),
            'exposerr': exposerr,
            'asn_candidate': member['asn_candidate']
        }

        self.update_validity(entry)
        members = self.current_product['members']
        members.append(entry)
        if exposerr not in _EMPTY:
            logger.warn('Member {} has error "{}"'.format(
                member['filename'],
                exposerr
            ))
            self.data['degraded_status'] = _DEGRADED_STATUS_NOTOK

        # Add entry to the short list
        self.members.add(entry[KEY])

    def _get_target(self):
        """Get string representation of the target

        Returns
        -------
        target: str
            The Level3 Product name representation
            of the target or source ID.
        """
        try:
            target = 's{0:0>5s}'.format(self.data['source_id'])
        except KeyError:
            target = 't{0:0>3s}'.format(self.data['target'])
        return target

    def _get_instrument(self):
        """Get string representation of the instrument

        Returns
        -------
        instrument: str
            The Level3 Product name representation
            of the instrument
        """
        instrument = self.constraints['instrument']['value']
        return instrument

    def _get_opt_element(self):
        """Get string representation of the optical elements

        Returns
        -------
        opt_elem: str
            The Level3 Product name representation
            of the optical elements.
        """
        opt_elem = ''
        join_char = ''
        try:
            value = self.constraints['opt_elem']['value']
        except KeyError:
            pass
        else:
            if value not in _EMPTY and value != 'clear':
                opt_elem = value
                join_char = '-'
        try:
            value = self.constraints['opt_elem2']['value']
        except KeyError:
            pass
        else:
            if value not in _EMPTY and value != 'clear':
                opt_elem = join_char.join(
                    [opt_elem, value]
                )
        if opt_elem == '':
            opt_elem = 'clear'
        return opt_elem

    def _get_exposure(self):
        """Get string representation of the exposure id

        Returns
        -------
        exposure: str
            The Level3 Product name representation
            of the exposure & activity id.

        Raises
        ------
        AssociationNotAConstraint
            No constraints produce this value
        """
        try:
            activity_id = self.constraints['activity_id']['value']
        except KeyError:
            raise AssociationNotAConstraint
        else:
            if activity_id not in _EMPTY:
                exposure = '{0:0>2s}'.format(activity_id)
        return exposure

    def _add_items(self, items, product_name=None, with_type=False):
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

        with_type: bool
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
            type_ = 'science'
            if with_type:
                item, type_ = item
            entry = {
                'expname': item,
                'exptype': type_
            }
            self.update_validity(entry)
            members.append(entry)
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


class Utility(object):
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
    def rename_to_level2b(level1b_name):
        """Rename a Level 1b Exposure to another level

        Parameters
        ----------
        level1b_name: str
            The Level 1b exposure name.

        Returns
        -------
        str
            The Level 2b name
        """
        match = re.match(_LEVEL1B_REGEX, level1b_name)
        if match is None or match.group('type') != '_uncal':
            logger.warn((
                'Member FILENAME="{}" is not a Level 1b name. '
                'Cannot transform to Level 2b.'
            ).format(
                level1b_name
            ))
            return level1b_name

        level2b_name = ''.join([
            match.group('path'),
            '_cal',
            match.group('extension')
        ])
        return level2b_name

    @staticmethod
    def get_candidate_list(value):
        """Parse the candidate list from a member value

        Parameters
        ----------
        value: str
            The value from the member to parse. Usually
            member['ASN_CANDIDATE']

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
                'inputs': ['targetid']
            },
        })

        # Check and continue initialization.
        super(AsnMixin_Target, self).__init__(*args, **kwargs)


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


class AsnMixin_Image(DMS_Level3_Base):
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

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'image3'
        super(AsnMixin_Image, self)._init_hook(member)


class AsnMixin_Spectrum(DMS_Level3_Base):
    """All things that are spectrum"""

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'spec3'
        super(AsnMixin_Spectrum, self)._init_hook(member)


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
