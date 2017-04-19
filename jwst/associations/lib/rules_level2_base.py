"""Base classes which define the Level2 Associations"""
import logging
from os.path import basename

import re

from jwst.associations import (
    Association,
    AssociationRegistry,
    libpath
)
from jwst.associations.association import getattr_from_list
from jwst.associations.lib.dms_base import (DMSBaseMixin, PRODUCT_NAME_DEFAULT)
from jwst.associations.lib.rules_level3_base import _EMPTY
from jwst.associations.lib.rules_level3_base import Utility as Utility_Level3

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = [
    'ASN_SCHEMA',
    'AsnMixin_Lv2Bkg',
    'AsnMixin_Lv2Image',
    'AsnMixin_Lv2Mode',
    'AsnMixin_Lv2Singleton',
    'AsnMixin_Lv2Spec',
    'DMSLevel2bBase',
    'Utility'
]

# The schema that these associations must adhere to.
ASN_SCHEMA = libpath('asn_schema_jw_level2b.json')

# Flag to exposure type
FLAG_TO_EXPTYPE = {
    'BACKGROUND': 'BACKGROUND',
}

# File templates
_DMS_POOLNAME_REGEX = 'jw(\d{5})_(\d{3})_(\d{8}[Tt]\d{6})_pool'
_LEVEL1B_REGEX = '(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'
_REGEX_LEVEL2A = '(?P<path>.+)(?P<type>_rate(ints)?)(?P<extension>\..+)'

# Key that uniquely identfies members.
KEY = 'expname'


class DMSLevel2bBase(DMSBaseMixin, Association):
    """Basic class for DMS Level2 associations."""

    # Set the validation schema
    schema_file = ASN_SCHEMA

    # Attribute values that are indicate the
    # attribute is not specified.
    INVALID_VALUES = _EMPTY

    def __init__(self, *args, **kwargs):

        # Keep the set of members included in this association
        self.members = set()

        # I am defined by the following constraints
        self.add_constraints({
            'program': {
                'value': None,
                'inputs': ['PROGRAM']
            },
        })

        # Initialize validity checks
        self.validity.update({
            'has_science': {
                'validated': False,
                'check': lambda entry: entry['exptype'] == 'SCIENCE'
            }
        })

        # Now, lets see if member belongs to us.
        super(DMSLevel2bBase, self).__init__(*args, **kwargs)

    def members_by_type(self, member_type):
        """Get list of members by their exposure type"""
        member_type = member_type.lower()
        try:
            members = self.current_product['members']
        except KeyError:
            result = []
        else:
            result = [
                member
                for member in members
                if member_type == member['exptype'].lower()
            ]

        return result

    def has_science(self, member, check_flags=None):
        """Only allow a single science in the association

        Parameters
        ----------
        member: dict
            The member in question

        check_flags: None or [key[, ...]]
            A list of extra keys to check for truthness in the member

        Returns
        -------
        bool
            True if member can be added
        """
        exptype = self.get_exptype(member, check_flags=check_flags)
        limit_reached = len(self.members_by_type('SCIENCE')) >= 1
        limit_reached = limit_reached and exptype == 'SCIENCE'
        return limit_reached

    def get_exptype(self, member, check_flags=None):
        """Get the exposure type for member

        Parameters
        ----------
        member: dict
            The member to be adding.


        check_flags: None or [key[, ...]]
            A list of extra keys to check for truthness in the member

        Returns
        -------
        exptype: str
        """
        exptype = Utility.get_exposure_type(member, default='SCIENCE')
        if check_flags:
            for flag in check_flags:
                try:
                    getattr_from_list(member, [flag], self.INVALID_VALUES)
                except KeyError:
                    continue
                else:
                    exptype = FLAG_TO_EXPTYPE[flag]
        return exptype

    def __eq__(self, other):
        """Compare equality of two assocaitions"""
        if isinstance(other, DMSLevel2bBase):
            result = self.data['asn_type'] == other.data['asn_type']
            result = result and (self.members == other.members)
            return result
        else:
            return NotImplemented

    def __ne__(self, other):
        """Compare inequality of two associations"""
        if isinstance(other, DMSLevel2bBase):
            return not self.__eq__(other)
        else:
            return NotImplemented

    def dms_product_name(self):
        """Define product name."""
        try:
            science = self.members_by_type('SCIENCE')[0]
        except IndexError:
            return PRODUCT_NAME_DEFAULT

        science_name = basename(science['expname']).lower()
        match = re.match(_REGEX_LEVEL2A, science_name)
        if match:
            return match.groupdict()['path']
        else:
            return science_name

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""
        self.data['target'] = member['TARGETID']
        self.data['program'] = str(member['PROGRAM'])
        self.data['asn_pool'] = basename(
            member.meta['pool_file']
        ).split('.')[0]
        self.data['constraints'] = '\n'.join(
            [cc for cc in self.constraints_to_text()]
        )
        self.data['asn_id'] = self.acid.id
        self.new_product(self.dms_product_name())

    def _add(self, member, check_flags=None):
        """Add member to this association.

        Parameters
        ----------
        member: dict
            The member to be adding.

        check_flags: None or [key[, ...]]
            A list of extra keys to check for truthness in the member
        """
        entry = {
            'expname': Utility.rename_to_level2a(member['FILENAME']),
            'exptype': self.get_exptype(member, check_flags=check_flags)
        }
        members = self.current_product['members']
        members.append(entry)
        self.update_validity(entry)

        # Add entry to the short list
        self.members.add(entry[KEY])

        # Update association state due to new member
        self.update_asn()

    def _add_items(self, items, meta=None):
        """ Force adding items to the association

        Parameters
        ----------
        items: [object[, ...]]
            A list of items to make members of the association.

        meta: dict
            A dict to be merged into the association meta information.
            The following are suggested to be assigned:
                - `asn_type`
                    The type of association.
                - `asn_rule`
                    The rule which created this association.
                - `asn_pool`
                    The pool from which the exposures came from
                - `program`
                    Originating observing program

        Notes
        -----
        This is a low-level shortcut into adding members, such as file names,
        to an association. All defined shortcuts and other initializations are
        by-passed, resulting in a potentially unusable association.
        """
        if meta is None:
            meta = {}
        for item in items:
            self.new_product()
            members = self.current_product['members']
            entry = {
                'expname': item,
                'exptype': 'SCIENCE'
            }
            members.append(entry)
            self.update_validity(entry)
            self.update_asn()
        self.data.update(meta)

    def update_asn(self):
        """Update association info based on current members"""
        self.current_product['name'] = self.dms_product_name()

    def __repr__(self):
        try:
            file_name, json_repr = self.ioregistry['json'].dump(self)
        except:
            return str(self.__class__)
        return json_repr

    def __str__(self):
        """Create human readable version of the association
        """

        result = ['Association {:s}'.format(self.asn_name)]

        # Parameters of the association
        result.append('    Parameters:')
        result.append('        Product type: {:s}'.format(self.data['asn_type']))
        result.append('        Rule:         {:s}'.format(self.data['asn_rule']))
        result.append('        Program:      {:s}'.format(self.data['program']))
        result.append('        Target:       {:s}'.format(self.data['target']))
        result.append('        Pool:         {:s}'.format(self.data['asn_pool']))

        for cc in self.constraints_to_text():
            result.append('        {:s}'.format(cc))

        # Products of the assocation
        for product in self.data['products']:
            result.append(
                '\t{} with {} members'.format(
                    product['name'],
                    len(product['members'])
                )
            )

        # That's all folks
        result.append('\n')
        return '\n'.join(result)


class Utility(object):
    """Utility functions that understand DMS Level 3 associations"""

    @staticmethod
    def rename_to_level2a(level1b_name):
        """Rename a Level 1b Exposure to another level

        Parameters
        ----------
        level1b_name: str
            The Level 1b exposure name.

        Returns
        -------
        str
            The Level 2a name
        """
        match = re.match(_LEVEL1B_REGEX, level1b_name)
        if match is None or match.group('type') != '_uncal':
            logger.warn((
                'Member FILENAME="{}" is not a Level 1b name. '
                'Cannot transform to Level 2a.'
            ).format(
                level1b_name
            ))
            return level1b_name

        level2a_name = ''.join([
            match.group('path'),
            '_rate',
            match.group('extension')
        ])
        return level2a_name

    @staticmethod
    def get_exposure_type(*args, **kwargs):
        return Utility_Level3.get_exposure_type(*args, **kwargs)

    @staticmethod
    def resequence(*args, **kwargs):
        return Utility_Level3.resequence(*args, **kwargs)

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
        lv2_asns = []
        for asn in associations:
            if isinstance(asn, DMSLevel2bBase):

                # Check validity
                if asn.is_valid:
                    lv2_asns.append(asn)

            else:
                finalized.append(asn)

        # Merge all the associations into common types
        merged_asns = Utility.merge_asns(lv2_asns)

        # Merge lists and return
        return finalized + merged_asns

    @staticmethod
    def merge_asns(asns):
        # Merge all the associations into common types
        merged_by_type = {}
        for asn in asns:
            try:
                current_asn = merged_by_type[asn['asn_type']]
            except KeyError:
                merged_by_type[asn['asn_type']] = asn
                current_asn = asn
            for product in asn['products']:
                merge_occured = False
                for current_product in current_asn['products']:
                    if product['name'] == current_product['name']:
                        member_names = set([
                            member['expname']
                            for member in product['members']
                        ])
                        current_member_names = [
                            member['expname']
                            for member in current_product['members']
                        ]
                        new_names = member_names.difference(current_member_names)
                        new_members = [
                            member
                            for member in product['members']
                            if member['expname'] in new_names
                        ]
                        current_product['members'].extend(new_members)
                        merge_occured = True
                if not merge_occured:
                    current_asn['products'].append(product)

        merged_asns = [
            asn
            for asn_type, asn in merged_by_type.items()
        ]
        return merged_asns


# ---------------------------------------------
# Mixins to define the broad category of rules.
# ---------------------------------------------
class AsnMixin_Lv2Image(DMSLevel2bBase):
    """Level 2 Image association base"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': (
                    'FGS_IMAGE'
                    '|MIR_IMAGE'
                    '|NIS_IMAGE'
                    '|NRC_IMAGE'
                ),
                'inputs': ['EXP_TYPE'],
                'force_unique': True,
            }
        })

        super(AsnMixin_Lv2Image, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        super(AsnMixin_Lv2Image, self)._init_hook(member)
        self.data['asn_type'] = 'image2'


class AsnMixin_Lv2Spec(DMSLevel2bBase):
    """Level 2 Spectral association base"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': (
                    'NRC_GRISM'
                    '|NRC_TSGRISM'
                    '|MIR_LRS-FIXEDSLIT'
                    '|MIR_LRS-SLITLESS'
                    '|MIR_MRS'
                    '|NRS_FIXEDSLIT'
                    '|NRS_IFU'
                    '|NRS_MSASPEC'
                    '|NRS_BRIGHTOBJ'
                    '|NIS_WFSS'
                    '|NIS_SOSS'
                ),
                'inputs': ['EXP_TYPE'],
                'force_unique': True,
            }
        })

        super(AsnMixin_Lv2Spec, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        super(AsnMixin_Lv2Spec, self)._init_hook(member)
        self.data['asn_type'] = 'spec2'


class AsnMixin_Lv2Mode(DMSLevel2bBase):
    """Fix the instrument configuration"""
    def __init__(self, *args, **kwargs):

        # I am defined by the following constraints
        self.add_constraints({
            'program': {
                'value': None,
                'inputs': ['PROGRAM']
            },
            'instrument': {
                'value': None,
                'inputs': ['INSTRUME']
            },
            'detector': {
                'value': None,
                'inputs': ['DETECTOR']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER', 'BAND']
            },
            'opt_elem2': {
                'value': None,
                'inputs': ['PUPIL', 'GRATING'],
                'required': False,
            },
            'subarray': {
                'value': None,
                'inputs': ['SUBARRAY'],
                'required': False,
            },
            'channel': {
                'value': None,
                'inputs': ['CHANNEL'],
                'required': False,
            }
        })

        super(AsnMixin_Lv2Mode, self).__init__(*args, **kwargs)


class AsnMixin_Lv2Singleton(DMSLevel2bBase):
    """Allow only single science exposure"""

    def __init__(self, *args, **kwargs):
        # I am defined by the following constraints
        self.add_constraints({
            'single_science': {
                'test': self.match_constraint,
                'value': 'False',
                'inputs': lambda member: str(
                    self.has_science(member)
                ),
            }
        })

        # Now, lets see if member belongs to us.
        super(AsnMixin_Lv2Singleton, self).__init__(*args, **kwargs)


class AsnMixin_Lv2Bkg(DMSLevel2bBase):
    """Acquire backgrounds"""

    def __init__(self, *args, **kwargs):

        # I am defined by the following constraints
        self.add_constraints({
            'background': {
                'inputs': ['ASN_CANDIDATE'],
                'value': '.+BACKGROUND.+',
                'force_unique': True,
                'is_acid': False,
                'required': False,
            },
            'single_science': {
                'test': self.match_constraint,
                'value': 'False',
                'inputs': lambda member: str(
                    self.has_science(member, check_flags=['BACKGROUND'])
                ),
            },
        })

        # Now, lets see if member belongs to us.
        super(AsnMixin_Lv2Bkg, self).__init__(*args, **kwargs)

    def _add(self, member, check_flags=None):
        if not check_flags:
            check_flags = []
        check_flags.append('BACKGROUND')
        super(AsnMixin_Lv2Bkg, self)._add(member, check_flags)
