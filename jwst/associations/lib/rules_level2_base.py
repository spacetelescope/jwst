"""Base classes which define the Level2 Associations"""
import logging
from os.path import (
    basename,
    splitext
)

import re

from jwst.associations import (
    Association,
    AssociationRegistry,
    libpath
)
from jwst.associations.lib.dms_base import (DMSBaseMixin, PRODUCT_NAME_DEFAULT)
from jwst.associations.lib.rules_level3_base import _EMPTY
from jwst.associations.lib.rules_level3_base import Utility as Utility_Level3

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = [
    'ASN_SCHEMA',
    'AsnMixin_Lv2Image',
    'AsnMixin_Lv2ImageNonScience',
    'AsnMixin_Lv2ImageScience',
    'AsnMixin_Lv2Mode',
    'AsnMixin_Lv2Singleton',
    'AsnMixin_Lv2Spec',
    'AsnMixin_Lv2SpecScience',
    'AsnMixin_Lv2Special',
    'DMSLevel2bBase',
    'Utility'
]

# The schema that these associations must adhere to.
ASN_SCHEMA = libpath('asn_schema_jw_level2b.json')

# Flag to exposure type
FLAG_TO_EXPTYPE = {
    'background': 'background',
}

# File templates
_DMS_POOLNAME_REGEX = 'jw(\d{5})_(\d{3})_(\d{8}[Tt]\d{6})_pool'
_LEVEL1B_REGEX = '(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'
_REGEX_LEVEL2A = '(?P<path>.+)(?P<type>_rate(ints)?)'

# Key that uniquely identfies items.
KEY = 'expname'

# Exposures that are always TSO
TSO_EXP_TYPES = (
    'mir_lrs-slitless',
    'nis_soss',
    'nrc_tsimage',
    'nrc_tsgrism',
    'nrs_brightobj'
)

# Exposures that get Level2b processing
IMAGE2_SCIENCE_EXP_TYPES = [
    'mir_image',
    'mir_lyot',
    'mir_4qpm',
    'nis_ami',
    'nis_image',
    'nrc_image',
    'nrc_coron',
    'nrc_tsimage',
]

IMAGE2_NONSCIENCE_EXP_TYPES = [
    'mir_coroncal',
    'mir_tacq',
    'nis_focus',
    'nis_tacq',
    'nis_taconfirm',
    'nrc_tacq',
    'nrc_taconfirm',
    'nrc_focus',
    'nrs_bota',
    'nrs_confirm',
    'nrs_focus',
    'nrs_image',
    'nrs_mimf',
    'nrs_taslit',
    'nrs_tacq',
    'nrs_taconfirm',
]

SPEC2_SCIENCE_EXP_TYPES = [
    'nrc_grism',
    'nrc_tsgrism',
    'mir_lrs-fixedslit',
    'mir_lrs-slitless',
    'mir_mrs',
    'nrs_fixedslit',
    'nrs_ifu',
    'nrs_msaspec',
    'nrs_brightobj',
    'nis_soss',
]


class DMSLevel2bBase(DMSBaseMixin, Association):
    """Basic class for DMS Level2 associations."""

    # Set the validation schema
    schema_file = ASN_SCHEMA

    # Attribute values that are indicate the
    # attribute is not specified.
    INVALID_VALUES = _EMPTY

    def __init__(self, *args, **kwargs):

        super(DMSLevel2bBase, self).__init__(*args, **kwargs)

        # Keep the set of members included in this association
        self.members = set()

        # I am defined by the following constraints
        self.add_constraints({
            'program': {
                'value': None,
                'inputs': ['program']
            },
            'is_tso': {
                'value': None,
                'inputs': ['tsovisit'],
                'required': False,
                'force_unique': True,
            }
        })

        # Initialize validity checks
        self.validity.update({
            'has_science': {
                'validated': False,
                'check': lambda member: member['exptype'] == 'science'
            }
        })

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

    def has_science(self, item):
        """Only allow a single science in the association

        Parameters
        ----------
        item: dict
            The item in question

        Returns
        -------
        bool
            True if item can be added
        """
        exptype = self.get_exposure_type(item)
        limit_reached = len(self.members_by_type('science')) >= 1
        limit_reached = limit_reached and exptype == 'science'
        return limit_reached

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
            science = self.members_by_type('science')[0]
        except IndexError:
            return PRODUCT_NAME_DEFAULT

        try:
            science_path, ext = splitext(science['expname'])
        except Exception:
            return PRODUCT_NAME_DEFAULT

        match = re.match(_REGEX_LEVEL2A, science_path)
        if match:
            return match.groupdict()['path']
        else:
            return science_path

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
        is_tso = self.constraints['is_tso']['value'] == 't'
        if not is_tso:
            is_tso = item['exp_type'] in TSO_EXP_TYPES
        member = {
            'expname': Utility.rename_to_level2a(
                item['filename'], is_tso=is_tso
            ),
            'exptype': self.get_exposure_type(item)
        }
        return member

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""
        self.data['target'] = item['targetid']
        self.data['program'] = '{:0>5s}'.format(item['program'])
        self.data['asn_pool'] = basename(
            item.meta['pool_file']
        ).split('.')[0]
        self.data['constraints'] = '\n'.join(
            [cc for cc in self.constraints_to_text()]
        )
        self.data['asn_id'] = self.acid.id
        self.new_product(self.dms_product_name())

    def _add(self, item):
        """Add item to this association.

        Parameters
        ----------
        item: dict
            The item to be adding.
        """
        member = self.make_member(item)
        members = self.current_product['members']
        members.append(member)
        self.update_validity(member)

        # Add member to the short list
        self.members.add(member[KEY])

        # Update association state due to new member
        self.update_asn()

    def _add_items(self, items, meta=None, product_name_func=None, **kwargs):
        """Force adding items to the association

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

        product_name_func: func
            Used if product name is 'undefined' using
            the class's procedures.

        Notes
        -----
        This is a low-level shortcut into adding members, such as file names,
        to an association. All defined shortcuts and other initializations are
        by-passed, resulting in a potentially unusable association.

        `product_name_func` is used to define the product names instead of
        the default methods. The call signature is:

            product_name_func(item, idx)

        where `item` is each item being added and `idx` is the count of items.

        """
        if meta is None:
            meta = {}
        for idx, item in enumerate(items, start=1):
            self.new_product()
            members = self.current_product['members']
            member = {
                'expname': item,
                'exptype': 'science'
            }
            members.append(member)
            self.update_validity(member)
            self.update_asn()

            # If a product name function is given, attempt
            # to use.
            if product_name_func is not None:
                try:
                    self.current_product['name'] = product_name_func(item, idx)
                except Exception:
                    logger.debug(
                        'Attempted use of product_name_func failed.'
                        ' Default product name used.'
                    )

        self.data.update(meta)
        self.sequence = next(self._sequence)

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


class Utility():
    """Utility functions that understand DMS Level 3 associations"""

    @staticmethod
    def rename_to_level2a(level1b_name, is_tso=False):
        """Rename a Level 1b Exposure to another level

        Parameters
        ----------
        level1b_name: str
            The Level 1b exposure name.

        is_tso: boolean
            Use 'rateints' instead of 'rate' as
            the suffix.

        Returns
        -------
        str
            The Level 2a name
        """
        match = re.match(_LEVEL1B_REGEX, level1b_name)
        if match is None or match.group('type') != '_uncal':
            logger.warn((
                'Item FILENAME="{}" is not a Level 1b name. '
                'Cannot transform to Level 2a.'
            ).format(
                level1b_name
            ))
            return level1b_name

        suffix = 'rate'
        if is_tso:
            suffix = 'rateints'
        level2a_name = ''.join([
            match.group('path'),
            '_',
            suffix,
            match.group('extension')
        ])
        return level2a_name

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

        return finalized + lv2_asns

    @staticmethod
    def merge_asns(associations):
        """merge level2 associations

        Parameters
        ----------
        associations: [asn(, ...)]
            Associations to search for merging.

        Returns
        -------
        associatons: [association(, ...)]
            List of associations, some of which may be merged.
        """
        others = []
        lv2_asns = []
        for asn in associations:
            if isinstance(asn, DMSLevel2bBase):
                lv2_asns.append(asn)
            else:
                others.append(asn)

        lv2_asns = Utility._merge_asns(lv2_asns)

        return others + lv2_asns

    @staticmethod
    def _merge_asns(asns):
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
                        new_names = member_names.difference(
                            current_member_names
                        )
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

# General constraints
# -------------------
class AsnMixin_Lv2Mode(DMSLevel2bBase):
    """Fix the instrument configuration"""
    def __init__(self, *args, **kwargs):

        # I am defined by the following constraints
        self.add_constraints({
            'program': {
                'value': None,
                'inputs': ['program']
            },
            'target': {
                'value': None,
                'inputs': ['targetid'],
            },
            'instrument': {
                'value': None,
                'inputs': ['instrume']
            },
            'detector': {
                'value': None,
                'inputs': ['detector']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['filter', 'band']
            },
            'opt_elem2': {
                'value': None,
                'inputs': ['pupil', 'grating'],
                'required': False,
            },
            'subarray': {
                'value': None,
                'inputs': ['subarray'],
                'required': False,
            },
            'channel': {
                'value': None,
                'inputs': ['channel'],
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
                'inputs': lambda item: str(
                    self.has_science(item)
                ),
            }
        })

        # Now, lets see if item belongs to us.
        super(AsnMixin_Lv2Singleton, self).__init__(*args, **kwargs)


class AsnMixin_Lv2Special(DMSLevel2bBase):
    """Process special exposures as science

    Spectral exposures that are marked as backgrounds, imprints, etc.,
    still get 2b processing just as normal science. However, no other
    exposures should get included into the association.

    """
    def __init__(self, *args, **kwargs):
        self.add_constraints({
            'is_special': {
                'value': None,
                'inputs': [
                    'bkgdtarg',
                    'is_imprt',
                    'is_psf'
                ],
                'force_unique': False,
            }
        })

        super(AsnMixin_Lv2Special, self).__init__(*args, **kwargs)

    def get_exposure_type(self, item, default='science'):
        """Override to force exposure type to always be science

        Parameters
        ----------
        item: dict
            The pool entry to determine the exposure type of

        default: str or None
            The default exposure type.
            If None, routine will raise LookupError

        Returns
        -------
        exposure_type: 'science'
            Always returns as science
        """
        return 'science'


# Image-like Exposures
# --------------------
class AsnMixin_Lv2Image(DMSLevel2bBase):
    """Level 2 Image association base"""

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(AsnMixin_Lv2Image, self)._init_hook(item)
        self.data['asn_type'] = 'image2'


class AsnMixin_Lv2ImageScience(DMSLevel2bBase):
    """Level 2 Image association base"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': '|'.join(IMAGE2_SCIENCE_EXP_TYPES),
                'inputs': ['exp_type'],
                'force_unique': True,
            }
        })

        super(AsnMixin_Lv2ImageScience, self).__init__(*args, **kwargs)


class AsnMixin_Lv2ImageNonScience(DMSLevel2bBase):
    """Process selected non-science exposures

    Exposures, such as target acquisitions,
    though considered non-science, still get 2b processing.

    """
    def __init__(self, *args, **kwargs):
        self.add_constraints({
            'non_science': {
                'value': '|'.join(IMAGE2_NONSCIENCE_EXP_TYPES),
                'inputs': ['exp_type'],
                'force_unique': False,
            }
        })

        super(AsnMixin_Lv2ImageNonScience, self).__init__(*args, **kwargs)

    def get_exposure_type(self, item, default='science'):
        """Override to force exposure type to always be science

        Parameters
        ----------
        item: dict
            The pool entry to determine the exposure type of

        default: str or None
            The default exposure type.
            If None, routine will raise LookupError

        Returns
        -------
        exposure_type: 'science'
            Always returns as science
        """
        return 'science'


# Spectral-like exposures
# -----------------------
class AsnMixin_Lv2Spec(DMSLevel2bBase):
    """Level 2 Spectral association base"""

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super(AsnMixin_Lv2Spec, self)._init_hook(item)
        self.data['asn_type'] = 'spec2'


class AsnMixin_Lv2SpecScience(DMSLevel2bBase):
    """Level 2 Spectral association base"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': '|'.join(SPEC2_SCIENCE_EXP_TYPES),
                'inputs': ['exp_type'],
                'force_unique': True,
            }
        })

        super(AsnMixin_Lv2SpecScience, self).__init__(*args, **kwargs)
