"""Base classes which define the Level3 Associations"""
from collections import defaultdict
import logging
from os.path import basename
import re

from jwst.associations import (
    Association,
    libpath
)
from jwst.associations.association import (
    evaluate,
    is_iterable
)
from jwst.associations.exceptions import AssociationNotAConstraint
from jwst.associations.lib.acid import ACID
from jwst.associations.lib.counter import Counter

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Start of the discovered association ids.
_DISCOVERED_ID_START = 3001

# Non-specified values found in DMS Association Pools
_EMPTY = (None, 'NULL')

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
_ASN_NAME_TEMPLATE_STAMP = 'jw{program}-{acid}_{stamp}_{type}_{sequence:03d}_asn'
_ASN_NAME_TEMPLATE = 'jw{program}-{acid}_{type}_{sequence:03d}_asn'
_LEVEL1B_REGEX = '(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'
_DMS_POOLNAME_REGEX = 'jw(\d{5})_(\d{8}[Tt]\d{6})_pool'

# Product name regex's
_REGEX_ACID_VALUE = '(o\d{3}|(c|a)\d{4})'


# Key that uniquely identfies members.
KEY = 'expname'

# Exposure EXP_TYPE to Association EXPTYPE mapping
_EXPTYPE_MAP = {
    'MIR_TACQ':      'TARGET_ACQUISTION',
    'NIS_TACQ':      'TARGET_ACQUISTION',
    'NIS_TACONFIRM': 'TARGET_ACQUISTION',
    'NRC_TACQ':      'TARGET_ACQUISTION',
    'NRC_TACONFIRM': 'TARGET_ACQUISTION',
    'NRS_AUTOFLAT':  'AUTOFLAT',
    'NRS_AUTOWAVE':  'AUTOWAVE',
    'NRS_CONFIRM':   'TARGET_ACQUISTION',
    'NRS_TACQ':      'TARGET_ACQUISTION',
    'NRS_TACONFIRM': 'TARGET_ACQUISTION',
    'NRS_TASLIT':    'TARGET_ACQUISTION',
}


class DMS_Level3_Base(Association):
    """Basic class for DMS Level3 associations."""

    # Attribute values that are indicate the
    # attribute is not specified.
    INVALID_VALUES = _EMPTY

    # Make sequences type-dependent
    _sequences = defaultdict(Counter)

    def __init__(self, *args, **kwargs):

        # Keep the set of members included in this association
        self.members = set()

        # Initialize discovered association ID
        self.discovered_id = Counter(_DISCOVERED_ID_START)

        # Initialize validity checks
        self.validity = {
            'has_science': {
                'validated': False,
                'check': lambda entry: entry['exptype'] == 'SCIENCE'
            }
        }

        # Let us see if member belongs to us.
        super(DMS_Level3_Base, self).__init__(*args, **kwargs)

    @property
    def is_valid(self):
        return all(test['validated'] for test in self.validity.values())

    @property
    def acid(self):
        """Association ID"""
        for _, constraint in self.constraints.items():
            if constraint.get('is_acid', False):
                value = re.sub('\\\\', '', constraint['value'])
                try:
                    acid = ACID(value)
                except ValueError:
                    pass
                else:
                    break
        else:
            id = 'a{:0>3}'.format(self.discovered_id.value)
            acid = ACID((id, 'DISCOVERED'))

        return acid

    @property
    def asn_name(self):
        program = self.data['program']
        version_id = self.version_id
        asn_type = self.data['asn_type']
        sequence = self.sequence

        if version_id:
            name = _ASN_NAME_TEMPLATE_STAMP.format(
                program=program,
                acid=self.acid.id,
                stamp=version_id,
                type=asn_type,
                sequence=sequence,
            )
        else:
            name = _ASN_NAME_TEMPLATE.format(
                program=program,
                acid=self.acid.id,
                type=asn_type,
                sequence=sequence,
            )
        return name.lower()

    @property
    def current_product(self):
        return self.data['products'][-1]

    def new_product(self, member):
        """Start a new product"""
        product = {
            'name': self.product_name(),
            'members': []
        }
        try:
            self.data['products'].append(product)
        except KeyError:
            self.data['products'] = [product]

    def product_name(self):
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

    def update_validity(self, entry):
        for test in self.validity.values():
            if not test['validated']:
                test['validated'] = test['check'](entry)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""
        super(DMS_Level3_Base, self)._init_hook(member)

        # Set which sequence counter should be used.
        self._sequence = self._sequences[self.data['asn_type']]

        self.schema_file = ASN_SCHEMA
        self.data['target'] = member['TARGETID']
        self.data['program'] = str(member['PROGRAM'])
        self.data['asn_pool'] = basename(
            member.meta['pool_file']
        ).split('.')[0]
        self.data['constraints'] = '\n'.join(
            [cc for cc in self.constraints_to_text()]
        )
        self.data['asn_id'] = self.acid.id
        self.new_product(member)

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
            exposerr = member['EXPOSERR']
        except KeyError:
            exposerr = None
        entry = {
            'expname': Utility.rename_to_level2b(member['FILENAME']),
            'exptype': Utility.get_exposure_type(member, default='SCIENCE'),
            'exposerr': exposerr,
            'asn_candidate': member['ASN_CANDIDATE']
        }

        self.update_validity(entry)
        members = self.current_product['members']
        members.append(entry)
        self.data['degraded_status'] = _DEGRADED_STATUS_OK
        if exposerr not in _EMPTY:
            self.data['degraded_status'] = _DEGRADED_STATUS_NOTOK
            logger.warn('Member {} has error "{}"'.format(
                member['FILENAME'],
                exposerr
            ))

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
            if value not in _EMPTY and value != 'CLEAR':
                opt_elem = value
                join_char = '-'
        try:
            value = self.constraints['opt_elem2']['value']
        except KeyError:
            pass
        else:
            if value not in _EMPTY and value != 'CLEAR':
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

    def __repr__(self):
        file_name, json_repr = self.to_json()
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
        counters = defaultdict(lambda : defaultdict(Counter))
        for asn in associations:
            asn.sequence = next(counters[asn.data['asn_id']][asn.data['asn_type']])

    @staticmethod
    def filter_discovered_only(
            associations,
            discover_ruleset,
            candidate_ruleset,
            keep_candidates=True,
    ):
        """Return only those associations that have multiple candidates

        Parameters
        ----------
        associations: iterable
            The list of associations to check. The list
            is that returned by the `generate` function.

        discover_ruleset: str
            The name of the ruleset that has the discover rules

        candidate_ruleset: str
            The name of the ruleset that finds just candidates

        keep_candidates: bool
            Keep explicit candidate associations in the list.

        Returns
        -------
        iterable
            The new list of just cross candidate associations.

        Notes
        -----
        This utility is only meant to run on associations that have
        been constructed. Associations that have been Association.dump
        and then Association.load will not return proper results.
        """
        # Split the associations along discovered/not discovered lines
        asn_by_ruleset = {
            candidate_ruleset: [],
            discover_ruleset: []
        }
        for asn in associations:
            asn_by_ruleset[asn.registry.name].append(asn)
        candidate_list = asn_by_ruleset[candidate_ruleset]
        discover_list = asn_by_ruleset[discover_ruleset]

        # Filter out the non-unique discovereds.
        for candidate in candidate_list:
            if len(discover_list) == 0:
                break
            unique_list = []
            for discover in discover_list:
                if discover.data['asn_type'] == candidate.data['asn_type'] and \
                   discover.members == candidate.members:
                    # This association is not unique. Ignore
                    pass
                else:
                    unique_list.append(discover)

            # Reset the discovered list to the new unique list
            # and try the next candidate.
            discover_list = unique_list

        if keep_candidates:
            discover_list.extend(candidate_list)
        return discover_list

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
    def get_exposure_type(member, default=None):
        """Determine the exposure type of a pool member

        Parameters
        ----------
        member: dict
            The pool entry to determine the exposure type of

        default: str or None
            The default exposure type.
            If None, routine will raise LookupError

        Returns
        -------
        exposure_type: str
            Exposure type. Can be one of
                'SCIENCE': Member contains science data
                'TARGET_AQUISITION': Member contains target acquisition data.
                'AUTOFLAT': NIRSpec AUTOFLAT
                'AUTOWAVE': NIRSpec AUTOWAVE

        Raises
        ------
        LookupError
            When `default` is None and an exposure type cannot be determined
        """
        result = default
        try:
            exp_type = member['EXP_TYPE']
        except KeyError:
            raise LookupError('Exposure type cannot be determined')

        result = _EXPTYPE_MAP.get(exp_type, default)
        return result

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
                'inputs': ['PROGRAM']
            },
            'instrument': {
                'value': None,
                'inputs': ['INSTRUME']
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
                'inputs': ['FILTER']
            },
            'opt_elem2': {
                'value': None,
                'inputs': ['PUPIL', 'GRATING'],
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
                'inputs': ['TARGETID']
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
                'value': 'MIRI',
                'inputs': ['INSTRUME']
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
                'value': 'NIRSPEC',
                'inputs': ['INSTRUME']
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
                'value': 'NIRISS',
                'inputs': ['INSTRUME']
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
                'value': 'NIRCAM',
                'inputs': ['INSTRUME']
            },
        })

        # Check and continue initialization.
        super(AsnMixin_NIRCAM, self).__init__(*args, **kwargs)


class AsnMixin_Image(DMS_Level3_Base):
    """All things that are in imaging mode"""

    def __init__(self, *args, **kwargs):

        self.add_constraints({
            'exp_type': {
                'value': 'NRC_IMAGE|MIR_IMAGE|NIS_IMAGE|FGS_IMAGE',
                'inputs': ['EXP_TYPE'],
                'force_unique': True,
            }
        })

        super(AsnMixin_Image, self).__init__(*args, **kwargs)


class AsnMixin_Spectrum(DMS_Level3_Base):
    """All things that are spectrum"""

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'spec'
        super(AsnMixin_Spectrum, self)._init_hook(member)


class AsnMixin_CrossCandidate(DMS_Level3_Base):
    """Basic constraints for Cross-Candidate associations"""

    def is_valid(self):
        candidates = set(
            member['asn_candidate_id']
            for product in self.data['products']
            for member in product['members']
        )
        return len(candidates) > 1
