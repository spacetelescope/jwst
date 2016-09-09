"""Base classes which define the Level3 Associations"""
from collections import namedtuple
import logging
from os.path import basename
import re

from jwst.associations import (
    Association,
    libpath
)
from jwst.associations.exceptions import AssociationNotAConstraint
from jwst.associations.lib.counter import Counter

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Start of the discovered association ids.
_DISCOVERED_ID_START = 3001

# Non-specified values found in DMS Association Pools
_EMPTY = (None, 'NULL', 'CLEAR')

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
_ASN_NAME_TEMPLATE = 'jw{program}-{acid}_{stamp}_{type}_{sequence:03d}_asn'
_LEVEL1B_REGEX = '(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'
_DMS_POOLNAME_REGEX = 'jw(\d{5})_(\d{8}[Tt]\d{6})_pool'

# Product name regex's
_REGEX_ACID_VALUE = '(o\d{3}|(c|a)\d{4})'
_REGEX_ACID_CONSTRAINT = '\(.*\'(?P<id>[a-z]\d{3,4})\'.*\,.*\'(?P<type>\w+)\'.*\)'


# Define the association's association candidate id.
ACID = namedtuple('ACID', ['id', 'type'])


class DMS_Level3_Base(Association):
    """Basic class for DMS Level3 associations."""

    def __init__(self, *args, **kwargs):

        self.candidates = set()

        # Initialize discovered association ID
        self.discovered_id = Counter(_DISCOVERED_ID_START)

        # Let us see if member belongs to us.
        super(DMS_Level3_Base, self).__init__(*args, **kwargs)

    @property
    def acid(self):
        """Association ID"""
        for _, constraint in self.constraints.items():
            if constraint.get('is_acid', False):
                value = re.sub('\\\\', '', constraint['value'])
                m = re.search(_REGEX_ACID_CONSTRAINT, value)
                acid = ACID(
                    id=m.groupdict()['id'],
                    type=m.groupdict()['type']
                )
                break
        else:
            id = 'a{:0>3}'.format(self.discovered_id.value)
            acid = ACID(id=id, type='DISCOVERED')

        return acid

    @property
    def asn_name(self):
        program = self.data['program']
        timestamp = self.timestamp
        asn_type = self.data['asn_type']
        sequence = self.sequence

        name = _ASN_NAME_TEMPLATE.format(
            program=program,
            acid=self.acid.id,
            stamp=timestamp,
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
        target = self._get_target_id()

        instrument = self._get_instrument()

        opt_elem = self._get_opt_element()

        try:
            exposure = self._get_exposure()
        except AssociationNotAConstraint:
            exposure = ''
        else:
            exposure = '-' + exposure

        product_name = 'jw{}-{}_{}_{}_{}_{{product_type}}{}.fits'.format(
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

        self.schema_file = ASN_SCHEMA
        self.data['targname'] = member['TARGETID']
        self.data['program'] = str(member['PROGRAM'])
        self.data['asn_pool'] = basename(
            member.meta['pool_file']
        ).split('.')[0]
        self.data['constraints'] = '\n'.join(
            [cc for cc in self.constraints_to_text()]
        )
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
            'exptype': member['PNTGTYPE'],
            'exposerr': exposerr,
            'asn_candidate': member['ASN_CANDIDATE']
        }
        members = self.current_product['members']
        members.append(entry)
        self.candidates.add(entry['asn_candidate'])
        self.data['degraded_status'] = _DEGRADED_STATUS_OK
        if exposerr not in _EMPTY:
            self.data['degraded_status'] = _DEGRADED_STATUS_NOTOK
            logger.warn('Member {} has error "{}"'.format(
                member['FILENAME'],
                exposerr
            ))

    def _get_target_id(self):
        """Get string representation of the target

        Returns
        -------
        target_id: str
            The Level3 Product name representation
            of the target or source ID.
        """
        try:
            target = 's{0:0>5s}'.format(self.data['source_id'])
        except KeyError:
            target = 't{0:0>3s}'.format(self.data['targname'])
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
        if self.constraints['opt_elem']['value'] not in _EMPTY:
            opt_elem = self.constraints['opt_elem']['value']
            join_char = '-'
        if self.constraints['opt_elem2']['value'] not in _EMPTY:
            opt_elem = join_char.join(
                [opt_elem, self.constraints['opt_elem2']['value']]
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
    def filter_cross_candidates(associations):
        """Return only those associations that have multiple candidates

        Parameters
        ----------
        associations: iterable
            The list of associations to check. The list
            is that returned by the `generate` function.

        Returns
        -------
        iterable
            The new list of just cross candidate associations.
        """
        result = [
            asn
            for asn in associations
            if len(asn.candidates) > 1
        ]
        return result

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

# ---------------------------------------------
# Mixins to define the broad category of rules.
# ---------------------------------------------


class AsnMixin_Unique_Config(DMS_Level3_Base):
    """Restrict to unique insturment configuration"""

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
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER']
            },
            'opt_elem2': {
                'value': None,
                'inputs': ['PUPIL']
            },
            'detector': {
                'value': '(?!NULL).+',
                'inputs': ['DETECTOR']
            },
        })

        super(AsnMixin_Unique_Config, self).__init__(*args, **kwargs)


class AsnMixin_Target(DMS_Level3_Base):
    """Constrain by Target"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'target_name': {
                'value': None,
                'inputs': ['TARGETID']
            },
        })

        # Check and continue initialization.
        super(AsnMixin_Target, self).__init__(*args, **kwargs)


class AsnMixin_MIRI(AsnMixin_Unique_Config):
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


class AsnMixin_NIRSPEC(AsnMixin_Unique_Config):
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


class AsnMixin_NIRISS(AsnMixin_Unique_Config):
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


class AsnMixin_NIRCAM(AsnMixin_Unique_Config):
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


class AsnMixin_Spectrum(AsnMixin_Unique_Config):
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
