"""Association Definitions: DMS Level2b product associations
"""
import logging
from os.path import basename
import re

from jwst.associations import (
    Association,
    libpath
)
from jwst.associations.lib.rules_level3_base import Utility as Utility_Level3

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# The schema that these associations must adhere to.
ASN_SCHEMA = libpath('asn_schema_jw_level2b.json')

# File templates
_DMS_POOLNAME_REGEX = 'jw(\d{5})_(\d{3})_(\d{8}[Tt]\d{6})_pool'
_LEVEL1B_REGEX = '(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'


class DMS_Level2b_Base(Association):
    """Basic class for DMS Level3 associations."""

    def __init__(self, *args, **kwargs):

        # I am defined by the following constraints
        self.add_constraints({
            'program': {
                'value': None,
                'inputs': ['PROGRAM']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER']
            },
            'detector': {
                'value': '(?!NULL).+',
                'inputs': ['DETECTOR']
            },
            'target': {
                'value': None,
                'inputs': ['TARGETID']
            },
        })

        # Now, lets see if member belongs to us.
        super(DMS_Level2b_Base, self).__init__(*args, **kwargs)

    @property
    def asn_name(self):
        template = 'jw_level2b_{}-{:03d}_{}_{}_asn'
        name = template.format(
            self.data['program'],
            self.sequence,
            self.data['asn_type'],
            self.version_id
        )
        return name.lower()

    @property
    def current_group(self):
        return self.data['groups'][-1]

    def new_group(self, member):
        """Start a new product"""
        group = {
            'members': []
        }
        try:
            self.data['groups'].append(group)
        except KeyError:
            self.data['groups'] = [group]

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""
        self.schema_file = ASN_SCHEMA
        self.data['target'] = member['TARGETID']
        self.data['program'] = str(member['PROGRAM'])
        self.data['asn_pool'] = basename(member.meta['pool_file']).split('.')[0]
        self.data['constraints'] = '\n'.join([cc for cc in self.constraints_to_text()])
        self.new_group(member)

    def _add(self, member):
        """Add member to this association."""
        entry = {
            'expname': Utility.rename_to_level2a(member['FILENAME']),
            'exptype': Utility.get_exposure_type(member, default='SCIENCE')
        }
        members = self.current_group['members']
        members.append(entry)

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
        result.append('\n    Groups:')
        for group in self.data['groups']:
            result.append('        Group:')
            for member in group['members']:
                result.append('            {:s}: {:s}'.format(member['expname'], member['exptype']))

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
            The Level 2b name
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


class Asn_MIRI_LRS_BKGNOD(DMS_Level2b_Base):
    """MIRI Background Nodding for LRS Fixed-slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'instrument': {
                'value': 'MIRI',
                'inputs': ['INSTRUME']
            },
            'patttype': {
                'value': 'POINT_SOURCE',
                'inputs': ['PATTTYPE']
            },
            'exp_type': {
                'value': 'MIR_LRS-FIXEDSLIT',
                'inputs': ['EXP_TYPE']
            },
        })

        # Check and continue initialization.
        super(Asn_MIRI_LRS_BKGNOD, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        self.data['asn_type'] = 'bkgnod'
        super(Asn_MIRI_LRS_BKGNOD, self)._init_hook(member)
