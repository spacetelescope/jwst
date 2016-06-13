"""Association Definitions: DMS Level3 product associations
"""
import logging
from os.path import basename
import re

from jwst.associations.association import (
    Association,
    getattr_from_list
)

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# The schema that these associations must adhere to.
_ASN_SCHEMA_LEVEL3 = 'asn_schema_jw_level3.json'

_DMS_POOLNAME_REGEX = 'jw(\d{5})_(\d{8}[Tt]\d{6})_pool'
_EMPTY = (None, 'NULL', 'CLEAR')
_DEGRADED_STATUS_OK = (
    'No known degraded exposures in association.'
)
_DEGRADED_STATUS_NOTOK = (
    'One or more members have an error associated with them.'
    '\nDetails can be found in the members.exposerr property.'
)

_LEVEL1B_REGEX = '(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)'


class DMS_Level3_Base(Association):
    """Basic class for DMS Level3 associations."""

    def __init__(self, *args, **kwargs):

        self.candidates = set()

        # Let us see if member belongs to us.
        super(DMS_Level3_Base, self).__init__(*args, **kwargs)

    @property
    def asn_name(self):
        template = 'jw{}_{}_{}_{:03d}_asn'
        program = self.data['program']
        asn_candidate_ids = self.constraints.get('asn_candidate_ids', None)
        if asn_candidate_ids is not None:
            program = '-'.join([
                program,
                '{0:0>4s}'.format(asn_candidate_ids['value'])
            ])
        timestamp = self.timestamp
        asn_type = self.data['asn_type']
        sequence = self.sequence

        name = template.format(
            program,
            timestamp,
            asn_type,
            sequence,
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
        program = self.data['program']

        asn_candidate_ids = self.constraints.get('asn_candidate_ids', None)
        if asn_candidate_ids is not None:
            id = asn_candidate_ids['value']
            if len(id) <= 3:
                id = 'o{0:0>3s}'.format(id)
            program = '-'.join([
                program,
                '{}'.format(id)
            ])
        try:
            target = 's{0:0>5s}'.format(self.data['source_id'])
        except KeyError:
            target = 't{0:0>3s}'.format(self.data['targname'])

        instrument = self.constraints['instrument']['value']
        asn_type = self.data['asn_type']

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

        exposure = ''
        try:
            activity_id = self.constraints['activity_id']['value']
        except KeyError:
            pass
        else:
            if activity_id not in _EMPTY:
                exposure = '-{0:0>2s}'.format(activity_id)

        product_name = 'jw{}_{}_{}_{}_{}{}.fits'.format(
            program,
            target,
            instrument,
            opt_elem,
            asn_type,
            exposure
        )

        return product_name.lower()

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""
        super(DMS_Level3_Base, self)._init_hook(member)

        self.schema_file = _ASN_SCHEMA_LEVEL3
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
            'asn_candidate_id': getattr_from_list(
                member,
                ['ASN_CANDIDATE_ID', 'OBS_NUM']
            )[1]
        }
        members = self.current_product['members']
        members.append(entry)
        self.candidates.add(entry['asn_candidate_id'])
        self.data['degraded_status'] = _DEGRADED_STATUS_OK
        if exposerr not in _EMPTY:
            self.data['degraded_status'] = _DEGRADED_STATUS_NOTOK
            logger.warn('Member {} has error "{}"'.format(
                member['FILENAME'],
                exposerr
            ))

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
            'target_acq': {
                'value': '(TARGET_ACQUISITION)?',
                'inputs': ['PNTGTYPE'],
                'force_unique': False,
            }
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

        self.data['asn_type'] = 'spectrum'
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


# --------------------------------
# Start of the User-level rules
# --------------------------------

# ----------------------------------
# Associations defined by Candidates

class Asn_Mosaic(
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """Association Candidate type of Mosaic"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'asn_candidate_id': {
                'value': None,
                'inputs': ['ASN_CANDIDATE_ID'],
            },
            'asn_candidate_type': {
                'value': 'MOSAIC',
                'inputs': ['ASN_CANDIDATE_TYPE'],
            },
            'wfsvisit': {
                'value': 'NULL',
                'inputs': ['WFSVISIT']
            }
        })

        # Now check and continue initialization.
        super(Asn_Mosaic, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'mosaic'
        super(Asn_Mosaic, self)._init_hook(member)


# Image associations
class Asn_Dither(
        AsnMixin_Image,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """Non-Association Candidate Dither Associations"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'pointing_type': {
                'value': 'SCIENCE',
                'inputs': ['PNTGTYPE']
            },
            'wfsvisit': {
                'value': 'NULL',
                'inputs': ['WFSVISIT'],
            },
            'asn_candidate_type': {
                'value': 'OBSERVATION',
                'inputs': ['ASN_CANDIDATE_TYPE'],
                'required': False,
            },
        })

        # Now check and continue initialization.
        super(Asn_Dither, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'dither'
        super(Asn_Dither, self)._init_hook(member)


class Asn_WFSCMB(
        AsnMixin_Image,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """Wavefront Sensing association"""
    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'wfsvisit': {
                'value': '(?!NULL).+',
                'inputs': ['WFSVISIT'],
            },
            'asn_candidate_type': {
                'value': '(?!OBSERVATION).+',
                'inputs': ['ASN_CANDIDATE_TYPE']
            },
            'asn_candidate_id': {
                'value': None,
                'inputs': ['ASN_CANDIDATE_ID']
            },
            'activity_id': {
                'value': None,
                'inputs': ['ACT_ID']
            }
        })

        # Now check and continue initialization.
        super(Asn_WFSCMB, self).__init__(*args, **kwargs)

    def _init_hook(self, member):
        """Post-check and pre-add initialization"""

        self.data['asn_type'] = 'wfscmb'
        super(Asn_WFSCMB, self)._init_hook(member)


# Spectrographic Associations
class Asn_MIRI_LRS_FIXEDSLIT(
        AsnMixin_Spectrum,
        AsnMixin_MIRI,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """MIRI LRS Fixed slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'patttype': {
                'value': None,
                'inputs': ['PATTTYPE'],
                'force_unique': True
            },
            'exp_type': {
                'value': 'MIR_LRS-FIXEDSLIT',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': 'P750L',
                'inputs': ['FILTER']
            },
            'subarray': {
                'value': 'FULL',
                'inputs': ['SUBARRAY']
            }
        })

        # Check and continue initialization.
        super(Asn_MIRI_LRS_FIXEDSLIT, self).__init__(*args, **kwargs)


class Asn_MIRI_LRS_SLITLESS(
        AsnMixin_Spectrum,
        AsnMixin_MIRI,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """MIRI LRS Slitless"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'MIR_LRS-SLITLESS',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': 'P750L',
                'inputs': ['FILTER']
            },
            'subarray': {
                'value': 'SUBPRISM',
                'inputs': ['SUBARRAY']
            }
        })

        # Check and continue initialization.
        super(Asn_MIRI_LRS_SLITLESS, self).__init__(*args, **kwargs)


class Asn_NIR_SO_SLITLESS(
        AsnMixin_Spectrum,
        AsnMixin_NIRISS,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """NIRISS Single-Object Slitless"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'detector': {
                'value': 'NIS',
                'inputs': ['DETECTOR']
            },
            'exp_type': {
                'value': 'NIS_SOSS',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': 'GR700XD',
                'inputs': ['PUPIL']
            },
            'subarray': {
                'value': 'FULL|SUBSTRIP256|SUBSTRIP80',
                'inputs': ['SUBARRAY'],
                'force_unique': True
            }
        })

        # Check and continue initialization.
        super(Asn_NIR_SO_SLITLESS, self).__init__(*args, **kwargs)


class Asn_NIR_FIXEDSLIT(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """NIRSPEC Fixed Slit"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'NRS_FIXEDSLIT',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER']
            },
            'grating': {
                'value': None,
                'inputs': ['GRATING']
            },
            'fixed_slit': {
                'value': None,
                'inputs': ['FXD_SLIT']
            },
            'subarray': {
                'value': None,
                'inputs': ['SUBARRAY']
            },
        })

        # Check and continue initialization.
        super(Asn_NIR_FIXEDSLIT, self).__init__(*args, **kwargs)


class Asn_NIR_MSA(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """NIRSPEC MSA"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'pointing_type': {
                'value': 'SCIENCE',
                'inputs': ['PNTGTYPE']
            },
            'exp_type': {
                'value': 'NRS_MSASPEC',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER']
            },
            'grating': {
                'value': None,
                'inputs': ['GRATING']
            },
        })

        # Check and continue initialization.
        super(Asn_NIR_MSA, self).__init__(*args, **kwargs)


class Asn_MIRI_MRS(
        AsnMixin_Spectrum,
        AsnMixin_MIRI,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """MIRI MRS (IFU)"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'MIR_MRS',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['BAND']
            },
        })

        # Check and continue initialization.
        super(Asn_MIRI_MRS, self).__init__(*args, **kwargs)


class Asn_NIR_IFU(
        AsnMixin_Spectrum,
        AsnMixin_NIRSPEC,
        AsnMixin_Target,
        AsnMixin_Unique_Config
):
    """NIRSPEC IFU"""

    def __init__(self, *args, **kwargs):

        # Setup for checking.
        self.add_constraints({
            'exp_type': {
                'value': 'NRS_IFU',
                'inputs': ['EXP_TYPE']
            },
            'opt_elem': {
                'value': None,
                'inputs': ['FILTER']
            },
            'grating': {
                'value': None,
                'inputs': ['GRATING']
            }
        })

        # Check and continue initialization.
        super(Asn_NIR_IFU, self).__init__(*args, **kwargs)
