"""Association attributes common to DMS-based Rules"""
from jwst.associations.lib.counter import Counter

from jwst.associations.exceptions import (
    AssociationNotAConstraint,
    AssociationNotValidError,
)
from jwst.associations.lib.acid import ACIDMixin
from jwst.associations.lib.constraint import (Constraint, AttrConstraint, SimpleConstraint)
from jwst.associations.lib.utilities import getattr_from_list


__all__ = ['Constraint_TargetAcq', 'Constraint_TSO', 'Constraint_WFSC', 'DMSBaseMixin']

# Default product name
PRODUCT_NAME_DEFAULT = 'undefined'

# DMS file name templates
_ASN_NAME_TEMPLATE_STAMP = 'jw{program}-{acid}_{stamp}_{type}_{sequence:05d}_asn'
_ASN_NAME_TEMPLATE = 'jw{program}-{acid}_{type}_{sequence:05d}_asn'

# Acquisition and Confirmation images
ACQ_EXP_TYPES = (
    'mir_tacq',
    'mir_taconfirm',
    'nis_taconfirm',
    'nis_tacq',
    'nrc_taconfirm',
    'nrc_tacq',
    'nrs_confirm',
    'nrs_msata',
    'nrs_taconfirm',
    'nrs_tacq',
    'nrs_taslit',
    'nrs_verify',
    'nrs_wata',
)

# Exposure EXP_TYPE to Association EXPTYPE mapping
EXPTYPE_MAP = {
    'mir_darkall':       'dark',
    'mir_darkimg':       'dark',
    'mir_darkmrs':       'dark',
    'mir_flatimage':     'flat',
    'mir_flatmrs':       'flat',
    'mir_flatimage-ext': 'flat',
    'mir_flatmrs-ext':   'flat',
    'mir_tacq':          'target_acquisition',
    'mir_taconfirm':     'target_acquisition',
    'nis_dark':          'dark',
    'nis_focus':         'engineering',
    'nis_lamp':          'engineering',
    'nis_tacq':          'target_acquisition',
    'nis_taconfirm':     'target_acquisition',
    'nrc_dark':          'dark',
    'nrc_flat':          'flat',
    'nrc_focus':         'engineering',
    'nrc_led':           'engineering',
    'nrc_tacq':          'target_acquisition',
    'nrc_taconfirm':     'target_acquisition',
    'nrs_autoflat':      'autoflat',
    'nrs_autowave':      'autowave',
    'nrs_confirm':       'target_acquisition',
    'nrs_dark':          'dark',
    'nrs_focus':         'engineering',
    'nrs_image':         'engineering',
    'nrs_lamp':          'engineering',
    'nrs_msata':         'target_acquisition',
    'nrs_tacq':          'target_acquisition',
    'nrs_taconfirm':     'target_acquisition',
    'nrs_taslit':        'target_acquisition',
    'nrs_wata':          'target_acquisition',
}

# Coronographic exposures
CORON_EXP_TYPES = [
    'mir_4qpm',
    'mir_lyot',
    'nrc_coron'
]

# Exposures that get Level2b processing
IMAGE2_SCIENCE_EXP_TYPES = [
    'fgs_image',
    'mir_4qpm',
    'mir_image',
    'mir_lyot',
    'nis_ami',
    'nis_image',
    'nrc_coron',
    'nrc_image',
    'nrs_mimf',
    'nrc_tsimage',
]

IMAGE2_NONSCIENCE_EXP_TYPES = [
    'mir_coroncal',
    'nis_focus',
    'nrc_focus',
    'nrs_focus',
    'nrs_image',
]
IMAGE2_NONSCIENCE_EXP_TYPES.extend(ACQ_EXP_TYPES)

SPEC2_SCIENCE_EXP_TYPES = [
    'mir_lrs-fixedslit',
    'mir_lrs-slitless',
    'mir_mrs',
    'nis_soss',
    'nis_wfss',
    'nrc_tsgrism',
    'nrc_wfss',
    'nrs_fixedslit',
    'nrs_ifu',
    'nrs_msaspec',
    'nrs_brightobj',
]

SPECIAL_EXPOSURE_MODIFIERS = {
    'background': ['bkgdtarg'],
    'imprint': ['is_imprt'],
    'psf': ['is_psf'],
}

# Exposures that are always TSO
TSO_EXP_TYPES = [
    'nrc_tsimage',
    'nrc_tsgrism',
    'nrs_brightobj'
]

# Define the valid optical paths vs detector for NIRSpect Fixed-slit Science
# Tuples are (SLIT, GRATING, FILTER, DETECTOR)
# All A-slits are represented by SLIT == 'a'.
NRS_FSS_VALID_OPTICAL_PATHS = (
    ('a', 'prism', 'clear',  'nrs1'),
    ('a', 'g395h', 'f290lp', 'nrs1'),
    ('a', 'g395h', 'f290lp', 'nrs2'),
    ('a', 'g235h', 'f170lp', 'nrs1'),
    ('a', 'g235h', 'f170lp', 'nrs2'),
    ('a', 'g140h', 'f100lp', 'nrs1'),
    ('a', 'g140h', 'f100lp', 'nrs2'),
    ('a', 'g140h', 'f070lp', 'nrs1'),
    ('a', 'g395m', 'f290lp', 'nrs1'),
    ('a', 'g235m', 'f170lp', 'nrs1'),
    ('a', 'g140m', 'f100lp', 'nrs1'),
    ('a', 'g140m', 'f070lp', 'nrs1'),
    ('s200b1', 'prism', 'clear',  'nrs2'),
    ('s200b1', 'g395h', 'f290lp', 'nrs2'),
    ('s200b1', 'g235h', 'f170lp', 'nrs2'),
    ('s200b1', 'g140h', 'f100lp', 'nrs2'),
    ('s200b1', 'g140h', 'f070lp', 'nrs1'),
    ('s200b1', 'g140h', 'f070lp', 'nrs2'),
    ('s200b1', 'g395m', 'f290lp', 'nrs2'),
    ('s200b1', 'g235m', 'f170lp', 'nrs2'),
    ('s200b1', 'g140m', 'f100lp', 'nrs2'),
    ('s200b1', 'g140m', 'f070lp', 'nrs1'),
    ('s200b1', 'g140m', 'f070lp', 'nrs2'),
)

NRS_FSS_VALID_LAMP_OPTICAL_PATHS = (
    ('prism', 'line4', 'nrs1'),
    ('prism', 'line4', 'nrs2'),
    ('prism', 'flat5', 'nrs1'),
    ('prism', 'flat5', 'nrs2'),
    ('prism', 'test', 'nrs1'),
    ('prism', 'test', 'nrs2'),
    ('g395h', 'flat3', 'nrs1'),
    ('g395h', 'flat3', 'nrs2'),
    ('g395h', 'line3', 'nrs1'),
    ('g395h', 'line3', 'nrs2'),
    ('g235h', 'flat2', 'nrs1'),
    ('g235h', 'flat2', 'nrs2'),
    ('g235h', 'line2', 'nrs1'),
    ('g235h', 'line2', 'nrs2'),
    ('g140h', 'flat1', 'nrs1'),
    ('g140h', 'flat1', 'nrs2'),
    ('g140h', 'ref', 'nrs1'),
    ('g140h', 'ref', 'nrs2'),
    ('g140h', 'line1', 'nrs1'),
    ('g140h', 'line1', 'nrs2'),
    ('g140h', 'flat4', 'nrs1'),
    ('g140h', 'flat4', 'nrs2'),
    ('g395m', 'flat3', 'nrs1'),
    ('g395m', 'flat3', 'nrs2'),
    ('g395m', 'line3', 'nrs1'),
    ('g395m', 'line3', 'nrs2'),
    ('g235m', 'flat2', 'nrs1'),
    ('g235m', 'flat2', 'nrs2'),
    ('g235m', 'line2', 'nrs1'),
    ('g235m', 'line2', 'nrs2'),
    ('g140m', 'flat1', 'nrs1'),
    ('g140m', 'flat1', 'nrs2'),
    ('g140m', 'line1', 'nrs1'),
    ('g140m', 'line1', 'nrs2'),
    ('g140m', 'flat4', 'nrs1'),
    ('g140m', 'flat4', 'nrs2'),
)

# Key that uniquely identifies members.
MEMBER_KEY = 'expname'

# Non-specified values found in DMS Association Pools
_EMPTY = (None, '', 'NULL', 'Null', 'null', '--', 'N', 'n',
          'F', 'f', 'FALSE', 'false', 'False', 'N/A', 'n/a')

# Degraded status information
_DEGRADED_STATUS_OK = (
    'No known degraded exposures in association.'
)
_DEGRADED_STATUS_NOTOK = (
    'One or more members have an error associated with them.'
    '\nDetails can be found in the member.exposerr attribute.'
)


class DMSBaseMixin(ACIDMixin):
    """Association attributes common to DMS-based Rules

    Attributes
    ----------
    sequence : int
        The sequence number of the current association
    """

    # Associations of the same type are sequenced.
    _sequence = Counter(start=1)

    def __init__(self, *args, **kwargs):
        super(DMSBaseMixin, self).__init__(*args, **kwargs)

        self._acid = None
        self._asn_name = None
        self.sequence = None
        if 'degraded_status' not in self.data:
            self.data['degraded_status'] = _DEGRADED_STATUS_OK
        if 'program' not in self.data:
            self.data['program'] = 'noprogram'

    @classmethod
    def create(cls, item, version_id=None):
        """Create association if item belongs

        Parameters
        ----------
        item : dict
            The item to initialize the association with.

        version_id : str or None
            Version_Id to use in the name of this association.
            If None, nothing is added.

        Returns
        -------
        (association, reprocess_list)
            2-tuple consisting of:

                - association : The association or, if the item does not
                  match this rule, None
                - [ProcessList[, ...]]: List of items to process again.
        """
        asn, reprocess = super(DMSBaseMixin, cls).create(item, version_id)
        if not asn:
            return None, reprocess
        asn.sequence = next(asn._sequence)
        return asn, reprocess

    @property
    def acid(self):
        """Association ID"""
        acid = self._acid
        if self._acid is None:
            acid = self.acid_from_constraints()
        return acid

    @property
    def asn_name(self):
        """The association name

        The name that identifies this association. When dumped,
        will form the basis for the suggested file name.

        Typically, it is generated based on the current state of
        the association, but can be overridden.
        """
        if self._asn_name:
            return self._asn_name

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

    @asn_name.setter
    def asn_name(self, name):
        """Override calculated association name"""
        self._asn_name = name

    @property
    def current_product(self):
        return self.data['products'][-1]

    @property
    def from_items(self):
        """The list of items that contributed to the association."""
        try:
            items = [
                member.item
                for product in self['products']
                for member in product['members']
            ]
        except KeyError:
            items = []
        return items

    @property
    def member_ids(self):
        """Set of all member ids in all products of this association"""
        member_ids = set(
            member[MEMBER_KEY]
            for product in self['products']
            for member in product['members']
        )
        return member_ids

    @property
    def validity(self):
        """Keeper of the validity tests"""
        try:
            validity = self._validity
        except AttributeError:
            self._validity = {}
            validity = self._validity
        return validity

    @validity.setter
    def validity(self, item):
        """Set validity dict"""
        self._validity = item

    def get_exposure_type(self, item, default='science'):
        """Determine the exposure type of a pool item

        Parameters
        ----------
        item : dict
            The pool entry to determine the exposure type of

        default : str or None
            The default exposure type.
            If None, routine will raise LookupError

        Returns
        -------
        exposure_type : str
            Exposure type. Can be one of

                - 'science': Item contains science data
                - 'target_acquisition': Item contains target acquisition data.
                - 'autoflat': NIRSpec AUTOFLAT
                - 'autowave': NIRSpec AUTOWAVE
                - 'psf': PSF
                - 'imprint': MSA/IFU Imprint/Leakcal

        Raises
        ------
        LookupError
            When `default` is None and an exposure type cannot be determined
        """
        return get_exposure_type(item, default=default, association=self)

    def is_member(self, new_member):
        """Check if member is already a member

        Parameters
        ----------
        new_member : Member
            The member to check for
        """
        try:
            current_members = self.current_product['members']
        except KeyError:
            return False

        for member in current_members:
            if member == new_member:
                return True
        return False

    def is_item_member(self, item):
        """Check if item is already a member of this association

        Parameters
        ----------
        item : dict
            The item to check for.

        Returns
        -------
        is_item_member : bool
            True if item is a member.
        """
        return item in self.from_items

    def is_item_tso(self, item, other_exp_types=None):
        """Is the given item TSO

        Determine whether the specific item represents
        TSO data or not. When used to determine naming
        of files, coronagraphic data will be included through
        the `other_exp_types` parameter.

        Parameters
        ----------
        item : dict
            The item to check for.

        other_exp_types: [str[,...]] or None
            List of other exposure types to consider TSO.

        Returns
        -------
        is_item_tso : bool
            Item represents a TSO exposure.
        """
        # If not a science exposure, such as target acquisitions,
        # then other TSO indicators do not apply.
        if item['pntgtype'] != 'science':
            return False

        # Target acquisitions are never TSO
        if item['exp_type'] in ACQ_EXP_TYPES:
            return False

        # Setup exposure list
        all_exp_types = TSO_EXP_TYPES.copy()
        if other_exp_types:
            all_exp_types += other_exp_types

        # Go through all other TSO indicators.
        try:
            is_tso = self.constraints['is_tso'].value == 't'
        except (AttributeError, KeyError):
            # No such constraint is defined. Just continue on.
            is_tso = False
        try:
            is_tso = is_tso or self.item_getattr(item, ['tsovisit'])[1] == 't'
        except KeyError:
            pass
        try:
            is_tso = is_tso or self.item_getattr(item, ['exp_type'])[1] in all_exp_types
        except KeyError:
            pass
        return is_tso

    def item_getattr(self, item, attributes):
        """Return value from any of a list of attributes

        Parameters
        ----------
        item : dict
            item to retrieve from

        attributes : list
            List of attributes

        Returns
        -------
        (attribute, value)
            Returns the value and the attribute from
            which the value was taken.

        Raises
        ------
        KeyError
            None of the attributes are found in the dict.
        """
        return item_getattr(item, attributes, self)

    def new_product(self, product_name=PRODUCT_NAME_DEFAULT):
        """Start a new product"""
        product = {
            'name': product_name,
            'members': []
        }
        try:
            self.data['products'].append(product)
        except (AttributeError, KeyError):
            self.data['products'] = [product]

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
        self.update_degraded_status()

    def update_degraded_status(self):
        """Update association degraded status"""

        if self.data['degraded_status'] == _DEGRADED_STATUS_OK:
            for product in self.data['products']:
                for member in product['members']:
                    try:
                        exposerr = member['exposerr']
                    except KeyError:
                        continue
                    else:
                        if exposerr not in _EMPTY:
                            self.data['degraded_status'] = _DEGRADED_STATUS_NOTOK
                            break

    def update_validity(self, entry):
        for test in self.validity.values():
            if not test['validated']:
                test['validated'] = test['check'](entry)

    @classmethod
    def reset_sequence(cls):
        cls._sequence = Counter(start=1)

    @classmethod
    def validate(cls, asn):
        super(DMSBaseMixin, cls).validate(asn)

        if isinstance(asn, DMSBaseMixin):
            result = False
            try:
                result = all(
                    test['validated']
                    for test in asn.validity.values()
                )
            except (AttributeError, KeyError):
                raise AssociationNotValidError('Validation failed')
            if not result:
                raise AssociationNotValidError(
                    'Validation failed validity tests.'
                )

        return True

    def _get_exposure(self):
        """Get string representation of the exposure id

        Returns
        -------
        exposure : str
            The Level3 Product name representation
            of the exposure & activity id.
        """
        exposure = ''
        try:
            activity_id = format_list(
                self.constraints['activity_id'].found_values
            )
        except KeyError:
            pass
        else:
            if activity_id not in _EMPTY:
                exposure = '{0:0>2s}'.format(activity_id)
        return exposure

    def _get_instrument(self):
        """Get string representation of the instrument

        Returns
        -------
        instrument : str
            The Level3 Product name representation
            of the instrument
        """
        instrument = format_list(self.constraints['instrument'].found_values)
        return instrument

    def _get_opt_element(self):
        """Get string representation of the optical elements

        Returns
        -------
        opt_elem : str
            The Level3 Product name representation
            of the optical elements.
        """
        # Retrieve all the optical elements
        opt_elems = []
        for opt_elem in ['opt_elem', 'opt_elem2', 'opt_elem3']:
            try:
                values = list(self.constraints[opt_elem].found_values)
            except KeyError:
                pass
            else:
                values.sort(key=str.lower)
                value = format_list(values)
                if value not in _EMPTY:
                    opt_elems.append(value)

        # Build the string. Sort the elements in order to
        # create data-independent results
        opt_elems.sort(key=str.lower)
        opt_elem = '-'.join(opt_elems)
        if opt_elem == '':
            opt_elem = 'clear'

        return opt_elem

    def _get_subarray(self):
        """Get string representation of the subarray

        Returns
        -------
        subarray : str
            The Level3 Product name representation
            of the subarray.
        """
        result = ''
        try:
            subarray = format_list(self.constraints['subarray'].found_values)
        except KeyError:
            subarray = None
        if subarray == 'full':
            subarray = None
        if subarray is not None:
            result = subarray

        return result

    def _get_target(self):
        """Get string representation of the target

        Returns
        -------
        target : str
            The Level3 Product name representation
            of the target or source ID.
        """
        target_id = format_list(self.constraints['target'].found_values)
        target = 't{0:0>3s}'.format(str(target_id))
        return target

    def _get_grating(self):
        """Get string representation of the grating in use

        Returns
        -------
        grating : str
            The Level3 Product name representation
            of the grating in use.
        """
        grating_id = format_list(self.constraints['grating'].found_values)
        grating = '{0:0>3s}'.format(str(grating_id))
        return grating


# -----------------
# Basic constraints
# -----------------
class DMSAttrConstraint(AttrConstraint):
    """DMS-focused attribute constraint

    Forces definition of invalid values
    """
    def __init__(self, **kwargs):

        if kwargs.get('invalid_values', None) is None:
            kwargs['invalid_values'] = _EMPTY

        super(DMSAttrConstraint, self).__init__(**kwargs)


class Constraint_TargetAcq(SimpleConstraint):
    """Select on target acquisition exposures

    Parameters
    ----------
    association:  ~jwst.associations.Association
        If specified, use the `get_exposure_type` method
        of the association rather than the utility version.
    """
    def __init__(self, association=None):
        if association is None:
            _get_exposure_type = get_exposure_type
        else:
            _get_exposure_type = association.get_exposure_type

        super(Constraint_TargetAcq, self).__init__(
            name='target_acq',
            value='target_acquisition',
            sources=_get_exposure_type
        )


class Constraint_TSO(Constraint):
    """Match on Time-Series Observations"""
    def __init__(self, *args, **kwargs):
        super(Constraint_TSO, self).__init__(
            [
                DMSAttrConstraint(
                    sources=['pntgtype'],
                    value='science'
                ),
                Constraint(
                    [
                        DMSAttrConstraint(
                            sources=['tsovisit'],
                            value='t',
                        ),
                        DMSAttrConstraint(
                            sources=['exp_type'],
                            value='|'.join(TSO_EXP_TYPES),
                        ),
                    ],
                    reduce=Constraint.any
                )
            ],
            name='is_tso'
        )


class Constraint_WFSC(Constraint):
    """Match on Wave Front Sensing and Control Observations"""
    def __init__(self, *args, **kwargs):
        super(Constraint_WFSC, self).__init__(
            [
                Constraint(
                    [
                        DMSAttrConstraint(
                            name='wfsc',
                            sources=['visitype'],
                            value='.+wfsc.+',
                            force_unique=True
                        )
                    ]
                )
            ]
        )


# #########
# Utilities
# #########
def format_list(alist):
    """Format a list according to DMS naming specs"""
    return '-'.join(alist)


def get_exposure_type(item, default='science', association=None):
    """Determine the exposure type of a pool item

    Parameters
    ----------
    item : dict
        The pool entry to determine the exposure type of

    default : str or None
        The default exposure type.
        If None, routine will raise LookupError



    Returns
    -------
    exposure_type : str
        Exposure type. Can be one of

        - 'science': Item contains science data
        - 'target_acquisition': Item contains target acquisition data.
        - 'autoflat': NIRSpec AUTOFLAT
        - 'autowave': NIRSpec AUTOWAVE
        - 'psf': PSF
        - 'imprint': MSA/IFU Imprint/Leakcal

    Raises
    ------
    LookupError
        When `default` is None and an exposure type cannot be determined
    """
    # Specify how attributes of the item are retrieved.
    def _item_attr(item, sources):
        """Get attribute value of an item

        This simplifies the call to `item_getattr`
        """
        source, value = item_getattr(item, sources, association=association)
        return value

    # Define default type.
    result = default

    # Retrieve pointing type. This decides the basic exposure type.
    # If the pointing is not science, we're done.
    try:
        result = _item_attr(item, ['pntgtype'])
    except KeyError:
        pass
    else:
        if result != 'science':
            return result

    # We have a science exposure. Refine further.
    #
    # Base type off of exposure type.
    try:
        exp_type = _item_attr(item, ['exp_type'])
    except KeyError:
        raise LookupError('Exposure type cannot be determined')

    result = EXPTYPE_MAP.get(exp_type, default)

    if result is None:
        raise LookupError('Cannot determine exposure type')

    # If result is not science, we're done.
    if result != 'science':
        return result

    # For `science` data, compare against special modifiers
    # to further refine the type.
    for special, source in SPECIAL_EXPOSURE_MODIFIERS.items():
        try:
            _item_attr(item, source)
        except KeyError:
            pass
        else:
            result = special
            break

    return result


def item_getattr(item, attributes, association=None):
    """Return value from any of a list of attributes

    Parameters
    ----------
    item : dict
        item to retrieve from

    attributes : list
        List of attributes

    Returns
    -------
    (attribute, value)
        Returns the value and the attribute from
        which the value was taken.

    Raises
    ------
    KeyError
        None of the attributes are found in the dict.
    """
    if association is None:
        invalid_values = _EMPTY
    else:
        invalid_values = association.INVALID_VALUES
    return getattr_from_list(
        item,
        attributes,
        invalid_values=invalid_values
    )


def nrsfss_valid_detector(item):
    """Check that a grating/filter combo can appear on the detector"""
    try:
        _, detector = item_getattr(item, ['detector'])
        _, filter = item_getattr(item, ['filter'])
        _, grating = item_getattr(item, ['grating'])
        _, slit = item_getattr(item, ['fxd_slit'])
    except KeyError:
        return False

    # Reduce all A slits to just 'a'.
    if slit != 's200b1':
        slit = 'a'

    return (slit, grating, filter, detector) in NRS_FSS_VALID_OPTICAL_PATHS


def nrsifu_valid_detector(item):
    """Check that a grating/filter combo can appear on the detector"""
    try:
        _, detector = item_getattr(item, ['detector'])
        _, filter = item_getattr(item, ['filter'])
        _, grating = item_getattr(item, ['grating'])
    except KeyError:
        return False

    # Just a checklist of paths:
    if grating in ['g395h', 'g235h']:
        return True
    elif grating in ['g395m', 'g235m', 'g140m'] and detector == 'nrs1':
        return True
    elif grating == 'prism' and filter == 'clear' and detector == 'nrs1':
        return True
    elif grating == 'g140h':
        if filter == 'f100lp':
            return True
        elif filter == 'f070lp' and detector == 'nrs1':
            return True

    # Nothing has matched. Not valid.
    return False


def nrslamp_valid_detector(item):
    """Check that a grating/lamp combo can appear on the detector"""
    try:
        _, detector = item_getattr(item, ['detector'])
        _, grating = item_getattr(item, ['grating'])
        _, opmode = item_getattr(item, ['opmode'])
    except KeyError:
        return False

    if opmode in ['msaspec']:
        # All settings can result in data on both detectors,
        # depending on which MSA shutters are open
        return True

    elif opmode in ['ifu']:
        # Need lamp value for IFU mode
        try:
            _, lamp = item_getattr(item, ['lamp'])
        except KeyError:
            return False

        # Just a checklist of paths:
        if grating in ['g395h', 'g235h']:
            # long-wave, high-res gratings result in data on both detectors
            return True
        elif grating in ['g395m', 'g235m', 'g140m'] and detector == 'nrs1':
            # all medium-res gratings result in data only on NRS1 detector
            return True
        elif grating == 'prism' and detector == 'nrs1':
            # prism results in data only on NRS1 detector
            return True
        elif grating == 'g140h':
            # short-wave, high-res grating results in data on both detectors,
            # except when lamp FLAT4 is in use (no data on NRS2)
            if not (detector == 'nrs2' and lamp == 'flat4'):
                return True

    elif opmode in ['fixedslit']:
        # All slits illuminated by lamps, regardless of grating or subarray
        try:
            _, lamp = item_getattr(item, ['lamp'])
        except KeyError:
            return False

        return (grating, lamp, detector) in NRS_FSS_VALID_LAMP_OPTICAL_PATHS

    # Nothing has matched. Not valid.
    return False
