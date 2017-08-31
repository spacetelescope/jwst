"""Association attributes common to DMS-based Rules"""
from .counter import Counter

from jwst.associations.association import getattr_from_list
from jwst.associations.exceptions import (
    AssociationNotValidError,
)
from jwst.associations.lib.acid import ACIDMixin

# Default product name
PRODUCT_NAME_DEFAULT = 'undefined'

# DMS file name templates
_ASN_NAME_TEMPLATE_STAMP = 'jw{program}-{acid}_{stamp}_{type}_{sequence:03d}_asn'
_ASN_NAME_TEMPLATE = 'jw{program}-{acid}_{type}_{sequence:03d}_asn'

# Exposure EXP_TYPE to Association EXPTYPE mapping
_EXPTYPE_MAP = {
    'mir_tacq':      'target_acquistion',
    'nis_tacq':      'target_acquistion',
    'nis_taconfirm': 'target_acquistion',
    'nrc_tacq':      'target_acquistion',
    'nrc_taconfirm': 'target_acquistion',
    'nrs_autoflat':  'autoflat',
    'nrs_autowave':  'autowave',
    'nrs_confirm':   'target_acquistion',
    'nrs_tacq':      'target_acquistion',
    'nrs_taconfirm': 'target_acquistion',
    'nrs_taslit':    'target_acquistion',
}

__all__ = ['DMSBaseMixin']


class DMSBaseMixin(ACIDMixin):
    """Association attributes common to DMS-based Rules"""

    # Associations of the same type are sequenced.
    _sequence = Counter(start=1)

    def __init__(self, *args, **kwargs):
        super(DMSBaseMixin, self).__init__(*args, **kwargs)

        self.sequence = None
        self.data.update({
            'program': 'none',
        })

    @classmethod
    def create(cls, member, version_id=None):
        """Create association if member belongs

        Parameters
        ----------
        member: dict
            The member to initialize the association with.

        version_id: str or None
            Version_Id to use in the name of this association.
            If None, nothing is added.

        Returns
        -------
        (association, reprocess_list)
            2-tuple consisting of:
            - association: The association or, if the member does not
                this rule, None
            - [ProcessList[, ...]]: List of members to process again.
        """
        asn, reprocess = super(DMSBaseMixin, cls).create(member, version_id)
        if not asn:
            return None, reprocess
        asn.sequence = next(asn._sequence)
        return asn, reprocess

    @property
    def acid(self):
        """Association ID"""
        return self.acid_from_constraints()

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

    @property
    def current_product(self):
        return self.data['products'][-1]

    def new_product(self, product_name=PRODUCT_NAME_DEFAULT):
        """Start a new product"""
        product = {
            'name': product_name,
            'members': []
        }
        try:
            self.data['products'].append(product)
        except KeyError:
            self.data['products'] = [product]

    def update_validity(self, entry):
        for test in self.validity.values():
            if not test['validated']:
                test['validated'] = test['check'](entry)

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

    @classmethod
    def reset_sequence(cls):
        cls._sequence = Counter(start=1)

    def get_exposure_type(self, member, default='science'):
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
                'PSF': PSF

        Raises
        ------
        LookupError
            When `default` is None and an exposure type cannot be determined
        """
        result = default

        # Look for specific attributes
        try:
            self.member_getattr(member, ['is_psf'])
        except KeyError:
            pass
        else:
            return 'psf'

        # Base type off of exposure type.
        try:
            exp_type = member['exp_type']
        except KeyError:
            raise LookupError('Exposure type cannot be determined')

        result = _EXPTYPE_MAP.get(exp_type, default)

        if result is None:
            raise LookupError('Cannot determine exposure type')
        return result

    def member_getattr(self, member, attributes):
        """Return value from any of a list of attributes

        Parameters
        ----------
        member: dict
            member to retrieve from

        attributes: list
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
        return getattr_from_list(
            member,
            attributes,
            invalid_values=self.INVALID_VALUES
        )
