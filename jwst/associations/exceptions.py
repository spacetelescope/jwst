__all__ = [
    'AssociationError',
    'AssociationNotAConstraint',
    'AssociationNotValidError',
]


class AssociationError(Exception):
    """Basic errors related to Associations"""


class AssociationNotAConstraint(AssociationError):
    """No matching constraint found"""


class AssociationNotValidError(AssociationError):
    """Given data structure is not a valid association"""
