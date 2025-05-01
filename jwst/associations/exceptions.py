__all__ = [
    "AssociationError",
    "AssociationNotAConstraintError",
    "AssociationNotValidError",
]


class AssociationError(Exception):
    """Basic errors related to Associations."""


class AssociationNotAConstraintError(AssociationError):
    """No matching constraint found."""


class AssociationNotValidError(AssociationError):
    """Given data structure is not a valid association."""
