"""Association exceptions"""


class AssociationError(Exception):
    """Basic failure of an association"""


class AssociationNotValidError(AssociationError):
    """Given data structure is not a valid association"""
