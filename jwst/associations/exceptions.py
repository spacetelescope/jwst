__all__ = [
    'AssociationError',
    'AssociationNotAConstraint',
    'AssociationNotValidError',
]


class AssociationError(Exception):
    """Basic errors related to Associations"""

    def __init__(self, message='No explanation given'):
        self.message = message

    def __str__(self):
        return 'Association Exception: {}'.format(self.message)


class AssociationNotAConstraint(AssociationError):
    """No matching constraint found"""


class AssociationNotValidError(AssociationError):
    """Given data structure is not a valid association"""
