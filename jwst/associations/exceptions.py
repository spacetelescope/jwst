__all__ = [
    'AssociationError',
    'AssociationNotAConstraint',
    'AssociationNotValidError',
    'AssociationProcessMembers',
]


class AssociationError(Exception):
    """Basic errors related to Associations"""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return 'Association Exception: {}'.format(self.message)


class AssociationNotAConstraint(AssertionError):
    """No matching constraint found"""


class AssociationNotValidError(AssociationError):
    """Given data structure is not a valid association"""


class AssociationProcessMembers(Exception):
    """Add members for processing

    Attributes
    ----------
    members: AssociationPool or [member, ...]
        List of members to reprocess.

    allowed_rules: [type(Association), ...]
        List of allowed rules to match members against.

    message: str
        Text explaining the exception.
    """
    def __init__(self, members, allowed_rules):
        self.members = members
        self.allowed_rules = allowed_rules
        self.message = '{} members to be reprocessed for rules {}'.format(
            len(members),
            allowed_rules
        )
