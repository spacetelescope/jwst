class AssociationError(Exception):
    """Basic errors related to Associations"""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return 'Association Exception: {}'.format(self.message)
