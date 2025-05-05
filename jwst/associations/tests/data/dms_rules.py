"""Avoid relative imports and mimic actual usage."""  # noqa: INP001

from jwst.associations import Association
from jwst.associations.association import finalize as general_asn_finalize
from jwst.associations.registry import RegistryMarker
from jwst.associations.lib.constraint import ConstraintTrue
from jwst.associations.lib.dms_base import DMSBaseMixin


@RegistryMarker.rule
class AsnDMSBase(DMSBaseMixin, Association):
    """Basic DMS rule."""

    def __init__(self, version_id=None):
        self.constraints = ConstraintTrue()
        super(AsnDMSBase, self).__init__(version_id=version_id)
        self.data["members"] = []

    def make_member(self, item):
        """
        Return item; placeholder method.

        Parameters
        ----------
        item : Member
            Member.

        Returns
        -------
        Member
            The item provided - this is a placeholder method.
        """
        return item

    def _add(self, item):
        self.data["members"].append(item)

    def finalize(self):
        """
        Perform finalization steps.

        Returns
        -------
        list(self)
            Return self in a list.
        """
        return [self]


# Use the generic finalization
RegistryMarker.callback("finalize")(general_asn_finalize)


class Utility:
    """Utility class that should not be part of the base utilities."""

    @staticmethod
    def not_valid_function():
        """Mark as not valid, a function that should not be part of the utilities."""


@RegistryMarker.utility
class ValidUtility:
    """Yes, valid."""

    @staticmethod
    def valid_function():
        """Yes, I'm good."""
