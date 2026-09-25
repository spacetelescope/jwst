"""
Define the I/O methods for Level 3 associations.

Particularly, load and store associations as JSON.
"""

import json as json_lib
import logging
from pathlib import Path

from jwst.associations.exceptions import AssociationNotValidError
from jwst.associations.lib.member import Member

# Configure logging
logger = logging.getLogger(__name__)

__all__ = ["AssociationEncoder", "json_asn_load", "json_asn_dump"]


class AssociationEncoder(json_lib.JSONEncoder):
    """JSON encoder to handle Associations and convert Member to dict."""

    def default(self, obj):
        """
        Convert Member to a simple dict.

        Parameters
        ----------
        obj : `~jwst.associations.lib.member.Member`
            If input is a Member object, return its data attribute.

        Returns
        -------
        dict or None
            Return the `~jwst.associations.lib.member.Member`
            data attribute, otherwise None.
        """
        if isinstance(obj, Member):
            return obj.data


def json_asn_load(serialized):
    """
    Unserialize an association from JSON.

    Parameters
    ----------
    serialized : str, dict, or file-like
        The JSON to read.

    Returns
    -------
    association : dict
        The association data.

    Raises
    ------
    jwst.associations.exceptions.AssociationNotValidError
        Cannot create or validate the association.
    """
    if isinstance(serialized, dict):  # No-op
        return serialized

    if isinstance(serialized, (str, Path)):
        loader = json_lib.loads
    else:
        # Presume a file object
        serialized.seek(0)
        loader = json_lib.load
    try:
        asn = loader(serialized)
    except Exception as err:
        logger.debug('Error unserializing: "%s"', repr(err))
        raise AssociationNotValidError(f"Container is not JSON: '{serialized}'") from err

    return asn


def json_asn_dump(asn):
    """
    Create JSON representation.

    Parameters
    ----------
    asn : `~jwst.associations.association.Association`
        The association to serialize.

    Returns
    -------
    asn_filename : str
        Suggested name for the JSON file.
        This is taken from ``asn_name`` attribute of
        the given association.

    serialized : str
        JSON serialization of the given association.
    """
    asn_filename = asn.asn_name
    if not asn_filename.endswith(".json"):
        asn_filename = asn_filename + ".json"
    serialized = json_lib.dumps(asn.data, cls=AssociationEncoder, indent=4, separators=(",", ": "))
    return asn_filename, serialized
