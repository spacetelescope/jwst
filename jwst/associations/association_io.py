"""
Define the I/O methods for Level 3 associations.

Particularly, load and store associations as JSON.
"""

import json as json_lib
import logging

from jwst.associations.exceptions import AssociationNotValidError
from jwst.associations.lib.member import Member

# Configure logging
logger = logging.getLogger(__name__)

__all__: list = []


class AssociationEncoder(json_lib.JSONEncoder):
    """JSON encoder to handle Associations and convert Member to dict."""

    def default(self, obj):
        """
        Convert Member to a simple dict.

        Parameters
        ----------
        obj : object
            The object - if a Member object, return its data attribute.

        Returns
        -------
        dict or None
            Return the Member data dictionary attribute, or
            None if obj is not a Member instance.
        """
        if isinstance(obj, Member):
            return obj.data


def json_asn_load(_cls, serialized):
    """
    Unserialize an association from JSON.

    Parameters
    ----------
    _cls : class
        The class from which further information will be gathered
        and possibly instantiated.

    serialized : str or file object
        The JSON to read

    Returns
    -------
    association : dict
        The association

    Raises
    ------
    AssociationNotValidError
        Cannot create or validate the association.
    """
    if isinstance(serialized, str):
        loader = json_lib.loads
    else:
        # Presume a file object
        serialized.seek(0)
        loader = json_lib.load
    try:
        asn = loader(serialized)
    except Exception as err:
        logger.debug(f'Error unserializing: "{err}"')
        raise AssociationNotValidError(f"Container is not JSON: '{serialized}'") from err

    return asn


def json_asn_dump(asn):
    """
    Create JSON representation.

    Parameters
    ----------
    asn : Association
        The association to serialize

    Returns
    -------
    name, str : tuple
        Tuple where the first item is the suggested
        Name for the JSON file.
        Second item is the string containing the JSON serialization.
    """
    asn_filename = asn.asn_name
    if not asn_filename.endswith(".json"):
        asn_filename = asn_filename + ".json"
    return (
        asn_filename,
        json_lib.dumps(asn.data, cls=AssociationEncoder, indent=4, separators=(",", ": ")),
    )
