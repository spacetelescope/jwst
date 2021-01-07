"""
Define the I/O methods for Level 3 associations
"""
import json as json_lib
import logging
import numpy as np
import yaml as yaml_lib

from .association import Association
from .exceptions import AssociationNotValidError
from .lib.member import Member

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = []


# Define JSON encoder to convert `Member` to `dict`
class AssociationEncoder(json_lib.JSONEncoder):
    """Encode to handle Associations"""
    def default(self, obj):

        # Convert Member to a simple dict
        if isinstance(obj, Member):
            return obj.data


@Association.ioregistry
class json():
    """Load and store associations as JSON"""

    @staticmethod
    def load(cls, serialized):
        """Unserialize an association from JSON

        Parameters
        ----------
        cls : class
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
            logger.debug('Error unserializing: "{}"'.format(err))
            raise AssociationNotValidError(
                'Containter is not JSON: "{}"'.format(serialized)
            )

        return asn

    @staticmethod
    def dump(asn):
        """Create JSON representation.

        Parameters
        ----------
        asn : Association
            The association to serialize

        Returns
        -------
        (name, str):
            Tuple where the first item is the suggested
            Name for the JSON file.
            Second item is the string containing the JSON serialization.
        """
        asn_filename = asn.asn_name
        if not asn_filename.endswith('.json'):
            asn_filename = asn_filename+'.json'
        return (
            asn_filename,
            json_lib.dumps(asn.data, cls=AssociationEncoder, indent=4, separators=(',', ': '))
        )


@Association.ioregistry
class yaml():
    """Load and store associations as YAML"""

    @staticmethod
    def load(cls, serialized):
        """Unserialize an association from YAML

        Parameters
        ----------
        cls : class
            The class from which further information will be gathered
            and possibly instantiated.

        serialized : str or file object
            The YAML to read

        Returns
        -------
        association : dict
            The association

        Raises
        ------
        AssociationNotValidError
            Cannot create or validate the association.
        """
        try:
            serialized.seek(0)
        except AttributeError:
            pass
        try:
            asn = yaml_lib.safe_load(serialized)
        except Exception as err:
            logger.debug('Error unserializing: "{}"'.format(err))
            raise AssociationNotValidError(
                'Container is not YAML: "{}"'.format(serialized)
            )
        return asn

    @staticmethod
    def dump(asn):
        """Create YAML representation.

         Parameters
        ----------
        asn : Association
            The association to serialize


        Returns
        -------
        (name, str):
            Tuple where the first item is the suggested
            Name for the YAML file.
            Second item is the string containing the YAML serialization.
        """
        asn_filename = asn.asn_name
        if not asn.asn_name.endswith('.yaml'):
            asn_filename = asn.asn_name+'.yaml'
        return (
            asn_filename,
            yaml_lib.dump(asn.data, default_flow_style=False)
        )


# Register YAML representers
def np_str_representer(dumper, data):
    """Convert numpy.str_ into standard YAML string"""
    return dumper.represent_scalar('tag:yaml.org,2002:str', str(data))
yaml_lib.add_representer(np.str_, np_str_representer)


def member_representer(dumper, member):
    """Convert a Member to its basic dict representation"""
    return dumper.represent_dict(member.data)
yaml_lib.add_representer(Member, member_representer)
