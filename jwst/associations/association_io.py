"""
Define the I/O methods for Level 3 associations
"""
import json as json_lib
import jsonschema
import logging
import numpy as np
import six
import yaml as yaml_lib

from .association import Association
from .exceptions import AssociationError

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@Association.ioregistry
class json(object):

    @staticmethod
    def load(cls, serialized):
        """Unserialize an assocation from JSON

        Parameters
        ----------
        cls: class
            The class from which further information will be gathered
            and possibly instantiated.

        serialized: str or file object
            The JSON to read

        Returns
        -------
        association: dict
            The association
        """
        if isinstance(serialized, six.string_types):
            loader = json_lib.loads
        else:
            loader = json_lib.load
        try:
            asn = loader(serialized)
        except Exception as err:
            logger.debug('Error unserializing: "{}"'.format(err))
            raise AssociationError(
                'Containter is not JSON: "{}"'.format(serialized)
            )

        return asn

    @staticmethod
    def dump(asn):
        """Create JSON representation.

        Parameters
        ----------
        asn: Association
            The association to serialize

        Returns
        -------
        (name, str):
            Tuple where the first item is the suggested
            base name for the JSON file.
            Second item is the string containing the JSON serialization.
        """
        return (
            asn.asn_name,
            json_lib.dumps(asn.data, indent=4, separators=(',', ': '))
        )


@Association.ioregistry
class yaml(object):
    """Load and store associations as YAML"""

    @staticmethod
    def load(cls, serialized):
        """Unserialize an assocation from YAML

        Parameters
        ----------
        cls: class
            The class from which further information will be gathered
            and possibly instantiated.

        serialized: str or file object
            The YAML to read

        Returns
        -------
        association: dict
            The association
        """
        try:
            asn = yaml_lib.load(serialized)
        except Exception as err:
            logger.debug('Error unserializing: "{}"'.format(err))
            raise AssociationError(
                'Container is not YAML: "{}"'.format(serialized)
            )
        return asn

    @staticmethod
    def dump(asn):
        """Create YAML representation.

         Parameters
        ----------
        asn: Association
            The association to serialize


        Returns
        -------
        (name, str):
            Tuple where the first item is the suggested
            base name for the YAML file.
            Second item is the string containing the YAML serialization.
        """
        return (
            asn.asn_name,
            yaml_lib.dump(asn.data, default_flow_style=False)
        )


# Register YAML representers
def np_str_representer(dumper, data):
    """Convert numpy.str_ into standard YAML string"""
    return dumper.represent_scalar('tag:yaml.org,2002:str', str(data))
yaml_lib.add_representer(np.str_, np_str_representer)
