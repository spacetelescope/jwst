"""Open a Level2 association given various inputs"""

from os.path import splitext

from ..associations import (
    AssociationRegistry,
    libpath,
    load_asn
)
from ..associations.asn_from_list import asn_from_list
from ..associations.lib.rules_level2_base import DMSLevel2bBase

__all__ = ['Level2Association']


DEFAULT_NAME = 'CreatedByLevel2Association'


class Level2Association(dict):
    """Read in or create a Level2 association

    Parameters
    ----------
    asn: dict or Association
        An already existing association

    Notes
    -----
    This class is normally not instantiated.
    the `open` method should be used as the factory
    method to read an association or create one from
    a string or `Datamodel` object, or a list of such
    objects.
    """

    default_lvl2asn_info = {
        'program': DEFAULT_NAME,
        'target': DEFAULT_NAME,
        'asn_pool': DEFAULT_NAME
    }

    level2b_registry = AssociationRegistry(
        definition_files=[libpath('rules_level2b.py')],
        include_default=False
    )


    @classmethod
    def open(cls, obj):
        """ Open object and return a Level2 association of it

        Parameters
        ----------
        obj: Association, str, Datamodel, [str[,...]], [Datamodel[,...]]
            The obj to return as an association

        Attributes
        ----------
        Along with the attributes belonging to a Level2 association, the
        following are added:

        filename: str
            The name of the association file, if such a file
            where passed in. Otherwise a default value is given.

        Returns
        -------
        association: DMSLevel2bBase
            An association created using given obj
        """
        try:
            with open(obj) as fp:
                pure_asn = load_asn(fp, registry=cls.level2b_registry)
        except Exception:
            if not isinstance(obj, list):
                obj = [obj]
            asn = asn_from_list(
                obj,
                rule=DMSLevel2bBase,
                meta=cls.default_lvl2asn_info,
                product_name_func=cls.model_product_name
            )
            asn.filename = DEFAULT_NAME
        else:
            asn = cls(pure_asn)
            asn.filename = obj

        return asn

    @staticmethod
    def model_product_name(model):
        """Product a model product name based on the model.

        Parameters
        ----------
        model: DataModel
            The model to get the name from
        """
        product_name, ext = splitext(model.meta.filename)
        return product_name
