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


class Level2Association(object):

    default_lvl2asn_info = {
        'program': 'lvl2asncreated',
        'target': 'lvl2asncreated',
        'asn_pool': 'lvl2asncreated',
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

        Returns
        -------
        association: DMSLevel2bBase
            An association created using given obj
        """
        try:
            with open(obj) as fp:
                asn = load_asn(fp, registry=cls.level2b_registry)
        except Exception:
            if not isinstance(obj, list):
                obj = [obj]
            asn = asn_from_list(
                obj,
                rule=DMSLevel2bBase,
                meta=cls.default_lvl2asn_info,
                product_name_func=cls.model_product_name
            )
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
