"""
JWST-specific Step and Pipeline base classes.
"""
from .step import Step
from .pipeline import Pipeline
from .. import datamodels


class JwstStep(Step):
    @classmethod
    def datamodels_open(cls, init, **kwargs):
        return datamodels.open(init, **kwargs)

    def load_as_level2_asn(self, obj):
        """Load object as an association

        Loads the specified object into a Level2 association.
        If necessary, prepend `Step.input_dir` to all members.

        Parameters
        ----------
        obj : object
            Object to load as a Level2 association

        Returns
        -------
        association : jwst.associations.lib.rules_level2_base.DMSLevel2bBase
            Association
        """
        # Prevent circular import:
        from ..associations.load_as_asn import LoadAsLevel2Asn
        from ..associations.lib.update_path import update_key_value

        asn = LoadAsLevel2Asn.load(obj, basename=self.output_file)
        update_key_value(asn, 'expname', (), mod_func=self.make_input_path)
        return asn

    def load_as_level3_asn(self, obj):
        """Load object as an association

        Loads the specified object into a Level3 association.
        If necessary, prepend `Step.input_dir` to all members.

        Parameters
        ----------
        obj : object
            Object to load as a Level3 association

        Returns
        -------
        association : jwst.associations.lib.rules_level3_base.DMS_Level3_Base
            Association
        """
        # Prevent circular import:
        from ..associations.load_as_asn import LoadAsAssociation
        from ..associations.lib.update_path import update_key_value

        asn = LoadAsAssociation.load(obj)
        update_key_value(asn, 'expname', (), mod_func=self.make_input_path)
        return asn


# JwstPipeline needs to inherit from Pipeline, but also
# be a subclass of JwstStep so that it will pass checks
# when constructing a pipeline using JwstStep class methods.
# It's important that Pipeline occur first so that it
# takes precedence in the resolution order.
class JwstPipeline(Pipeline, JwstStep):
    pass
