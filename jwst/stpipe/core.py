"""
JWST-specific Step and Pipeline base classes.
"""
from collections.abc import Sequence
from stdatamodels import DataModel

from .. import __version_commit__, __version__
from .. import datamodels

from stpipe import crds_client
from stpipe import Step
from stpipe import Pipeline
from ..lib.suffix import remove_suffix
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class JwstStep(Step):

    spec = """
    output_ext = string(default='.fits')  # Output file type
    """

    @classmethod
    def _datamodels_open(cls, init, **kwargs):
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

    def finalize_result(self, result, reference_files_used):
        if isinstance(result, DataModel):
            result.meta.calibration_software_revision = __version_commit__ or 'RELEASE'
            result.meta.calibration_software_version = __version__

            if len(reference_files_used) > 0:
                for ref_name, filename in reference_files_used:
                    if hasattr(result.meta.ref_file, ref_name):
                        getattr(result.meta.ref_file, ref_name).name = filename
                result.meta.ref_file.crds.sw_version = crds_client.get_svn_version()
                result.meta.ref_file.crds.context_used = crds_client.get_context_used(result.crds_observatory)
                if self.parent is None:
                    log.info(f"Results used CRDS context: {result.meta.ref_file.crds.context_used}")

    def record_step_status(self, datamodel, cal_step, success=True):
        """Record whether or not a step completed in meta.cal_step

        Parameters
        ----------
        datamodel : `~jwst.datamodels.JwstDataModel` instance
            This is the datamodel or container of datamodels to modify in place

        cal_step : str
            The attribute in meta.cal_step for recording the status of the step

        success : bool
            If True, then 'COMPLETE' is recorded.  If False, then 'SKIPPED'
        """
        if success:
            status = 'COMPLETE'
        else:
            status = 'SKIPPED'
            self.skip = True

        if isinstance(datamodel, Sequence):
            for model in datamodel:
                model.meta.cal_step._instance[cal_step] = status
        else:
            datamodel.meta.cal_step._instance[cal_step] = status

        # TODO: standardize cal_step naming to point to the official step name

    def remove_suffix(self, name):
        return remove_suffix(name)


# JwstPipeline needs to inherit from Pipeline, but also
# be a subclass of JwstStep so that it will pass checks
# when constructing a pipeline using JwstStep class methods.
class JwstPipeline(Pipeline, JwstStep):
    def finalize_result(self, result, reference_files_used):
        if isinstance(result, DataModel):
            log.info(f"Results used CRDS context: {crds_client.get_context_used(result.crds_observatory)}")
