"""JWST-specific Step and Pipeline base classes."""

from functools import wraps
import logging
import warnings

from stdatamodels.jwst.datamodels import JwstDataModel
from stdatamodels.jwst import datamodels
from stpipe import crds_client, Step, Pipeline

from jwst import __version_commit__, __version__
from ._cal_logs import _LOG_FORMATTER
from ..lib.suffix import remove_suffix


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class JwstStep(Step):
    """A JWST pipeline step."""

    spec = """
    output_ext = string(default='.fits')  # Output file type
    """  # noqa: E501

    _log_records_formatter = _LOG_FORMATTER

    @classmethod
    def _datamodels_open(cls, init, **kwargs):
        return datamodels.open(init, **kwargs)

    def load_as_level2_asn(self, obj):
        """
        Load object as an association.

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
        update_key_value(asn, "expname", (), mod_func=self.make_input_path)
        return asn

    def load_as_level3_asn(self, obj):
        """
        Load object as an association.

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
        update_key_value(asn, "expname", (), mod_func=self.make_input_path)
        return asn

    def finalize_result(self, result, reference_files_used):
        """
        Update the result with the software version and reference files used.

        Parameters
        ----------
        result : `~jwst.datamodels.DataModel`
            The output data model to be updated.
        reference_files_used : list of tuple
            The names and file paths of reference files used.
        """
        if isinstance(result, JwstDataModel):
            result.meta.calibration_software_revision = __version_commit__ or "RELEASE"
            result.meta.calibration_software_version = __version__

            if len(reference_files_used) > 0:
                for ref_name, filename in reference_files_used:
                    if hasattr(result.meta.ref_file, ref_name):
                        getattr(result.meta.ref_file, ref_name).name = filename
                result.meta.ref_file.crds.sw_version = crds_client.get_svn_version()
                result.meta.ref_file.crds.context_used = crds_client.get_context_used(
                    result.crds_observatory
                )
                if self.parent is None:
                    log.info(f"Results used CRDS context: {result.meta.ref_file.crds.context_used}")

            if self.class_alias:
                if not hasattr(result, "cal_logs"):
                    result.cal_logs = {}
                setattr(result.cal_logs, self.class_alias, self._log_records)

    def remove_suffix(self, name):
        """
        Remove the suffix if a known suffix is already in name.

        Parameters
        ----------
        name : str
            The name to remove the suffix from.

        Returns
        -------
        name : str
            The name with the suffix removed.
        """
        return remove_suffix(name)

    @wraps(Step.run)
    def run(self, *args, **kwargs):
        """
        Run the step.

        Parameters
        ----------
        *args
            Arguments passed to `stpipe.Step.run`.
        **kwargs
            Keyword arguments passed to `stpipe.Step.run`.

        Returns
        -------
        result : Any
            The step output
        """
        result = super().run(*args, **kwargs)
        if not self.parent:
            log.info(f"Results used jwst version: {__version__}")
        return result

    @wraps(Step.__call__)
    def __call__(self, *args, **kwargs):  # numpydoc ignore=RT01
        """Deprecated method. Use `run` instead."""  # noqa: D401
        if not self.parent:
            warnings.warn(
                "Step.__call__ is deprecated. It is equivalent to Step.run "
                "and is not recommended. See "
                "https://jwst-pipeline.readthedocs.io/en/latest/jwst/"
                "user_documentation/running_pipeline_python.html"
                "#advanced-use-pipeline-run-vs-pipeline-call for more details.",
                UserWarning,
                stacklevel=2,
            )
        return super().__call__(*args, **kwargs)


class JwstPipeline(Pipeline, JwstStep):
    """
    A JWST pipeline.

    JwstPipeline needs to inherit from Pipeline, but also
    be a subclass of JwstStep so that it will pass checks
    when constructing a pipeline using JwstStep class methods.
    """

    def finalize_result(self, result, _reference_files_used):
        """
        Update the result with the software version and reference files used.

        Parameters
        ----------
        result : `~jwst.datamodels.DataModel`
            The output data model to be updated.
        _reference_files_used : list of tuple
            The names and file paths of reference files used.
        """
        if isinstance(result, JwstDataModel):
            log.info(
                "Results used CRDS context: "
                f"{crds_client.get_context_used(result.crds_observatory)}"
            )

            if self.class_alias:
                if not hasattr(result, "cal_logs"):
                    result.cal_logs = {}

                # remove the step logs as they're captured by the pipeline log
                for _, step in self.step_defs.items():
                    if hasattr(result.cal_logs, step.class_alias):
                        delattr(result.cal_logs, step.class_alias)

                setattr(result.cal_logs, self.class_alias, self._log_records)
