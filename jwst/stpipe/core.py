"""JWST-specific Step and Pipeline base classes."""

from functools import wraps
import logging
from pathlib import Path

from stdatamodels.jwst.datamodels import JwstDataModel, read_metadata
from stdatamodels.jwst import datamodels
from stpipe import crds_client, Step, Pipeline

from jwst import __version_commit__, __version__
from jwst.datamodels import ModelLibrary, ModelContainer
from ._cal_logs import _LOG_FORMATTER
from jwst.lib.suffix import remove_suffix


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

    @classmethod
    def _get_crds_parameters(cls, dataset):
        """
        Get CRDS parameters for the given dataset.

        If the input dataset is a filename, achieve this by lazy-loading its metadata.

        Parameters
        ----------
        dataset : str
            The name of the dataset.

        Returns
        -------
        dict
            A dictionary of CRDS parameters.
        str
            The name of the observatory.
        """
        crds_observatory = "jwst"

        # list or container: just set to zeroth model
        # this is what stpipe does internally for ModelContainer already
        if isinstance(dataset, (list, tuple, ModelContainer)):
            if len(dataset) == 0:
                raise ValueError(f"Input dataset {dataset} is empty")
            dataset = dataset[0]

        # Already open: use the model's method to get CRDS parameters
        if isinstance(dataset, (ModelLibrary, JwstDataModel)):
            return (
                dataset.get_crds_parameters(),
                crds_observatory,
            )

        # If we get here, we had better have a filename
        if isinstance(dataset, str):
            dataset = Path(dataset)
        if not isinstance(dataset, Path):
            raise TypeError(f"Cannot get CRDS parameters for {dataset} of type {type(dataset)}")

        # for associations, open as ModelLibrary, which supports lazy-loading
        if dataset.suffix.lower() == ".json":
            model = ModelLibrary(dataset, asn_n_members=1, asn_exptypes=["science"])
            return (model.get_crds_parameters(), crds_observatory)

        # for all other cases, use read_metadata directly to lazy-load
        return (read_metadata(dataset, flatten=True), crds_observatory)

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
        from jwst.associations.load_as_asn import LoadAsLevel2Asn
        from jwst.associations.lib.update_path import update_key_value

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
        from jwst.associations.load_as_asn import LoadAsAssociation
        from jwst.associations.lib.update_path import update_key_value

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
