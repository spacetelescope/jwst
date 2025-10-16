"""JWST-specific Step and Pipeline base classes."""

import logging
from pathlib import Path

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import JwstDataModel, read_metadata
from stpipe import Pipeline, crds_client
from stpipe import Step as _Step

from jwst import __version__, __version_commit__
from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.lib.suffix import remove_suffix
from jwst.stpipe._cal_logs import _LOG_FORMATTER

log = logging.getLogger(__name__)

__all__ = ["JwstStep", "JwstPipeline"]


class JwstStep(_Step):
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

    @staticmethod
    def get_stpipe_loggers():
        """
        Get the names of loggers to configure.

        Returns
        -------
        loggers : tuple of str
            Tuple of log names to configure.
        """
        # Specify the log names for any dependencies whose
        # loggers we want to configure and for the special "py.warnings"
        # logger which is the source of warning log messages for python warnings
        return ("jwst", "stcal", "stdatamodels", "stpipe", "tweakwcs", "py.warnings")

    def load_as_level2_asn(self, obj):
        """
        Load object as an association.

        Loads the specified object into a Level2 association.
        If necessary, prepend ``Step.input_dir`` to all members.

        Parameters
        ----------
        obj : object
            Object to load as a Level2 association

        Returns
        -------
        association : object
            Association from ``jwst.associations.lib.rules_level2_base.DMSLevel2bBase``
        """
        # Prevent circular import:
        from jwst.associations.lib.update_path import update_key_value
        from jwst.associations.load_as_asn import LoadAsLevel2Asn

        asn = LoadAsLevel2Asn.load(obj, basename=self.output_file)
        update_key_value(asn, "expname", (), mod_func=self.make_input_path)
        return asn

    def load_as_level3_asn(self, obj):
        """
        Load object as an association.

        Loads the specified object into a Level3 association.
        If necessary, prepend ``Step.input_dir`` to all members.

        Parameters
        ----------
        obj : object
            Object to load as a Level3 association

        Returns
        -------
        association : object
            Association from ``jwst.associations.lib.rules_level3_base.DMS_Level3_Base``
        """
        # Prevent circular import:
        from jwst.associations.lib.update_path import update_key_value
        from jwst.associations.load_as_asn import LoadAsAssociation

        asn = LoadAsAssociation.load(obj)
        update_key_value(asn, "expname", (), mod_func=self.make_input_path)
        return asn

    def prepare_output(self, init, make_copy=None, open_models=True, open_as_type=None, **kwargs):
        """
        Open the input data as a model, making a copy if necessary.

        If the input data is a filename or path, it is opened
        and the open model is returned.

        If it is a list of models, it is opened as a ModelContainer.
        In this case, or if the input is a simple datamodel or a
        ModelContainer, a deep copy of the model/container is returned,
        in order to avoid modifying the input models.

        If the input is a ModelLibrary, it is simply returned, in order
        to avoid making unnecessary copies for performance-critical
        use cases.

        All copies are skipped if this step has a parent (i.e. it is
        called as part of a pipeline).

        Set make_copy explicitly to True or False to override the above
        behavior.

        Parameters
        ----------
        init : str, list, JwstDataModel, ModelContainer, or ModelLibrary
            Input data to open.
        make_copy : bool or None
            If True, a copy of the input will always be made.
            If False, a copy will never be made.  If None, a copy is
            conditionally made, depending on the input and whether the
            step is called in a standalone context.
        open_models : bool
            If True and the input is a filename or list of filenames,
            then datamodels.open will be called to open the input.
            If False, the input is returned as is.
        open_as_type : class or None
            If provided, the input will be opened as the specified class
            before returning. Intended for use with simple datamodel input
            only: container types and associations should be handled directly
            in the calling code.
        **kwargs
            Additional keyword arguments to pass to datamodels.open. Used
            only if the input is a str or list.

        Returns
        -------
        model : JwstDataModel or ModelContainer or ModelLibrary
            The opened datamodel(s).

        Raises
        ------
        TypeError
            If make_copy=True and the input is a type that cannot be copied.
        """
        # Check whether input contains datamodels
        copy_needed = False
        if isinstance(init, list):
            is_datamodel = [isinstance(m, datamodels.JwstDataModel) for m in init]
            if any(is_datamodel):
                # Make the list into a ModelContainer, since it contains models
                init = ModelContainer(init)
                copy_needed = True
        elif isinstance(init, (datamodels.JwstDataModel, ModelContainer)):
            copy_needed = True

        # Input might be a filename or path.
        # In that case, open it if desired.
        if not isinstance(init, (datamodels.JwstDataModel, ModelLibrary, ModelContainer)):
            if open_models:
                if open_as_type is not None:
                    # It is assumed the provided class is appropriate for the input.
                    input_models = open_as_type(init, **kwargs)
                else:
                    input_models = datamodels.open(init)
            else:
                # Return the filename or path -
                # the calling code will handle opening it as needed.
                input_models = init
        elif isinstance(init, datamodels.JwstDataModel):
            # Simple data model: update the datamodel type if needed
            if open_as_type is not None and type(init) is not open_as_type:
                # This will make a shallow copy.
                input_models = open_as_type(init, **kwargs)
            else:
                # Otherwise use the init model directly
                input_models = init
        else:
            # ModelContainer or ModelLibrary: use the init model directly.
            input_models = init

        # Make a deep copy if needed
        if make_copy is None:
            make_copy = copy_needed and self.parent is None
        if make_copy:
            try:
                input_models = input_models.copy()
            except AttributeError:
                # This should only happen if make_copy is explicitly set to
                # True and the input is a string or a ModelLibrary.
                raise TypeError(
                    f"Copy is not possible for input type {type(input_models)}"
                ) from None

        return input_models

    def finalize_result(self, result, reference_files_used):
        """
        Update the result with the software version and reference files used.

        Parameters
        ----------
        result : `~stdatamodels.DataModel`
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

    def run(self, *args):
        """
        Run the step.

        Parameters
        ----------
        *args
            Arguments passed to `stpipe.Step.run`.

        Returns
        -------
        result : Any
            The step output
        """
        result = super().run(*args)
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
        result : `~stdatamodels.DataModel`
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
