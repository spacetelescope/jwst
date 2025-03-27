"""Utilities for working with JWST pipeline steps."""

import importlib.util
import inspect
import logging
import os
import re
from collections.abc import Sequence
from functools import wraps
from importlib import import_module

from jwst import datamodels

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Step classes that are not user-api steps
NON_STEPS = [
    "EngDBLogStep",
    "FunctionWrapper",
    "JwstPipeline",
    "JwstStep",
    "Pipeline",
    "Step",
    "SystemCall",
]

NOT_SET = "NOT SET"
COMPLETE = "COMPLETE"
SKIPPED = "SKIPPED"


def all_steps():
    """
    List all classes subclassed from Step.

    Returns
    -------
    steps : dict
        Key is the classname, value is the class
    """
    from jwst.stpipe import Step

    jwst = import_module("jwst")
    jwst_fpath = os.path.split(jwst.__file__)[0]

    steps = {}
    for module in load_local_pkg(jwst_fpath):
        more_steps = {
            klass_name: klass
            for klass_name, klass in inspect.getmembers(
                module, lambda o: inspect.isclass(o) and issubclass(o, Step)
            )
            if klass_name not in NON_STEPS
        }
        steps.update(more_steps)

    return steps


def load_local_pkg(fpath):
    """
    Make a generator to list all modules under fpath.

    Parameters
    ----------
    fpath : str
        File path to the package to load.

    Yields
    ------
    module
        The next module found in the package.
    """
    package_fpath, package = os.path.split(fpath)
    package_fpath_len = len(package_fpath) + 1
    try:
        for module_fpath in folder_traverse(
            fpath, basename_regex=r"[^_].+\.py$", path_exclude_regex="(tests)|(regtest)"
        ):
            folder_path, fname = os.path.split(module_fpath[package_fpath_len:])
            module_path = folder_path.split("/")
            module_path.append(os.path.splitext(fname)[0])  # noqa: PTH122
            module_path = ".".join(module_path)
            try:
                spec = importlib.util.spec_from_file_location(module_path, module_fpath)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
            except Exception as err:
                logger.debug(f'Cannot load module "{module_path}": {str(err)}')
            else:
                yield module
    except Exception as err:
        logger.debug(f'Cannot complete package loading: Exception occurred: "{str(err)}"')


def folder_traverse(folder_path, basename_regex=".+", path_exclude_regex="^$"):
    """
    Traverse folder and generate full file paths for each file found.

    Parameters
    ----------
    folder_path : str
        The folder to traverse

    basename_regex : str
        Regular expression that must match
        the `basename` part of the file path.

    path_exclude_regex : str
        Regular expression to exclude a path.

    Yields
    ------
    file_path : str
        The full path to the next file
    """
    basename_regex = re.compile(basename_regex)
    path_exclude_regex = re.compile(path_exclude_regex)
    for root, _dirs, files in os.walk(folder_path):
        if path_exclude_regex.search(root):
            continue
        for file in files:
            if basename_regex.match(file):
                yield os.path.join(root, file)  # noqa: PTH118


def record_step_status(datamodel, cal_step, success=True):
    """
    Record whether or not a step completed in meta.cal_step.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel`, `~jwst.datamodels.ModelContainer`,
        `~jwst.datamodels.ModelLibrary`, str, or Path instance
        This is the datamodel or container of datamodels to modify in place

    cal_step : str
        The attribute in meta.cal_step for recording the status of the step

    success : bool
        If True, then 'COMPLETE' is recorded.  If False, then 'SKIPPED'
    """
    if success:
        status = COMPLETE
    else:
        status = SKIPPED

    if isinstance(datamodel, Sequence):
        for model in datamodel:
            model.meta.cal_step._instance[cal_step] = status  # noqa: SLF001
    elif isinstance(datamodel, datamodels.ModelLibrary):
        with datamodel:
            for model in datamodel:
                model.meta.cal_step._instance[cal_step] = status  # noqa: SLF001
                datamodel.shelve(model)
    else:
        datamodel.meta.cal_step._instance[cal_step] = status  # noqa: SLF001

    # TODO: standardize cal_step naming to point to the official step name


def query_step_status(datamodel, cal_step):
    """
    Query the status of a step in meta.cal_step.

    Parameters
    ----------
    datamodel : `~jwst.datamodels.JwstDataModel` or `~jwst.datamodels.ModelContainer` instance
        The datamodel or container of datamodels to check

    cal_step : str
        The attribute in meta.cal_step to check

    Returns
    -------
    status : str
        The status of the step in meta.cal_step, typically 'COMPLETE' or 'SKIPPED'

    Notes
    -----
    In principle, a step could set the COMPLETE status for only some subset
    of models, so checking the zeroth model instance may not always be correct.
    However, this is not currently done in the pipeline. This function should be
    updated to accommodate that use-case as needed.
    """
    if isinstance(datamodel, Sequence):
        return getattr(datamodel[0].meta.cal_step, cal_step, NOT_SET)
    else:
        return getattr(datamodel.meta.cal_step, cal_step, NOT_SET)


def invariant_filename(save_model_func):
    """
    Restore meta.filename after save_model.

    Parameters
    ----------
    save_model_func : function
        The save_model function to wrap

    Returns
    -------
    save_model : function
        The wrapped save_model function
    """

    @wraps(save_model_func)
    def save_model(model, **kwargs):
        try:
            filename = model.meta.filename
        except AttributeError:
            filename = None

        result = save_model_func(model, **kwargs)

        if filename:
            model.meta.filename = filename

        return result

    return save_model
