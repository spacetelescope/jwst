"""
Implementation of Step methods related to ASDF config files.  This module
will eventually fully replace config_parser.py, but we'll need to maintain
both until we replace configobj with traitlets.
"""
from copy import deepcopy
from datetime import datetime

from .utilities import get_fully_qualified_class_name

import asdf


_CONFIG_SCHEMA_URI = "http://stsci.edu/schemas/stpipe/step_config-1.0.0"
_LEGACY_CONFIG_SCHEMA_URI = "http://stsci.edu/schemas/stpipe/step_config-0.1.0"

# The config schema validates that this substring is not present:
_TEMPLATE_PLACEHOLDER = "<SPECIFY>"

_META_TEMPLATE = {
    "author": _TEMPLATE_PLACEHOLDER,
    "date": _TEMPLATE_PLACEHOLDER,
    "description": f"Parameters for calibration step {_TEMPLATE_PLACEHOLDER}",
    "instrument": {
        "name": _TEMPLATE_PLACEHOLDER,
    },
    "origin": _TEMPLATE_PLACEHOLDER,
    "pedigree": _TEMPLATE_PLACEHOLDER,
    "reftype": _TEMPLATE_PLACEHOLDER,
    "telescope": _TEMPLATE_PLACEHOLDER,
    "useafter": _TEMPLATE_PLACEHOLDER,
}


class StepConfig:
    """
    Step configuration container.

    Parameters
    ----------
    class_name : str
        Fully-qualified Step subclass name.
    name : str
        Nickname
    parameters : dict
        Parameters indexed by parameter name.
    steps : list of StepConfig
        List of sub-step configs.
    """
    def __init__(self, class_name, name, parameters, steps):
        self._class_name = class_name
        self._name = name
        self._parameters = parameters
        self._steps = steps

    @property
    def class_name(self):
        return self._class_name

    @property
    def name(self):
        return self._name

    @property
    def parameters(self):
        return self._parameters

    @property
    def steps(self):
        return self._steps

    def __eq__(self, other):
        if not isinstance(other, StepConfig):
            return False

        return (
            self.class_name == other.class_name and
            self.name == other.name and
            self.parameters == other.parameters and
            self.steps == other.steps
        )

    def _to_tree(self):
        tree = {
            "class": self.class_name,
            "name": self.name,
            "parameters": self.parameters,
        }
        if len(self.steps) > 0:
            tree["steps"] = [s._to_tree() for s in self.steps]
        return tree

    def to_asdf(self, include_metadata=False):
        """
        Convert this config to an AsdfFile, which may be
        used to write the config to a file.

        Parameters
        ----------
        include_metadata : bool, optional
            Set to True to include metadata that is required
            for submission to CRDS.

        Returns
        -------
        asdf.AsdfFile
        """
        result = asdf.AsdfFile(self._to_tree())

        if include_metadata:
            meta = deepcopy(_META_TEMPLATE)
            meta["date"] = meta["date"].replace(_TEMPLATE_PLACEHOLDER, datetime.utcnow().replace(microsecond=0).isoformat())
            meta["description"] = meta["description"].replace(_TEMPLATE_PLACEHOLDER, self.class_name)
            result["meta"] = meta

        _validate_asdf(result, _CONFIG_SCHEMA_URI)

        return result

    @classmethod
    def from_asdf(cls, asdf_file):
        """
        Create a StepConfig from an open AsdfFile.

        Parameters
        ----------
        asdf_file : asdf.AsdfFile

        Returns
        -------
        stpipe.config.StepConfig

        Raises
        ------
        asdf.ValidationError
            If the file does not validate against any of the
            recognized schemas.
        """
        try:
            _validate_asdf(asdf_file, _CONFIG_SCHEMA_URI)
            return cls._from_tree(asdf_file.tree)
        except asdf.ValidationError as e:
            # Maybe it's an old-style ASDF file:
            try:
                _validate_asdf(asdf_file, _LEGACY_CONFIG_SCHEMA_URI)
                return cls._from_legacy_tree(asdf_file.tree)
            except asdf.ValidationError:
                # Raise the original error so that we encourage
                # use of the new config file format.
                raise e

    @classmethod
    def _from_tree(cls, tree):
        class_name = tree["class"]
        name = tree["name"]
        parameters = deepcopy(tree["parameters"])
        steps = [cls._from_tree(s) for s in tree.get("steps", [])]
        return StepConfig(class_name, name, parameters, steps)

    @classmethod
    def _from_legacy_tree(cls, tree):
        def _inner(parameters, name):
            class_name = parameters.pop("class", None)
            steps = [_inner(s, n) for n, s in parameters.pop("steps", {}).items()]
            return StepConfig(class_name, name, parameters, steps)

        parameters = deepcopy(tree["parameters"])
        name = parameters.pop("name")
        return _inner(parameters, name)


def export_config(step):
    """
    Export a step's current parameters to a StepConfig object.

    Parameters
    ----------
    step : stpipe.Step

    Returns
    -------
    stpipe.config.StepConfig
    """
    class_name = get_fully_qualified_class_name(step)
    name = step.name

    parameters = step.get_pars()
    # The Pipeline class includes step parameters, but we're
    # going to collect those ourselves so we can also get
    # ahold of the name and class name.
    parameters.pop("steps", None)

    steps = [
        export_config(getattr(step, step_name))
        for step_name, _ in getattr(step, "step_defs", {}).items()
    ]

    return StepConfig(class_name, name, parameters, steps)


def _validate_asdf(asdf_file, schema_uri):
    """
    Validate an ASDF file against the config schema.

    Parameters
    ----------
    asdf_file : asdf.AsdfFile
        Open config file.
    schema_uri : str
        StepConfig schema URI (must be registered with the ASDF entry point).

    Raises
    ------
    asdf.ValidationError
    """
    # TODO: We should add a method on AsdfFile to facilitate this:
    schema = asdf.schema.load_schema(schema_uri, asdf_file.resolver)
    tagged_tree = asdf.yamlutil.custom_tree_to_tagged_tree(asdf_file.tree, asdf_file)
    asdf.schema.validate(tagged_tree, asdf_file, schema)
