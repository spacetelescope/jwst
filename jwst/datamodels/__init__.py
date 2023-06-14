# The jwst.datamodels submodule was moved to stdatamodels.jwst.datamodels
# https://github.com/spacetelescope/jwst/pull/7439

import importlib
from inspect import ismodule
import sys

from stdatamodels.jwst.datamodels.util import open

from .container import ModelContainer
from .source_container import SourceModelContainer

import stdatamodels.jwst.datamodels

# Import everything defined in stdatamodels.jwst.datamodels.__all__
from stdatamodels.jwst.datamodels import * # noqa: F403

# Define __all__ to include stdatamodels.jwst.datamodels.__all__
__all__ = [
    'open',
    'ModelContainer', 'SourceModelContainer',
] + stdatamodels.jwst.datamodels.__all__


# Modules that are not part of stdatamodels
_jwst_modules = ["container", "source_container"]

# Models that are not part of stdatamodels
_jwst_models = ["ModelContainer", "SourceModelContainer"]

# Deprecated modules in stdatamodels
_deprecated_modules = ['drizproduct', 'multiprod']

# Deprecated models in stdatamodels
_deprecated_models = ['DrizProductModel', 'MultiProductModel', 'MIRIRampModel']

# Import all submodules from stdatamodels.jwst.datamodels
for attr in dir(stdatamodels.jwst.datamodels):
    if attr[0] == '_':
        continue
    if attr in _jwst_models or attr in _deprecated_modules or attr in _deprecated_models:
        continue
    obj = getattr(stdatamodels.jwst.datamodels, attr)
    if ismodule(obj):
        # Make the submodule available locally
        locals()[attr] = obj
        # Add the submodule to sys.modules so that a call
        # to jwst.datamodels.dqflags will return the submodule
        # stdatamodels.jwst.datamodels.dqflags
        sys.modules[f"jwst.datamodels.{attr}"] = obj

# Add a few submodules to sys.modules without exposing them locally
for _submodule_name in ['schema', 'schema_editor', 'validate']:
    _submodule = importlib.import_module(f"stdatamodels.jwst.datamodels.{_submodule_name}")
    sys.modules[f"jwst.datamodels.{_submodule_name}"] = _submodule
