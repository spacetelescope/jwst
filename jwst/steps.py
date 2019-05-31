import importlib
from pkgutil import iter_modules

import jwst

submodules = [mod for _, mod, ispkg in iter_modules(jwst.__path__) if ispkg]
# timeconversion requires non-standard dependencies and enviroment settings
submodules.remove('timeconversion')
for module in submodules:
    importlib.import_module("jwst." + module)
