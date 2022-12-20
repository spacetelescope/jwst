#!/usr/bin/env python

import argparse
import pkgutil
from importlib import import_module

import jwst


def main():
    """Import all packages/modules in jwst to test they are importable

    This script is used in a Travis build that doesn't install any test
    dependencies from the package. This is to verify that all the modules
    are at least importable with only the dependencies listed in
    install_requires in setup.py.

    This is to prevent adding code to the runtime codebase that have
    dependencies that are only pulled in via the test dependencies.  So for
    example, accidentally adding a runtime dependency on pytest.
    """

    parser = argparse.ArgumentParser(
        description='Check if install_requires is up-to-date')
    parser.parse_args()

    package = jwst
    prefix = package.__name__ + '.'

    no_check = ['test', 'time']

    for importer, module_name, ispkg in pkgutil.walk_packages(package.__path__, prefix=prefix):
        if not any(word in module_name for word in no_check):
            import_module(module_name)


if __name__ == '__main__':
    main()
