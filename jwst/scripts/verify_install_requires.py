#!/usr/bin/env python

import argparse

import pytest

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
        description="Check if install_requires is up-to-date"
    )
    parser.parse_args()

    pytest.main(
        [
            f"{jwst.__path__}/tests/test_import.py",
            f"-c {jwst.__path__}/tests/empty_config/pytest.ini",
        ]
    )


if __name__ == "__main__":
    main()
