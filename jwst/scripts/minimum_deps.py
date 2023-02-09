#!/usr/bin/env python

"""
Generate a requirements-min.txt file based on install_requires
"""
import warnings

import requests
from contextlib import suppress
from importlib_metadata import requires
from packaging.version import parse, InvalidVersion
from packaging.requirements import Requirement


def get_minimum_version(requirement):
    """Return minimum version available on PyPi for a given version specification"""
    if not requirement.specifier:
        warnings.warn(
            f'No version specifier for {requirement.name} in '
            'install_requires.  Using lowest available version on PyPi.',
            stacklevel=2,
        )

    content = requests.get(
        f'https://pypi.python.org/pypi/{requirement.name}/json'
    ).json()
    versions = []
    for v in content['releases'].keys():
        with suppress(InvalidVersion):
            versions.append(parse(v))

    versions = sorted(versions)

    for version in versions:
        if version in requirement.specifier:
            # If the requirement does not list any version, the lowest will be
            # returned
            return version

    # If the specified version does not exist on PyPi, issue a warning
    # and return the lowest available version
    warnings.warn(
        f'Exact version specified in {requirement} not found '
        'on PyPi.  Using lowest available version.',
        stacklevel=2,
    )
    return versions[0]


def write_minimum_requirements_file(filename: str = None, extras: list = None):
    """Write out a requirements-min.txt file for minimum package versions"""
    if filename is None:
        filename = 'requirements-min.txt'

    if extras is None:
        extras = []

    with open(filename, 'w') as fd:
        for r in requires('jwst'):
            requirement = Requirement(r)

            if requirement.marker is None or any(requirement.marker.evaluate({'extra': e}) for e in extras):
                if requirement.url is None:
                    version = get_minimum_version(requirement)
                    fd.write(f'{requirement.name}=={version}\n')
                else:
                    fd.write(f'{requirement}\n')


if __name__ == '__main__':
    write_minimum_requirements_file()
