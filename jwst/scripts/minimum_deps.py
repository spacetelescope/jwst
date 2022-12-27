#!/usr/bin/env python

"""
Generate a requirements-min.txt file based on install_requires
"""
import warnings

import pkg_resources
import requests


def get_minimum_version(requirement):
    """Return minimum version available on PyPi for a given version specification"""
    if not requirement.specs:
        warnings.warn(
            f'No version specifier for {requirement.name} in '
            'install_requires.  Using lowest available version on PyPi.',
            stacklevel=2,
        )
    content = requests.get(
        f'https://pypi.python.org/pypi/{requirement.name}/json'
    ).json()
    versions = sorted(
        [pkg_resources.parse_version(v) for v in content['releases'].keys()]
    )
    for version in versions:
        if version in requirement:
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


def write_minimum_requirements_file(filename: str = None):
    if filename is None:
        filename = 'requirements-min.txt'

    """Write out a requirements-min.txt file for minimum package versions"""
    dist = pkg_resources.get_distribution('jwst')

    with open(filename, 'w') as fd:
        for requirement in dist.requires():
            if requirement.url is None:
                version = get_minimum_version(requirement)
                fd.write(f'{requirement.name}=={version}\n')
            else:
                fd.write(f'{requirement}\n')


if __name__ == '__main__':
    write_minimum_requirements_file()
