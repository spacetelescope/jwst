import sys
import importlib
from glob import glob
from os.path import basename
from setuptools import setup, find_packages, Extension, _install_setup_requires


try:
    from distutils.config import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.cfg'])

# Get some config values
metadata = dict(conf.items('metadata'))
PACKAGENAME = metadata['name']

options = dict(conf.items('options'))
SETUP_REQUIRES = [
    s for s in map(str.strip, options['setup_requires'].splitlines()) if s
]
INSTALL_REQUIRES = [
    s for s in map(str.strip, options['install_requires'].splitlines()) if s
]

scripts = [s for s in glob("scripts/*") if basename(s) != "__pycache__"]

# package_data values are glob patterns relative to each specific subpackage.
package_data = {
    "": [
        "*.asdf",
        "*.cfg",
        "tests/data/*.csv",
        "tests/data/*.ecsv",
        "tests/data/*.fits",
        "tests/data/**/*.fits",
        "*.json",
        "tests/data/*.json",
        "tests/data/**/*.json",
        "tests/data/*.txt",
        "*.yaml",
    ],

    "jwst.fits_generator": [
        "templates/*.inc",
        "templates/*.txt",
        "tests/okfile/*.prop",
    ],

    "jwst.lib": ["tests/data/siaf.db"],

    # Include the rules .py files in associations test data
    "jwst.associations": ["tests/data/*.py"],

    # Include C extensions
    "jwst.lib.src": ["*.c"],

    # Include the transforms schemas
    "jwst.transforms": ["schemas/stsci.edu/jwst_pipeline/*.yaml"],
    "jwst.stpipe.resources": ["schemas/*.yaml"],
}

# Install packages required for this setup to proceed:
_install_setup_requires(dict(setup_requires=SETUP_REQUIRES))
for dep_pkg in SETUP_REQUIRES:
    try:
        importlib.import_module(dep_pkg)
    except ImportError:
        print(f"{dep_pkg} is required in order to install '{PACKAGENAME}'.\n"
              f"Please install {dep_pkg} first.",
              file=sys.stderr)
        exit(1)

# Setup C module include directories
import numpy  # noqa: E402
include_dirs = [numpy.get_include()]

# Setup C module macros
define_macros = [('NUMPY', '1')]

# Handle MSVC `wcsset` redefinition
if sys.platform == 'win32':
    define_macros += [
        ('_CRT_SECURE_NO_WARNING', None),
        ('__STDC__', 1)
    ]

setup(
    use_scm_version=True,
    setup_requires=SETUP_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    scripts=scripts,
    packages=find_packages(),
    package_data=package_data,
    ext_modules=[
        Extension(
            'jwst.lib.winclip',
            ['jwst/lib/src/winclip.c'],
            include_dirs=include_dirs,
            define_macros=define_macros
        )
    ],
)
