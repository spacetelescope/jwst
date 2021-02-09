from os.path import basename
from setuptools import setup, find_packages, Extension
from glob import glob


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
    # Include the transforms schemas
    "jwst.transforms": ["schemas/stsci.edu/jwst_pipeline/*.yaml"],
    # Include C extensions
    "jwst.grism_lib": ["*.c"],
    "jwst.stpipe.resources": ["schemas/*.yaml"],
}

setup(
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    scripts=scripts,
    packages=find_packages(),
    package_data=package_data,
    ext_modules=[Extension('jwst.grism_lib.polyclip_c', ['jwst/grism_lib/polyclip_c.c'])],
)
