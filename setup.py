import numpy
from setuptools import setup, find_packages, Extension

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
        "*.cat",
        "*.hdr",
    ],
    "jwst.fits_generator": [
        "templates/*.inc",
        "templates/*.txt",
        "tests/okfile/*.prop",
    ],
    "jwst.lib": [
        "tests/data/*.asdf",
        "tests/data/*.db",
        "tests/data/*.ecsv",
        "tests/data/*.fits",
    ],
    # Include the rules .py files in associations test data
    "jwst.associations": ["tests/data/*.py"],
    # Include C extensions
    "jwst.lib.src": ["*.c"],
    "jwst.cube_build.src": ["*.c"],
    "jwst.straylight.src": ["*.c"],
    # Include the transforms schemas
    "jwst.transforms": ["resources/schemas/stsci.edu/jwst_pipeline/*.yaml"],
    "jwst.stpipe.resources": ["schemas/*.yaml"],
}

# Setup C module include directories
include_dirs = [numpy.get_include()]

# Setup C module macros
define_macros = [("NUMPY", "1")]

setup(
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    packages=find_packages(),
    package_data=package_data,
    ext_modules=[
        Extension(
            "jwst.lib.winclip",
            ["jwst/lib/src/winclip.c"],
            include_dirs=include_dirs,
            define_macros=define_macros,
        ),
        Extension(
            "jwst.cube_build.cube_match_internal",
            [
                "jwst/cube_build/src/cube_match_internal.c",
                "jwst/cube_build/src/cube_utils.c",
            ],
            include_dirs=include_dirs,
            define_macros=define_macros,
        ),
        Extension(
            "jwst.cube_build.cube_match_sky_pointcloud",
            [
                "jwst/cube_build/src/cube_match_sky_pointcloud.c",
                "jwst/cube_build/src/cube_utils.c",
                "jwst/cube_build/src/cube_dq_utils.c",
            ],
            include_dirs=include_dirs,
            define_macros=define_macros,
        ),
        Extension(
            "jwst.cube_build.cube_match_sky_driz",
            [
                "jwst/cube_build/src/cube_match_sky_driz.c",
                "jwst/cube_build/src/cube_utils.c",
                "jwst/cube_build/src/cube_dq_utils.c",
            ],
            include_dirs=include_dirs,
            define_macros=define_macros,
        ),
        Extension(
            "jwst.cube_build.blot_median",
            ["jwst/cube_build/src/blot_median.c"],
            include_dirs=include_dirs,
            define_macros=define_macros,
        ),
        Extension(
            "jwst.straylight.calc_xart",
            ["jwst/straylight/src/calc_xart.c"],
            include_dirs=include_dirs,
            define_macros=define_macros,
        ),
    ],
)
