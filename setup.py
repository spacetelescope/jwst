import sysconfig

import numpy as np
from setuptools import Extension, setup

FREE_THREADED_PYTHON = sysconfig.get_config_var("Py_GIL_DISABLED") == 1

# Setup C module include directories
include_dirs = [np.get_include()]

# Setup C module macros
define_macros = [
    ("NUMPY", "1"),
]
if not FREE_THREADED_PYTHON:
    define_macros.append(("Py_LIMITED_API", 0x030B0000))  # PY_VERSION_HEX for 3.11

SETUPTOOLS_OPTIONS = {}
if not FREE_THREADED_PYTHON:
    SETUPTOOLS_OPTIONS["bdist_wheel"] = {"py_limited_api": "cp311"}

setup(
    # importing these extension modules is tested in `.github/workflows/build.yml`;
    # when adding new modules here, make sure to add them to the `test_command` entry there
    ext_modules=[
        Extension(
            "jwst.lib.winclip",
            ["jwst/lib/src/winclip.c"],
            include_dirs=include_dirs,
            define_macros=define_macros,
            py_limited_api=not FREE_THREADED_PYTHON,
        ),
        Extension(
            "jwst.cube_build.cube_match_internal",
            [
                "jwst/cube_build/src/cube_match_internal.c",
                "jwst/cube_build/src/cube_utils.c",
            ],
            include_dirs=include_dirs,
            define_macros=define_macros,
            py_limited_api=not FREE_THREADED_PYTHON,
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
            py_limited_api=not FREE_THREADED_PYTHON,
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
            py_limited_api=not FREE_THREADED_PYTHON,
        ),
        Extension(
            "jwst.cube_build.blot_median",
            ["jwst/cube_build/src/blot_median.c"],
            include_dirs=include_dirs,
            define_macros=define_macros,
            py_limited_api=not FREE_THREADED_PYTHON,
        ),
        Extension(
            "jwst.straylight.calc_xart",
            ["jwst/straylight/src/calc_xart.c"],
            include_dirs=include_dirs,
            define_macros=define_macros,
            py_limited_api=not FREE_THREADED_PYTHON,
        ),
    ],
    options=SETUPTOOLS_OPTIONS,
)
