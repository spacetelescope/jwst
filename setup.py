import numpy as np
from setuptools import setup, Extension

# Setup C module include directories
include_dirs = [np.get_include()]

# Setup C module macros
define_macros = [("NUMPY", "1")]

setup(
    # importing these extension modules is tested in `.github/workflows/build.yml`;
    # when adding new modules here, make sure to add them to the `test_command` entry there
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
