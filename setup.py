import sys
import sysconfig
from typing import TypedDict

import numpy as np
from setuptools import Extension, setup

GIL_ENABLED_PYTHON = sysconfig.get_config_var("Py_GIL_DISABLED") != 1


class CompilerOptions(TypedDict):
    """Compiler options for the C extension modules."""

    include_dirs: list[str]
    define_macros: list[tuple[str, str | None]]
    libraries: list[str]
    extra_compile_args: list[str]
    py_limited_api: bool


cfg: CompilerOptions = {
    "include_dirs": [np.get_include()],
    "define_macros": [
        ("NUMPY", "1"),
    ],
    "libraries": [],
    "extra_compile_args": [],
    "py_limited_api": GIL_ENABLED_PYTHON,
}

if GIL_ENABLED_PYTHON:
    cfg["define_macros"].append(("Py_LIMITED_API", "0x030C0000"))  # PY_VERSION_HEX for 3.12

SETUPTOOLS_OPTIONS = {}
if GIL_ENABLED_PYTHON:
    SETUPTOOLS_OPTIONS["bdist_wheel"] = {"py_limited_api": "cp312"}

if sys.platform == "win32":
    cfg["define_macros"].extend(
        [
            ("WIN32", None),
            ("__STDC__", "1"),
            ("_CRT_SECURE_NO_WARNINGS", None),
        ]
    )
else:
    cfg["libraries"].append("m")
    cfg["extra_compile_args"] = [
        "-std=c11",
        "-O3",
        "-Wall",
        "-Wextra",
        "-Werror",
        "-Wpedantic",
        "-Wno-unused-parameter",
        "-Wundef",
        # Extra flags for stricter warnings can be enabled if needed
        # "-Wconversion",
        # "-Wsign-conversion",
        # "-Wcast-qual",
        # "-Wcast-align",
        # "-Wshadow",
        # "-Wstrict-prototypes",
        # "-Wmissing-prototypes",
    ]

setup(
    # importing these extension modules is tested in `.github/workflows/build.yml`;
    # when adding new modules here, make sure to add them to the `test_command` entry there
    ext_modules=[
        Extension(
            "jwst.lib.winclip",
            ["jwst/lib/src/winclip.c"],
            **cfg,
        ),
        Extension(
            "jwst.cube_build.cube_match_internal",
            [
                "jwst/cube_build/src/cube_match_internal.c",
                "jwst/cube_build/src/cube_utils.c",
            ],
            **cfg,
        ),
        Extension(
            "jwst.cube_build.cube_match_sky_pointcloud",
            [
                "jwst/cube_build/src/cube_match_sky_pointcloud.c",
                "jwst/cube_build/src/cube_utils.c",
            ],
            **cfg,
        ),
        Extension(
            "jwst.cube_build.cube_match_sky_driz",
            [
                "jwst/cube_build/src/cube_match_sky_driz.c",
                "jwst/cube_build/src/cube_utils.c",
            ],
            **cfg,
        ),
        Extension(
            "jwst.cube_build.blot_median",
            ["jwst/cube_build/src/blot_median.c"],
            **cfg,
        ),
        Extension(
            "jwst.straylight.calc_xart",
            ["jwst/straylight/src/calc_xart.c"],
            **cfg,
        ),
    ],
    options=SETUPTOOLS_OPTIONS,
)
