"""
This module supports the entry points for ASDF support for the `jwst.transforms`.
"""

import sys

if sys.version_info < (3, 9):
    import importlib_resources
else:
    import importlib.resources as importlib_resources


from asdf.resource import DirectoryResourceMapping
from jwst import transforms


def get_resource_mappings():
    """
    Get the `jwst.transforms` resource mappings, that is the schemas for the transforms.

    This method is registered with the `asdf.resource_mappings` entry point for
    the `jwst_pipeline`.

    Returns
    -------
    list of the `asdf.resource.ResourceMapping` instances containing the `jwst.transforms`
    schemas.
    """
    resources_root = importlib_resources.files(transforms) / "resources"
    if not resources_root.is_dir():
        raise RuntimeError(f"Missing resources directory: {resources_root=}")

    return [
        DirectoryResourceMapping(
            resources_root / "schemas" / "stsci.edu" / "jwst_pipeline",
            "http://stsci.edu/schemas/jwst_pipeline/",
        ),
        DirectoryResourceMapping(
            resources_root / "manifests",
            "asdf://stsci.edu/jwst_pipeline/manifests/",
        ),
    ]


def get_extensions():
    """
    Get the jwst.transforms extension.
    This method is registered with the asdf.extensions entry point.

    Returns
    -------
    list of asdf.extension.Extension
    """
    from . import extensions
    return extensions.TRANSFORM_EXTENSIONS
