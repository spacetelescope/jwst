"""
This module supports the entry points for ASDF support for the `jwst.datamodels`.
"""

import sys

if sys.version_info < (3, 9):
    import importlib_resources
else:
    import importlib.resources as importlib_resources


from asdf.resource import DirectoryResourceMapping
from jwst import datamodels


def get_resource_mappings():
    """
    Get the `jwst.datamodels` resource mappings, that is the schemas for the datamodels.

    This method is registered with the `asdf.resource_mappings` entry point for
    the `jwst_datamodel`.

    Returns
    -------
    list of the `asdf.resource.ResourceMapping` instances containing the `jwst.datamodels`
    schemas.
    """
    resources_root = importlib_resources.files(datamodels)
    if not resources_root.is_dir():
        raise RuntimeError(f"Missing resources directory: {resources_root=}")

    return [
        DirectoryResourceMapping(
            resources_root / "schemas",
            "http://stsci.edu/schemas/jwst_datamodel/",
        ),
        DirectoryResourceMapping(
            resources_root / "metaschema",
            "http://stsci.edu/schemas/fits-schema/",
        )
    ]
