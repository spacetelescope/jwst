from pathlib import Path
import sys

if sys.version_info < (3, 9):
    import importlib_resources
else:
    import importlib.resources as importlib_resources


from asdf.resource import DirectoryResourceMapping
from jwst import transforms


def get_resource_mappings():
    resources_root = importlib_resources.files(transforms) / "resources"
    if not resources_root.is_dir():
        # In an editable install, the resources directory will exist off the
        # repository root:
        resources_root = Path(__file__).parent.parent.parent / "resources"
        if not resources_root.is_dir():
            raise RuntimeError("Missing resources directory")

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
