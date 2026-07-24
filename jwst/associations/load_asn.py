"""Load an Association from a file or object."""

import logging
from inspect import isclass
from pathlib import Path

from astropy.utils.decorators import deprecated_renamed_argument

from jwst.associations import Association, AssociationRegistry
from jwst.associations.exceptions import AssociationNotValidError

log = logging.getLogger(__name__)

__all__ = ["load_asn"]


@deprecated_renamed_argument("fmt", None, since="2.1")
def load_asn(
    serialized,
    fmt=None,  # noqa: ARG001
    first=True,
    validate=True,
    registry=AssociationRegistry,
    **kwargs,
):
    """
    Load an Association from a file or object.

    Parameters
    ----------
    serialized : object
        The serialized form of the association.

    fmt : str or None
        The format to force. If None, try all available.

        .. version-deprecated:: 2.1
            Only JSON format is supported now.

    first : bool
        A serialization potentially matches many rules.
        Only return the first successful load.

    validate : bool
        Validate against the class's defined schema, if any.

    registry : `~jwst.associations.AssociationRegistry` or None
        The `~jwst.associations.AssociationRegistry` to use.
        If None, no registry is used.
        Can be passed just a registry class instead of instance.

    **kwargs : dict
        Other arguments to pass to the ``load`` methods defined
        in the `~jwst.associations.lib.keyvalue_registry.KeyValueRegistry`

    Returns
    -------
    `~jwst.associations.Association`
        The loaded association.

    Raises
    ------
    jwst.associations.exceptions.AssociationNotValidError
        Cannot create or validate the association.

    Notes
    -----
    While the serialized object must be in
    JSON format, the input can be either a string or
    a file object containing the string.

    If no registry is specified, the default
    :meth:`~jwst.associations.Association.load` method is used.
    """
    fname = getattr(serialized, "name", None)
    if fname is not None:
        suffix = Path(fname).suffix.replace(".", "")
        if suffix != "json":
            msg = (
                f"File extension '{suffix}' is not recognized as JSON. "
                "Please ensure association files have a .json extension."
            )
            raise AssociationNotValidError(msg)

    return _do_load(serialized, first=first, validate=validate, registry=registry, **kwargs)


def _do_load(serialized, first=True, validate=True, registry=AssociationRegistry, **kwargs):
    if registry is None:
        asn = Association.load(serialized, validate=validate)
    else:
        if isclass(registry):
            registry = registry()
        asn = registry.load(serialized, first=first, validate=validate, **kwargs)
    return asn
