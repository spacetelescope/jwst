"""Load an Association from a file or object."""

import logging
import warnings
from inspect import isclass
from pathlib import Path

from jwst.associations import Association, AssociationRegistry
from jwst.associations.exceptions import AssociationNotValidError

log = logging.getLogger(__name__)

__all__ = ["load_asn"]


def load_asn(
    serialized, fmt=None, first=True, validate=True, registry=AssociationRegistry, **kwargs
):
    """
    Load an Association from a file or object.

    Parameters
    ----------
    serialized : object
        The serialized form of the association.
    fmt : str or None
        The format to force. If None, try all available.
    first : bool
        A serialization potentially matches many rules.
        Only return the first successful load.
    validate : bool
        Validate against the class's defined schema, if any.
    registry : AssociationRegistry or None
        The `~jwst.associations.AssociationRegistry` to use.
        If None, no registry is used.
        Can be passed just a registry class instead of instance.
    **kwargs : dict
        Other arguments to pass to the ``load`` methods defined
        in the `~jwst.associations.lib.keyvalue_registry.KeyValueRegistry`

    Returns
    -------
    Association object
        The loaded association.

    Raises
    ------
    AssociationNotValidError
        Cannot create or validate the association.

    Notes
    -----
    The serialized object can be in any format
    supported by the registered I/O routines. For example, for
    JSON format, the input can be either a string or
    a file object containing the string.

    If no registry is specified, the default
    :meth:`~jwst.associations.Association.load` method is used.
    """
    fname = getattr(serialized, "name", None)
    if fname is not None:
        suffix = Path(fname).suffix.replace(".", "")
        if suffix in ("yaml", "yml"):
            msg = (
                "Support for associations as YAML files is deprecated. "
                "Please use JSON format with extension .json instead."
            )
            warnings.warn(
                msg,
                DeprecationWarning,
                stacklevel=2,
            )
            log.warning(msg)
        elif suffix != "json":
            msg = (
                f"File extension '{suffix}' is not recognized as JSON. "
                "Attempting to load anyway, but this behavior is deprecated and will be removed "
                "in a future release. Please ensure association files have a .json extension."
            )
            warnings.warn(
                msg,
                DeprecationWarning,
                stacklevel=2,
            )
            log.warning(msg)

    try:
        asn = _do_load(
            serialized, fmt="json", first=first, validate=validate, registry=registry, **kwargs
        )
        if fmt == "yaml":
            raise AssociationNotValidError(  # noqa: TRY301
                "Association file is valid JSON, but YAML format was forced."
            )
    except AssociationNotValidError:
        # if JSON load fails, try YAML load to preserve deprecated behavior
        # ignore deprecation warning coming from yaml.load since a more specific deprecationwarning
        # will already be emitted somewhere in this function
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message="Support for associations as YAML files is deprecated",
            )
            asn = _do_load(
                serialized, fmt="yaml", first=first, validate=validate, registry=registry, **kwargs
            )
        # if we successfully loaded as YAML, emit a warning about deprecation of YAML support
        if suffix == "json":
            msg = (
                "Association file has json suffix but is invalid JSON or is force-loaded as YAML. "
                "Attempting to load as YAML, but YAML support is deprecated and will be removed "
                "in a future release. In the past, invalid JSON "
                "(e.g. due to extra trailing commas) was sometimes quietly loaded as YAML, "
                "and this behavior will no longer be supported. "
                "Please fix any JSON formatting issues in the association file."
            )
            warnings.warn(
                msg,
                DeprecationWarning,
                stacklevel=2,
            )
            log.warning(msg)
    return asn


def _do_load(
    serialized, fmt=None, first=True, validate=True, registry=AssociationRegistry, **kwargs
):
    if registry is None:
        asn = Association.load(serialized, fmt=fmt, validate=validate)
    else:
        if isclass(registry):
            registry = registry()
        asn = registry.load(serialized, fmt=fmt, first=first, validate=validate, **kwargs)
    return asn
