"""Load an Association from a file or object."""

from inspect import isclass

from . import Association, AssociationRegistry


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
    validate : bool
        Validate against the class' defined schema, if any.
    first : bool
        A serialization potentially matches many rules.
        Only return the first successful load.
    registry : AssociationRegistry or None
        The `AssociationRegistry` to use.
        If None, no registry is used.
        Can be passed just a registry class instead of instance.
    kwargs : dict
        Other arguments to pass to the `load` methods defined
        in the `Association.IORegistry`

    Returns
    -------
    Association object

    Raises
    ------
    AssociationNotValidError
        Cannot create or validate the association.

    Notes
    -----
    The `serialized` object can be in any format
    supported by the registered I/O routines. For example, for
    `json` and `yaml` formats, the input can be either a string or
    a file object containing the string.

    If no registry is specified, the default `Association.load`
    method is used.
    """
    if registry is None:
        return Association.load(serialized, fmt=fmt, validate=validate)

    if isclass(registry):
        registry = registry()
    return registry.load(serialized, fmt=fmt, first=first, validate=validate, **kwargs)
