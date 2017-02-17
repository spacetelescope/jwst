"""Load an Association from a file or object"""
from . import AssociationRegistry


def load_asn(
        serialized,
        format=None,
        first=True,
        validate=True,
        registry=None,
        **kwargs
):
    """Load an Association from a file or object

    Parameters
    ----------
    serialized: object
        The serialized form of the association.

    format: str or None
        The format to force. If None, try all available.

    validate: bool
        Validate against the class' defined schema, if any.

    first: bool
        A serialization potentially matches many rules.
        Only return the first succesful load.

    registry: AssociationRegistry or None
        The `AssociationRegistry` to use.
        If None, the default registry is used.

    kwargs: dict
        Other arguments to pass to the `load` methods defined
        in the `Association.IORegistry`

    Returns
    -------
    The Association object

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
    """
    if registry is None:
        registry = AssociationRegistry()
    return registry.load(
        serialized,
        format=format,
        first=first,
        validate=validate,
        **kwargs
    )
