"""Create an association from a list"""

from .lib.rules_level3_base import DMS_Level3_Base


def asn_from_list(items, rule=DMS_Level3_Base, **kwargs):
    """Creat an association from a list

    Parameters
    ----------
    items: [object [, ...]]
        List of items to add.

    rule: `Association` rule
        The association rule to use.

    kwargs: dict
        Other named parameters required or pertinent to adding
        the items to the association.

    Returns
    -------
    association: `Association`-based instance
        The association with the items added.

    Notes
    -----
    This is a lower-level tool for artificially creating
    an association. As such, the association created may not be valid.
    It is presume the user knows what they are doing.
    """

    asn = rule()
    asn._add_items(items, **kwargs)
    return asn
