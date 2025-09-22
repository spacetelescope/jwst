"""Association edit operations."""

import json
import os
import warnings
from pathlib import Path

from jwst.associations import Association, AssociationNotValidError
from jwst.lib import suffix

__all__ = ["AsnFileWarning", "add", "reader", "remove", "writer"]


class AsnFileWarning(Warning):
    """Generic association file warning."""

    pass


def add(asn, filenames, exptype):
    """
    Add new filenames to association.

    Parameters
    ----------
    asn : object
        An association object.
    filenames : list of str
        The filenames to be added to the association.
    exptype : str
        The exposure type of the filenames to be added.

    Returns
    -------
    object
        The modified association object.
    """
    for filename in filenames:
        path = _path(filename)
        expname = Path(path).name
        member = {"expname": expname, "exptype": exptype}

        for product in asn["products"]:
            product["members"].append(member)

    return asn


def reader(association_file):
    """
    Read the association file.

    Parameters
    ----------
    association_file : str
        An association filename.

    Returns
    -------
    ~jwst.associations.association.Association
        The association object.
    """
    association_path = _path(association_file)
    asn_format = association_path.suffix[1:]
    if asn_format != "json":
        raise OSError("This is not an association file: " + association_file)
    try:
        with association_path.open() as fd:
            serialized = fd.read()
            asn = Association.load(serialized, format=asn_format)
    except AssociationNotValidError as err:
        raise OSError("Cannot read association file: " + association_file) from err
    return asn


def remove(asn, filenames, ignore):
    """
    Remove the named files from the association.

    Parameters
    ----------
    asn : object
        An association object.
    filenames : list of str
        The filenames to be removed from the association.
    ignore : bool
        Ignore the filename suffix when matching filenames if True.

    Returns
    -------
    ~jwst.associations.association.Association
        The modified association object.
    """
    not_found = []
    for filename in filenames:
        found = _lookup(asn, filename, ignore_suffix=ignore)

        if len(found) == 0:
            not_found.append(filename)
        else:
            for i, j in found[::-1]:
                del asn["products"][i]["members"][j]
                if len(asn["products"][i]["members"]) == 0:
                    del asn["products"][i]

    if len(not_found) > 0:
        errmsg = "Filenames not found: " + ",".join(not_found)
        warnings.warn(errmsg, AsnFileWarning, stacklevel=1)

    return asn


def _lookup(asn, filename, ignore_suffix=False):
    """
    Look up the locations where a file is found in an association.

    Parameters
    ----------
    asn : ~jwst.associations.association.Association
        The input association object.
    filename : str or Path
        The filename to find in the association.

    Returns
    -------
    list
        The list of product-member tuple values where the
        provided filename is present in the association.
    """
    found = []
    path = _path(filename)
    basename = path.name
    root = path.stem
    if ignore_suffix:
        search_key, _ = suffix.remove_suffix(root)
    else:
        search_key = basename

    for i, product in enumerate(asn["products"]):
        for j, member in enumerate(product["members"]):
            expname = member.get("expname")
            if expname:
                if ignore_suffix:
                    root = Path(expname).stem
                    match_key, _ = suffix.remove_suffix(root)
                else:
                    match_key = expname
                if search_key == match_key:
                    found.append((i, j))

    return found


def _path(filename):
    """
    Command line filename processing.

    Parameters
    ----------
    filename : str or Path
        The filename.

    Returns
    -------
    str
        Full path including any username and environment
        variables included.
    """
    return Path(os.path.expandvars(filename)).expanduser().resolve()


def writer(asn, output_file):
    """
    Write the association out to disk.

    Parameters
    ----------
    asn : object
        An association object

    output_file : str or Path
        The filename of the association

    Raises
    ------
    ValueError
        The output filename was not a json file
    """
    output_file = Path(output_file)
    if output_file.suffix != ".json":
        raise ValueError(f"Only json format supported for output: {output_file}")

    serialized = json.dumps(asn, indent=4, separators=(",", ": "))

    in_place = Path.exists(output_file)
    if in_place:
        temp_file = _rename(output_file)

    try:
        with Path.open(output_file, "w") as fd:
            fd.write(serialized)
    except Exception:
        if in_place:
            Path.rename(temp_file, output_file)
        raise
    if in_place:
        Path.unlink(temp_file)


def _rename(output_file):
    """
    Rename output file to prevent overwriting.

    Parameters
    ----------
    output_file : Path
        The filename of the association

    Returns
    -------
    Path
        The renamed filename Path.
    """
    trial = 0
    while 1:
        trial += 1
        temp_file = output_file.with_name(f"{output_file.name}.sv{trial}")
        if not Path.exists(temp_file):
            Path.rename(output_file, temp_file)
            return temp_file
