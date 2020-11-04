"""asn_gather: Copy or Move data that is listed in an association"""
from pathlib import Path
import subprocess

__all__ = ['asn_gather']


def asn_gather(source_asn_path, destination=None, exp_types=None, copy=True,
               recurse=False, member_root=None, shellcmd='rsync -Pur --no-perms --chmod=ugo=rwX'):
    """Copy/Move members of an association from one location to another

    The association is copied into the destination, re-written such that the member
    list points to the new location of the members. If `copy` is `False`, the member
    list will simply point back to the original location of the members.

    If members cannot be copied, warnings will be generated and the member information
    in the copied association will be left untouched.

    Parameters
    ----------
    source_asn_path : str, pathlib.Path
        The association to gather.

    destination : str, pathlib.Path, or None
        The folder to place the association and its members.
        If None, the current working directory is used.

    exp_types : [str[,...]] or None
        List of exposure types to gather.
        If None, all are gathered.

    copy : bool
        Copy the members to the destination. Otherwise, just copy
        the association itself with the members pointing to the original
        location.

    recurse : bool
        If members appear to come from other associations, attempt
        to gather.

    member_root : str, pathlib.Path, or None
        The folder where the members are found.
        If None, the folder path of the requested association is used.

    Returns
    -------
    dest_asn : pathlib.Path
        The association.
    """
    from .load_as_asn import LoadAsAssociation
    from . import Association

    source_asn_path = Path(source_asn_path)
    source_folder = source_asn_path.parent
    if destination is None:
        dest_folder = Path('./')
    else:
        dest_folder = Path(destination)

    # Create the associations
    source_asn = LoadAsAssociation.load(source_asn_path)
    dest_asn = LoadAsAssociation.load(source_asn_path)

    # Create the new association
    dest_asn['products'] = []
    for src_product in source_asn['products']:
        members = [
            {'expname': src_member['expname'], 'exptype': src_member['exptype']}
            for src_member in src_product['members']
            if exp_types is None or src_member['exptype'] in exp_types
        ]
        if members:
            product = {
                'name': src_product['name'],
                'members': members
            }
            dest_asn['products'].append(product)

    if not dest_asn['products']:
        raise RuntimeError('No products could be gathered.')

    # Save new association.
    dest_path = dest_folder / source_asn_path.name
    _, serialized = dest_asn.dump()
    with open(dest_path, 'w') as fh:
        fh.write(serialized)

    # That's all folks
    return dest_path


def from_cmdline():
    """Collect asn_gather arguments from the commandline"""
