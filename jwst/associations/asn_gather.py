"""asn_gather: Copy data that is listed in an association"""
import logging
from pathlib import Path
import subprocess

__all__ = ['asn_gather']


# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
LogLevels = [logging.WARNING, logging.INFO, logging.DEBUG]


def asn_gather(association, destination=None, exp_types=None, exclude_types=None,
               source_folder=None,
               shellcmd='rsync -urv --no-perms --chmod=ugo=rwX'):
    """Copy members of an association from one location to another

    The association is copied into the destination, re-written such that the member
    list points to the new location of the members.

    Parameters
    ----------
    association : str, pathlib.Path
        The association to gather.

    destination : str, pathlib.Path, or None
        The folder to place the association and its members.
        If None, the current working directory is used.

    exp_types : [str[,...]] or None
        List of exposure types to gather.
        If None, all are gathered.

    exclude_types : [str[,...]] or None
        List of exposure types to exclude.

    source_folder : str or None
       Folder where the members originate from.
       If None, the folder of the association is presumed. 

    shellcmd : str
        The shell command to use to do the copying of the
        individual members.

    Returns
    -------
    dest_asn : pathlib.Path
        The copied association.
    """
    from .load_as_asn import LoadAsAssociation
    from . import Association

    exclude_types = exclude_types if exclude_types is not None else []
    source_asn_path = Path(association)
    if source_folder is None:
        source_folder = source_asn_path.parent
    else:
        source_folder = Path(source_folder)
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
            src_member
            for src_member in src_product['members']
            if src_member['exptype'] not in exclude_types and \
            (exp_types is None or src_member['exptype'] in exp_types)
        ]
        if members:
            product = {
                'name': src_product['name'],
                'members': members
            }
            dest_asn['products'].append(product)

    if not dest_asn['products']:
        raise RuntimeError('No products could be gathered.')

    # Copy the members.
    shellcmd_args = shellcmd.split(' ')
    for product in dest_asn['products']:
        for member in product['members']:
            src_path = Path(member['expname'])
            logger.info(f'*** Copying member {src_path.name}')
            if str(src_path.parent).startswith('.'):
                src_path = source_folder / src_path
            dest_path = dest_folder / src_path.name
            process_args = shellcmd_args + [str(src_path), str(dest_path)]
            logger.debug(f'Shell command in use: {process_args}')
            result = subprocess.run(process_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
            logger.debug(result.stdout.decode())
            logger.info('...done')
            member['expname'] = dest_path.name

    # Save new association.
    dest_path = dest_folder / source_asn_path.name
    _, serialized = dest_asn.dump()
    logger.info(f'Copying the association file itself {dest_path}')
    with open(dest_path, 'w') as fh:
        fh.write(serialized)

    # That's all folks
    return dest_path


def from_cmdline(args=None):
    """Collect asn_gather arguments from the commandline

    Parameters
    ----------
    args : [str[,...]]
        List of arguments to parse

    Returns : dict
        Dict of the arguments and their values.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Gather an association to a new location'
    )

    parser.add_argument(
        'association',
        help='Association to gather.'
    )
    parser.add_argument(
        'destination',
        help='Folder to copy the association to.'
    )
    parser.add_argument(
        '-s', '--source', dest='source_folder',
        help='Folder where the members currently reside. Default is the folder where the association resides.'
    )
    parser.add_argument(
        '-t', '--exp-types', default=None, dest='exp_types',
        action='append',
        help='Exposure types to gather. If not specified, all exposure types are used.'
    )
    parser.add_argument(
        '-x', '--exclude-types', default=None, dest='exclude_types',
        action='append',
        help='Exposure types to exclude.'
    )
    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help='Increase verbosity. Specifying multiple times adds more output.'
    )
    parser.add_argument(
        '-c', '--cmd', dest='shellcmd',
        default='rsync -urv --no-perms --chmod=ugo=rwX',
        help=('Shell command to use to perform the copy. Specify as a single string. '
              'Default: "%(default)s"')
    )

    parsed = parser.parse_args(args)

    # Set output detail.
    level = LogLevels[min(len(LogLevels)-1, parsed.verbose)]
    logger.setLevel(level)

    # That's all folks.
    gather_args = vars(parsed)
    del gather_args['verbose']
    return gather_args
