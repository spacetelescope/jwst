import os
import sys
import json
import warnings
import os.path as op

from ..lib import suffix
from . import Association, AssociationNotValidError

#-----------------------------------------------------------------------
# Externally callable functions

class AsnFileWarning(Warning):
    pass

def add(asn, filenames, exptype):
    """
    Add new filenames to association

    Parameters
    ----------
    asn : object
        An association object

    filenames : list of str
        The filenames to be added to the association

    exptype : str
        The exposure type of the filenames to be added


    Returns
    -------
    The modified association object


    Raises
    ------
    ValueError
        If a filename to be added was not found in the filesystem


     """
    for filename in filenames:
        path = _path(filename)
        expname = op.basename(path)
        member = {'expname': expname, 'exptype': exptype}

        for product in asn['products']:
            product['members'].append(member)

    return asn


def reader(association_file):
    """
    Read the association file or die trying

    Parameters
    ----------
    association_file : str
        An association filename


    Returns
    -------
    The association object


    Raises
    ------
    IOError
        An error occurred when reading the association file


    """
    association_path = _path(association_file)
    asn_format = association_path.split('.')[-1]
    if asn_format != 'json':
        raise IOError("This is not an association file: " +
                       association_file)
    try:
        with open(association_path) as fd:
            serialized = fd.read()
            asn = Association.load(serialized, format=asn_format)
    except AssociationNotValidError:
        raise IOError("Cannot read association file: " +
                       association_file)
    return asn


def remove(asn, filenames, ignore):
    """
    Remove the named files from the association

    Parameters
    ----------
    asn : object
        An association object

    filenames : list of str
        The filenames to be removed the association

    ignore : bool
        Ignore the filename suffix when matching filenames?


    Returns
    -------
    The modified association object


    Raises
    ------
    ValueError
        If a filename to be removed was not found in the association


    """
    not_found = []
    for filename in filenames:
        found = _lookup(asn, filename, ignore_suffix=ignore)

        if len(found) == 0:
            not_found.append(filename)
        else:
            for i, j in found[::-1]:
                del asn['products'][i]['members'][j]
                if len(asn['products'][i]['members']) == 0:
                    del asn['products'][i]


    if len(not_found) > 0:
        errmsg = "Filenames not found: " + ','.join(not_found)
        warnings.warn(errmsg, AsnFileWarning)

    return asn


def writer(asn, output_file):
    """
    Write the association out to disk

    Parameters
    ----------
    asn : object
        An association object

    output_file : str
        The filename of the association

    Raises
    ------
    ValueError
        The output filename was not a json file


    """

    asn_format = output_file.split('.')[-1]
    if asn_format != 'json':
        raise ValueError('Only json format supported for output: ' +
                         output_file)

    serialized = json.dumps(asn, indent=4, separators=(',', ': '))

    in_place = op.exists(output_file)
    if in_place:
        temp_file = _rename(output_file)

    try:
        fd = open(output_file, 'w')
        fd.write(serialized)
        fd.close()
    except:
        if in_place:
            os.rename(temp_file, output_file)
        raise
    if in_place:
        os.remove(temp_file)

#-----------------------------------------------------------------------
# Internal functions

def _lookup(asn, filename, ignore_suffix=False):
    """
    Look up the locations a file is found in an association
    """

    found = []
    path = _path(filename)
    basename = op.basename(path)
    (root, ext) = op.splitext(basename)
    if ignore_suffix:
        search_key, separator = suffix.remove_suffix(root)
    else:
        search_key = basename

    for i, product in enumerate(asn['products']):
        for j, member in  enumerate(product['members']):
            expname = member.get('expname')
            if expname:
                if ignore_suffix:
                    (root, ext) = op.splitext(expname)
                    match_key, separator = suffix.remove_suffix(root)
                else:
                    match_key = expname
                if search_key == match_key:
                    found.append((i, j))

    return found


def _path(filename):
    """
    Command line filename processing
    """

    return op.abspath(op.expanduser(op.expandvars(filename)))


def _rename(output_file):
    """
    Rename output file to prevent overwriting
    """

    trial = 0
    while 1:
        trial += 1
        temp_file = "%s.sv%02d" % (output_file, trial)
        if not op.exists(temp_file):
            os.rename(output_file, temp_file)
            return temp_file


