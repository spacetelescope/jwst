#
#  Create the set of JWST Level 1b FITS data
#
#
#  Base directory is /user/rij/jwst/generated_data
#
#  Walk the directory trees below this, looking for proposal files
#
#  For each proposal file
#
#    get the exposures
#    create the level1b data

import os

from . import proposalparser
from . import create_dms_data


def get_proposals(base_directory='.'):
    '''Returns a list of all proposals (files that end in .prop)
    that are below the baseDirectory.  Returns a list of
    (directory, proposalfilename) tuples'''

    proposal_list = []
    for dirname, subdirlist, filelist in os.walk(base_directory):
        for name in filelist:
            if name.endswith('.prop') and not name.startswith('.'):
                proposal_list.append((dirname, name))
    return proposal_list

def write_observation_identifiers(id):
    '''Write out an ObservationIdentifiers.dat file'''
    filename = 'ObservationIdentifiers.dat'
    #
    #  Delete file if it exists
    #
    try:
        os.remove(filename)
    except:
        pass
##     try:
    f1 = open(filename, 'w')
    lines = []
    lines.append('<<file obsid>>\n')
    lines.append('<<header primary>>\n')
    lines.append("PROGRAM   = '%s'\n" % id[0:5])
    lines.append("OBSERVTN  = '%s'\n" % id[5:8])
    lines.append("VISIT     = '%s'\n" % id[8:11])
    lines.append("VISITGRP  = '%s'\n" % id[11:13])
    lines.append("SEQ_ID    = '%s'\n" % id[13:14])
    lines.append("ACT_ID    = '%s'\n" % id[14:16])
    lines.append("EXPOSURE  = '%s'\n" % id[16:21])
    f1.writelines(lines)
    f1.close()
##     except:
##         print("Problem writing to file %s" % filename)
    return filename

def remove_observation_identifiers(obsid):
    '''Remove the observation identifiers file'''
    try:
        os.remove(obsid)
    except:
        pass
    return

def pre_clean():
    #
    #  Clean up the current directory by moving existing Level 1b files
    #  to a subdirectory named 'previous'
    #
    #  Start by emptying the contents of this directory
    try:
        os.chdir('previous')
        previous_files = os.listdir('.')
        for file in previous_files:
            try:
                os.remove(file)
            except:
                pass
        os.chdir('..')
    except:
        os.mkdir('previous')
    filelist = os.listdir('.')
    for file in filelist:
        if file.startswith('jw') and file.endswith('.fits'):
            os.rename(file, 'previous/%s' % file)
    return

def run(base_directory='.', level='1b'):
    '''Do it'''

    proposals = get_proposals(base_directory=base_directory)

    for proposal in proposals:

        directory = proposal[0]
        proposal_file = proposal[1]

        os.chdir(directory)
        #
        #  Clean up everything by moving existing output files to a 'previous'
        #  subdirectory
        pre_clean()

        detectors = proposalparser.get_detectors(proposal_file)

        for detector in detectors:
            base = detector['base']
            subarray = detector['subarray']
            exp_type = detector['exp_type']
            if exp_type == '':
                exp_type = 'UNKNOWN'
            obsid = detector['id']
            obsidfile = write_observation_identifiers(obsid)
            subarray_argument = subarray
            print('Creating Level %s data' % level)
            create_dms_data.create_dms(base,
                                       level=level,
                                       parfile=obsidfile,
                                       subarray=subarray_argument,
                                       exp_type=exp_type)
            remove_observation_identifiers(obsidfile)
