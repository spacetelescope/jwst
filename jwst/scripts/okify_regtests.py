#!/usr/bin/env python

"""
Regression test okifying script.

Requires JFrog CLI (https://jfrog.com/getcli/) configured with credentials
that have write access to the jwst-pipeline repository.
"""
import json
import os
import shutil
import subprocess
import tempfile
from argparse import ArgumentParser
from contextlib import contextmanager
from glob import glob
from pathlib import Path

import asdf
import readchar
from colorama import Fore

ARTIFACTORY_REPO = 'jwst-pipeline-results'
SPECFILE_SUFFIX = '_okify.json'
RTDATA_SUFFIX = '_rtdata.asdf'
TERMINAL_WIDTH = shutil.get_terminal_size((80, 20)).columns


def parse_args():
    parser = ArgumentParser(description='Okify regression test results')
    parser.add_argument(
        'build_number',
        help='Jenkins build number for JWST builds',
        metavar='build-number',
    )
    parser.add_argument(
        '--job-name',
        help='Jenkins job name under [RT] (default: JWST)',
        default='JWST',
        metavar='job-name',
    )
    parser.add_argument(
        '--dry-run', action='store_true', help='pass the --dry-run flag to JFrog CLI'
    )

    return parser.parse_args()


def artifactory_copy(specfile, dry_run=False):
    jfrog_args = []

    if dry_run:
        jfrog_args.append('--dry-run')

    args = list(['jfrog', 'rt', 'cp'] + jfrog_args + [f'--spec={specfile}'])
    subprocess.run(args, check=True)


def artifactory_folder_copy(specfile, dry_run=False):
    """Copy a folder after removing target folder"""
    jfrog_args = []
    if dry_run:
        jfrog_args.append('--dry-run')

    # Since two different jfrog operations are required, need to read in
    # the spec to perform the delete.
    with open(specfile) as fh:
        spec = json.load(fh)
    pattern = spec['files'][0]['pattern'] + '/'
    target = spec['files'][0]['target']

    # Remove the target
    folder = target + Path(pattern).stem
    args = ['jfrog', 'rt', 'del', folder, '--quiet=true'] + jfrog_args
    subprocess.run(args, check=True)

    # Copy pattern to parent of target.
    args = ['jfrog', 'rt', 'cp', '--spec', specfile] + jfrog_args
    subprocess.run(args, check=True)


def artifactory_dispatch(okify_op, specfile, dry_run):
    """Perform the indicated artifactory operation

    Parameters
    ----------
    okify_op : str
        The operation to perform:
        - 'file_copy': Copy individual files
        - 'folder_copy': Replace whole folders

    specfile : str
        The full path to the jfrog spec file

    dry_run : bool
        True to just show what would be done
    """
    if okify_op == 'file_copy':
        artifactory_copy(os.path.abspath(specfile), dry_run=dry_run)
    elif okify_op == 'folder_copy':
        artifactory_folder_copy(os.path.abspath(specfile), dry_run=dry_run)
    else:
        raise RuntimeError(f'Unknown artifactory operation: {okify_op}')


def artifactory_get_breadcrumbs(build_number, job_name, suffix):
    """Download specfiles or other breadcrumb from Artifactory associated with
    a build number and return a list of their locations on the local file system

    An example search would be:

    jfrog rt search jwst-pipeline-results/*/*_okify.json --props='build.number=540;build.name=RT :: JWST'
    """
    build_name = f'RT :: {job_name}'

    # Retrieve all the okify specfiles for failed tests.
    args = list(
        ['jfrog', 'rt', 'dl']
        + [f'{ARTIFACTORY_REPO}/*/*{suffix}']
        + [f'--build={build_name}/{build_number}']
        + ['--flat']
    )
    subprocess.run(args, check=True, capture_output=True)

    return sorted(glob(f'*{suffix}'))


def artifactory_get_build_artifacts(build_number, job_name):
    specfiles = artifactory_get_breadcrumbs(build_number, job_name, SPECFILE_SUFFIX)
    asdffiles = artifactory_get_breadcrumbs(build_number, job_name, RTDATA_SUFFIX)

    if len(specfiles) != len(asdffiles):
        raise RuntimeError('Different number of _okify.json and _rtdata.asdf files')

    for a, b in zip(specfiles, asdffiles):
        if a.replace(SPECFILE_SUFFIX, '') != b.replace(RTDATA_SUFFIX, ''):
            raise RuntimeError('The _okify.json and _rtdata.asdf files are not matched')

    return specfiles, asdffiles


@contextmanager
def pushd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def main():
    args = parse_args()

    build = args.build_number
    name = args.job_name

    # Create and chdir to a temporary directory to store specfiles
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f'Downloading test logs to {tmpdir}')
        with pushd(tmpdir):
            # Retrieve all the okify specfiles for failed tests.
            specfiles, asdffiles = artifactory_get_build_artifacts(build, name)

            number_failed_tests = len(specfiles)

            print(f'{number_failed_tests} failed tests to okify')

            for i, (specfile, asdffile) in enumerate(zip(specfiles, asdffiles)):

                # Print traceback and OKify info for this test failure
                with asdf.open(asdffile) as af:
                    okify_op = af.tree['okify_op']
                    traceback = af.tree['traceback']
                    remote_results_path = af.tree['remote_results_path']
                    output = af.tree['output']
                    truth_remote = af.tree['truth_remote']
                    try:
                        test_name = af.tree['test_name']
                    except KeyError:
                        test_name = 'test_name'

                remote_results = os.path.join(
                    remote_results_path, os.path.basename(output)
                )

                test_number = i + 1

                print(
                    f'{Fore.RED}'
                    + (f' {test_name} '.center(TERMINAL_WIDTH, '—'))
                    + f'{Fore.RESET}'
                )
                print(traceback)
                print(f'{Fore.RED}' + ('—' * TERMINAL_WIDTH) + f'{Fore.RESET}')
                print(f'{Fore.GREEN}OK: {remote_results}')
                print(f'--> {truth_remote}{Fore.RESET}')
                print(
                    f'{Fore.RED}'
                    + (
                        f'[ test {test_number} of {number_failed_tests} ]'.center(
                            TERMINAL_WIDTH, '—'
                        )
                    )
                    + f'{Fore.RESET}'
                )

                # Ask if user wants to okify this test
                while True:
                    print(
                        f'{Fore.GREEN}\'o\' to okify{Fore.RESET}, '
                        f'{Fore.CYAN}\'s\' to skip{Fore.RESET}, '
                        f'{Fore.MAGENTA}\'q\' to quit{Fore.RESET}: '
                    )
                    # Get the keyboard character input without pressing return
                    result = readchar.readkey()
                    if result not in ['o', 's', 'q']:
                        print('Unrecognized command, try again')
                    else:
                        break
                if result == 'q':
                    break
                elif result == 's':
                    pass
                else:
                    artifactory_dispatch(
                        okify_op, os.path.abspath(specfile), dry_run=args.dry_run
                    )
                    print('')


if __name__ == '__main__':
    main()
