import copy
from datetime import datetime
import os
from difflib import unified_diff
from io import StringIO

from ci_watson.artifactory_helpers import (
    get_bigdata,
    BigdataError,
    generate_upload_schema,
)

from astropy.io import fits
from astropy.io.fits import FITSDiff, HDUDiff


TODAYS_DATE = datetime.now().strftime("%Y-%m-%d")

def compare_outputs(outputs, raise_error=True, ignore_keywords=[],
                    ignore_hdus=[], ignore_fields=[], rtol=0.0, atol=0.0,
                    input_path=[], docopy=True, results_root=None,
                    verbose=True):
    """
    Compare output with "truth" using appropriate
    diff routine; namely:
    * ``fitsdiff`` for FITS file comparisons.
    * ``unified_diff`` for ASCII products.
    Only after all elements of ``outputs`` have been
    processed will the method report any success or failure, with
    failure of any one comparison *not* preventing the rest of the
    comparisons to be performed.
    Parameters
    ----------
    outputs : list of tuple or dict
        This list defines what outputs from running the test will be
        compared.  Three distinct types of values as list elements
        are supported:
        * 2-tuple : ``(test output filename, truth filename)``
        * 3-tuple : ``(test output filename, truth filename, HDU names)``
        * dict : ``{'files': (output, truth), 'pars': {key: val}}``
        If filename contains extension such as ``[hdrtab]``,
        it will be interpreted as specifying comparison of just that HDU.
    raise_error : bool
        Raise ``AssertionError`` if difference is found.
    ignore_keywords : list of str
        List of FITS header keywords to be ignored by
        ``FITSDiff`` and ``HDUDiff``.
    ignore_hdus : list of str
        List of FITS HDU names to ignore by ``FITSDiff``.
        This is only available for ``astropy>=3.1``.
    ignore_fields : list of str
        List FITS table column names to be ignored by
        ``FITSDiff`` and ``HDUDiff``.
    rtol, atol : float
        Relative and absolute tolerance to be used by
        ``FITSDiff`` and ``HDUDiff``.
    input_path : list or tuple
        A series of sub-directory names under :func:`get_bigdata_root`
        that leads to the path of the 'truth' files to be compared
        against. If not provided, it assumes that 'truth' is in the
        working directory. For example, with :func:`get_bigdata_root`
        pointing to ``/grp/test_data``, a file at::
            /grp/test_data/pipeline/dev/ins/test_1/test_a.py
        would require ``input_path`` of::
            ["pipeline", "dev", "ins", "test_1"]
    docopy : bool
        If `True`, 'truth' will be copied to output directory before
        comparison is done.
    results_root : str or `None`
        If not `None`, for every failed comparison, the test output
        is automatically renamed to the given 'truth' in the output
        directory and :func:`generate_upload_schema` will be called
        to generate a JSON scheme for Artifactory upload.
        If you do not need this functionality, use ``results_root=None``.
    verbose : bool
        Print extra info to screen.
    Returns
    -------
    creature_report : str
        Report from FITS or ASCII comparator.
        This is part of error message if ``raise_error=True``.
    Examples
    --------
    There are multiple use cases for this method, specifically
    related to how ``outputs`` are defined upon calling this method.
    The specification of the ``outputs`` can be any combination of the
    following patterns:
    1. 2-tuple inputs::
           outputs = [('file1.fits', 'file1_truth.fits')]
       This definition indicates that ``file1.fits`` should be compared
       as a whole with ``file1_truth.fits``.
    2. 2-tuple inputs with extensions::
           outputs = [('file1.fits[hdrtab]', 'file1_truth.fits[hdrtab]')]
       This definition indicates that only the HDRTAB extension from
       ``file1.fits`` will be compared to the HDRTAB extension from
       ``file1_truth.fits``.
    3. 3-tuple inputs::
           outputs = [('file1.fits', 'file1_truth.fits', ['primary', 'sci'])]
       This definition indicates that only the PRIMARY and SCI extensions
       should be compared between the two files. This creates a temporary
       ``HDUList`` object comprising only the given extensions for comparison.
    4. Dictionary of inputs and parameters::
           outputs = [{'files': ('file1.fits', 'file1_truth.fits'),
                       'pars': {'ignore_keywords': ['ROOTNAME']}}]
        This definition indicates that ROOTNAME will be ignored during
        the comparison between the files specified in ``'files'``.
        Any input parameter for ``FITSDiff`` or ``HDUDiff`` can be specified
        as part of the ``'pars'`` dictionary.
        In addition, the input files listed in ``'files'`` can also include
        an extension specification, such as ``[hdrtab]``, to limit the
        comparison to just that extension.
    This example from an actual test definition demonstrates
    how multiple input defintions can be used at the same time::
        outputs = [
            ('jw99999_nircam_f140m-maskbar_psfstack.fits',
             'jw99999_nircam_f140m-maskbar_psfstack_ref.fits'
            ),
            ('jw9999947001_02102_00002_nrcb3_a3001_crfints.fits',
             'jw9999947001_02102_00002_nrcb3_a3001_crfints_ref.fits'
            ),
            {'files': ('jw99999_nircam_f140m-maskbar_i2d.fits',
                       'jw99999_nircam_f140m-maskbar_i2d_ref.fits'),
             'pars': {'ignore_hdus': ['HDRTAB']},
            {'files': ('jw99999_nircam_f140m-maskbar_i2d.fits',
                       'jw99999_nircam_f140m-maskbar_i2d_ref.fits',
                       ['primary','sci','dq']),
             'pars': {'rtol': 0.000001}
            },
            {'files': ('jw99999_nircam_f140m-maskbar_i2d.fits[hdrtab]',
                       'jw99999_nircam_f140m-maskbar_i2d_ref.fits[hdrtab]'),
             'pars': {'ignore_keywords': ['NAXIS1', 'TFORM*'],
                      'ignore_fields': ['COL1', 'COL2']}
            }]
    .. note:: Each ``outputs`` entry in the list gets interpreted and processed
              separately.
    """
    __tracebackhide__ = True
    default_kwargs = {'rtol': rtol, 'atol': atol,
                      'ignore_keywords': ignore_keywords,
                      'ignore_fields': ignore_fields,
                      'ignore_hdus': ignore_hdus}

    all_okay = True
    creature_report = ''
    updated_outputs = []  # To track outputs for Artifactory JSON schema

    for entry in outputs:
        diff_kwargs = copy.deepcopy(default_kwargs)
        extn_list = None
        num_entries = len(entry)

        if isinstance(entry, dict):
            entry_files = entry['files']
            actual = entry_files[0]
            desired = entry_files[1]
            if len(entry_files) > 2:
                extn_list = entry_files[2]
            diff_kwargs.update(entry.get('pars', {}))
        elif num_entries == 2:
            actual, desired = entry
        elif num_entries == 3:
            actual, desired, extn_list = entry
        else:
            all_okay = False
            creature_report += '\nERROR: Cannot handle entry {}\n'.format(
                entry)
            continue

        # TODO: Use regex?
        if actual.endswith(']'):
            if extn_list is not None:
                all_okay = False
                creature_report += (
                    '\nERROR: Ambiguous extension requirements '
                    'for {} ({})\n'.format(actual, extn_list))
                continue
            actual_name, actual_extn = actual.split('[')
            actual_extn = actual_extn.replace(']', '')
        else:
            actual_name = actual
            actual_extn = None

        if desired.endswith(']'):
            if extn_list is not None:
                all_okay = False
                creature_report += (
                    '\nERROR: Ambiguous extension requirements '
                    'for {} ({})\n'.format(desired, extn_list))
                continue
            desired_name, desired_extn = desired.split('[')
            desired_extn = desired_extn.replace(']', '')
        else:
            desired_name = desired
            desired_extn = None

        actual = os.path.abspath(actual)

        # Get "truth" image
        try:
            os.makedirs('truth', exist_ok=True)
            os.chdir('truth')
            desired = get_bigdata(*input_path, desired_name, docopy=docopy)
            desired = os.path.abspath(desired)
            os.chdir('..')
        except BigdataError:
            all_okay = False
            creature_report += '\nERROR: Cannot find {} in {}\n'.format(
                desired_name, input_path)
            continue

        if desired_extn is not None:
            desired_name = desired
            desired = "{}[{}]".format(desired, desired_extn)

        if verbose:
            print("\nComparing:\n {}\n {}".format(actual, desired))

        if actual.endswith('.fits') and desired.endswith('.fits'):
            # Build HDULists for comparison based on user-specified extensions
            if extn_list is not None:
                with fits.open(actual) as f_act:
                    with fits.open(desired) as f_des:
                        actual_hdu = fits.HDUList(
                            [f_act[extn] for extn in extn_list])
                        actual_hdu.filename = lambda: os.path.basename(actual)
                        desired_hdu = fits.HDUList(
                            [f_des[extn] for extn in extn_list])
                        desired_hdu.filename = lambda: os.path.basename(desired)
                        fdiff = FITSDiff(actual_hdu, desired_hdu,
                                         **diff_kwargs)
                        creature_report += '\na: {}\nb: {}\n'.format(
                            actual, desired)  # diff report only gives hash
            # Working with FITS files...
            else:
                fdiff = FITSDiff(actual, desired, **diff_kwargs)

            creature_report += fdiff.report()

            if not fdiff.identical:
                all_okay = False
                # Only keep track of failed results which need to
                # be used to replace the truth files (if OK).
                updated_outputs.append((actual, desired))

        elif actual_extn is not None or desired_extn is not None:
            if 'ignore_hdus' in diff_kwargs:  # pragma: no cover
                diff_kwargs.pop('ignore_hdus')  # Not applicable

            # Specific element of FITS file specified
            with fits.open(actual_name) as f_act:
                with fits.open(desired_name) as f_des:
                    actual_hdu = f_act[actual_extn]
                    desired_hdu = f_des[desired_extn]
                    fdiff = HDUDiff(actual_hdu, desired_hdu, **diff_kwargs)

            creature_report += 'a: {}\nb: {}\n'.format(actual, desired)
            creature_report += fdiff.report()

            if not fdiff.identical:
                all_okay = False
                # Only keep track of failed results which need to
                # be used to replace the truth files (if OK).
                updated_outputs.append((actual_name, desired_name))

        else:
            # ASCII-based diff
            with open(actual) as afile:
                actual_lines = afile.readlines()
            with open(desired) as dfile:
                desired_lines = dfile.readlines()

            udiff = unified_diff(actual_lines, desired_lines,
                                 fromfile=actual, tofile=desired)
            udiffIO = StringIO()
            udiffIO.writelines(udiff)
            udiff_report = udiffIO.getvalue()
            udiffIO.close()

            if len(udiff_report) == 0:
                creature_report += ('\na: {}\nb: {}\nNo differences '
                                    'found.\n'.format(actual, desired))
            else:
                all_okay = False
                creature_report += udiff_report
                # Only keep track of failed results which need to
                # be used to replace the truth files (if OK).
                updated_outputs.append((actual, desired))

    if not all_okay and results_root is not None:  # pragma: no cover
        schema_pattern, tree, testname = generate_upload_params(
            results_root, updated_outputs, verbose=verbose)
        generate_upload_schema(schema_pattern, tree, testname)

    if not all_okay and raise_error:
        raise AssertionError(os.linesep + creature_report)

    return creature_report


def generate_upload_params(results_root, updated_outputs, verbose=True):
    """
    Generate pattern, target, and test name for :func:`generate_upload_schema`.
    This uses ``BUILD_TAG`` and ``BUILD_MATRIX_SUFFIX`` on Jenkins CI to create
    meaningful Artifactory target path. They are optional for local runs.
    Other attributes like user, time stamp, and test name are also
    automatically determined.
    In addition to renamed outputs, ``*.log``is also inserted into the
    ``schema_pattern``.
    Parameters
    ----------
    results_root : str
        See :func:`compare_outputs` for more info.
    updated_outputs : list
        List containing tuples of ``(actual, desired)`` of failed
        test output comparison to be processed.
    verbose : bool
        Print extra info to screen.
    Returns
    -------
    schema_pattern, tree, testname
        Analogous to ``pattern``, ``target``, and ``testname`` that are
        passed into :func:`generate_upload_schema`, respectively.
    """
    import getpass

    # Create instructions for uploading results to artifactory for use
    # as new comparison/truth files
    testname = os.path.split(os.path.abspath(os.curdir))[1]

    # Meaningful test dir from build info.
    # TODO: Organize results by day test was run. Could replace with git-hash
    whoami = getpass.getuser() or 'nobody'
    user_tag = 'NOT_CI_{}'.format(whoami)
    build_tag = os.environ.get('BUILD_TAG', user_tag)
    build_matrix_suffix = os.environ.get('BUILD_MATRIX_SUFFIX', '0')
    subdir = '{}_{}_{}'.format(TODAYS_DATE, build_tag, build_matrix_suffix)
    tree = os.path.join(results_root, subdir, testname) + os.sep
    schema_pattern = []

    # Write out JSON file to enable retention of different results.
    # Also rename outputs as new truths.
    for test_result, truth in updated_outputs:
        schema_pattern.append(test_result)
        if verbose:
            print("\nFailed comparison:")
            print("    {}".format(test_result))
            print("    {}".format(truth))

    return schema_pattern, tree, testname
