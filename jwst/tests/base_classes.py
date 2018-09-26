"""Regression test helpers."""
import getpass
import os
import sys
from io import StringIO
import shutil
import datetime
from difflib import unified_diff

import pytest
from astropy.io import fits
from astropy.io.fits import FITSDiff, HDUDiff
from astropy.utils.data import conf

from ci_watson.artifactory_helpers import get_bigdata, generate_upload_schema


# Base classes for actual tests.
# NOTE: Named in a way so pytest will not pick them up here.
class BaseTest(object):
    prevdir = os.getcwd()
    use_ftp_crds = False
    timeout = 30  # seconds
    tree = ''
    copy_local = True

    # Numpy default for allclose comparison
    rtol = 1e-7
    atol = 0

    # To be defined by instrument/test
    input_repo = ''  # e.g., 'drizzlepac' or 'jwst-pipeline'
    results_root = None  # e.g., 'drizzlepac-results'
    input_loc = ''  # root directory for 'input' files
    ref_loc = []    # root path for 'truth' files: ['test1','truth'] or ['test3']
    ignore_keywords = []
    ignore_table_keywords = []
    ignore_fields = []
    ignore_hdus = []

    # To be defined by individual test
    subdir = ''

    @pytest.fixture(autouse=True)
    def setup_class(self, tmpdir, envopt):
        """
        Run test in own dir so we can keep results separate from
        other tests.
        """
        # create working directory specified for the test
        if not tmpdir.ensure(self.subdir, dir=True):
            p = tmpdir.mkdir(self.subdir).strpath
        else:
            p = tmpdir.join(self.subdir).strpath
        os.chdir(p)

        # This controls astropy.io.fits timeout
        conf.remote_timeout = self.timeout

        # Update tree to point to correct environment
        if envopt == "None":
            envopt = ''
        self.tree = envopt

        # Configure environment for tests
        self.set_environ()

    def teardown_class(self):
        """Reset path and variables."""
        conf.reset('remote_timeout')
        os.chdir(self.prevdir)
        if self.use_ftp_crds and self.prevref is not None:
            os.environ[self.refstr] = self.prevref

    def set_environ(self):
        """Set environment variables for test environment to be used

           Each class/test will define what needs to be set, if anything.
           For example, control use of Astrometry updates for HST with:
               os.environ['ASTROMETRY_STEP_CONTROL'] = 'OFF'
        """
        pass

    def get_data_path(self, *args):
        """
        Return path to remote source of input data
        """
        path = os.path.join(self.input_repo, self.tree, self.input_loc, *args)

        return path

    def get_data(self, *args, **kwargs):
        """
        Download `filename` into working directory using
        `artifactory_helpers/get_bigdata()`.
        This will then return the full path to the local copy of the file.
        """
        # If user has specified action for copy_local, apply it with
        # default behavior being whatever was defined in the base class.
        copy_local = kwargs.get('copy_local', self.copy_local)
        local_file = get_bigdata(self.tree,
                                 self.input_loc,
                                 *args)
        """
        local_file = get_bigdata(self.tree,
                                 self.input_loc,
                                 *args,
                                 repo=self.input_repo,
                                 copy_local=copy_local)

        """
        return local_file

    def raw_from_asn(self, asn_file):
        """
        Return the list of  member exposures specified in the
        association file.

        .. WARNING::
        This method needs to be defined by each subclass.
        """
        msg = "Sub-class needs to define this method."
        raise NotImplementedError(msg)


    def get_input_file(self, *args, refsep='$', **kwargs):
        """
        Download or copy input file (e.g., RAW) into the working directory.
        Can be differentiated from 'get_data' by overloading with version for
        the test that performs additional steps, such as identifying and
        downloading calibration ref files from CRDS (HST does not do this
        automatically).
        """
        # If user has specified action for copy_local, apply it with
        # default behavior being whatever was defined in the base class.
        copy_local = kwargs.get('copy_local', self.copy_local)

        self.get_data(*args, copy_local=copy_local)

    def compare_outputs(self, outputs, raise_error=True, **kwargs):
        """
        Compare output with "truth" using appropriate
        diff routine; namely,
            ``fitsdiff`` for FITS file comparisons
            ``unified_diff`` for ASCII products.

        Only after all elements of `outputs` have been
        processed will the method report any success or failure, with
        failure of any one comparison *not* preventing the rest of the
        comparisons to be performed.

        Parameters
        ----------
        outputs : list of tuple or dicts
            This list defines what outputs from running the test will be
            compared.  Three distinct types of values as list elements
            are supported::

              - 2-tuple : (test output filename, truth filename)
              - 3-tuple : (test output filename, truth filename, HDU names)
              - dict : {'files':[], 'pars':()}

            If filename contains extension such as '[hdrtab]' (no quotes),
            it will be interpreted as specifying comparison of just that HDU.

        raise_error : bool
            Raise ``AssertionError`` if difference is found.

        kwargs : keyword-value pairs
            These user-specified inputs will use these values, not the
            values from the same-named class attributes, for the
            parameters that control the operation of the diff
            functions used in the comparison; namely, FITSDiff and HDUDiff.
            The currently supported attributes which can be overidden,
            along with type of values accepted, includes::

              - ignore_keywords : list
              - ignore_hdus : list
              - ignore_fields : list
              - rtol : float
              - atol : float

            NOTE:  Setting these values does NOT change the values of
            the class attributes for these parameters.

        Returns
        -------
        report : str
            Report from ``fitsdiff``.
            This is part of error message if ``raise_error=True``.

        Syntax
        ------
        There are multiple use cases for this method, specifically
        related to how `outputs` are defined upon calling this method.
        The specification of the `outputs` can be any combination of the
        following patterns.

        1. 2-tuple inputs
            >>> outputs = [('file1.fits', 'file1_truth.fits')]

            This definition indicates that `file1.fits` should be compared
            as a whole with `file1_truth.fits`.

        2. 2-tuple inputs with extensions
            >>> outputs = [('file1.fits[hdrtab]',
                            'file1_truth.fits[hdrtab]')]

            This definition indicates that only the HDRTAB extension from
            `file1.fits` will be compared to the HDRTAB extension from
            `file1_truth.fits`.

        3.  3-tuple inputs
            >>> outputs = [('file1.fits', 'file1_truth.fits',
                            ['primary','sci','err','groupdq', 'pixeldq'])]

            This definition indicates that only the extensions specified
            in the list as the 3rd element of the tuple should be compared
            between the two files.  This will cause a temporary FITS
            HDUList object comprising only those extensions specified in
            the list to be generated for each file and those HDUList objects
            will then be compared.

        4.  dictionary of inputs and parameters
            >>> outputs = {'files':('file1.fits', 'file1_truth.fits'),
                           'pars':{'ignore_keywords':self.ignore_keywords+['ROOTNAME']}
                          }

            This definition indicates that all keywords defined by self.ignore_keywords
            along with ROOTNAME will be ignored during the comparison between the
            files specified in 'files'.  Any input parameter for FITSDiff
            or HDUDiff can be specified as part of the `pars` dictionary.
            In addition, the input files listed in `files` can also include
            an extension specification, such as '[hdrtab]', to limit the
            comparison to just that extension.

        Example:
        This example from an actual test definition demonstrates
        how multiple input defintions can be used at the same time.::

            outputs = [( # Compare psfstack product
                        'jw99999-a3001_t1_nircam_f140m-maskbar_psfstack.fits',
                        'jw99999-a3001_t1_nircam_f140m-maskbar_psfstack_ref.fits'
                       ),
                       (
                        'jw9999947001_02102_00002_nrcb3_a3001_crfints.fits',
                        'jw9999947001_02102_00002_nrcb3_a3001_crfints_ref.fits'
                       ),
                       {'files':( # Compare i2d product
                                'jw99999-a3001_t1_nircam_f140m-maskbar_i2d.fits',
                                'jw99999-a3001_t1_nircam_f140m-maskbar_i2d_ref.fits'
                         ),
                         'pars': {'ignore_hdus':self.ignore_hdus+['HDRTAB']}
                       },
                       {'files':( # Compare the HDRTAB in the i2d product
                        'jw99999-a3001_t1_nircam_f140m-maskbar_i2d.fits[hdrtab]',
                        'jw99999-a3001_t1_nircam_f140m-maskbar_i2d_ref.fits[hdrtab]'
                       ),
                        'pars': {'ignore_keywords':
                                 self.ignore_keywords+['NAXIS1', 'TFORM*'],
                                 'ignore_fields':self.ignore_keywords}
                       }
                      ]
        .. NOTE::
        Note that each entry in the list gets interpreted and processed
        separately.
        """
        all_okay = True
        creature_report = ''
        # Create instructions for uploading results to artifactory for use
        # as new comparison/truth files
        testpath, testname = os.path.split(os.path.abspath(os.curdir))
        # organize results by day test was run...could replace with git-hash
        whoami = getpass.getuser() or 'nobody'
        dt = datetime.datetime.now().strftime("%d%b%YT")
        ttime = datetime.datetime.now().strftime("%H_%M_%S")
        user_tag = 'NOT_CI_{}_{}'.format(whoami, ttime)
        build_tag = os.environ.get('BUILD_TAG',  user_tag)
        build_suffix = os.environ.get('BUILD_MATRIX_SUFFIX', 'standalone')
        testdir = "{}_{}_{}".format(testname, build_tag, build_suffix)

        # Parse any user-specified kwargs
        ignore_keywords = kwargs.get('ignore_keywords', self.ignore_keywords)
        ignore_hdus = kwargs.get('ignore_hdus', self.ignore_hdus)
        ignore_fields = kwargs.get('ignore_fields', self.ignore_fields)
        rtol = kwargs.get('rtol', self.rtol)
        atol = kwargs.get('atol', self.atol)

        updated_outputs = []
        extn_list = None
        for entry in outputs:
            num_entries = len(entry)
            if isinstance(entry, dict):
                actual = entry['files'][0]
                desired = entry['files'][1]
                diff_pars = entry['pars']
                ignore_keywords = diff_pars.get('ignore_keywords', ignore_keywords)
                ignore_hdus = diff_pars.get('ignore_hdus', ignore_hdus)
                ignore_fields = diff_pars.get('ignore_fields', ignore_fields)
                rtol = diff_pars.get('rtol', rtol)
                atol = diff_pars.get('atol', atol)
            elif num_entries == 2:
                actual, desired = entry
            elif num_entries == 3:
                actual, desired, extn_list = entry

            if desired.endswith(']'):
                desired_name, desired_extn = desired.split('[')
                desired_extn = desired_extn.replace(']','')
            else:
                desired_name = desired
                desired_extn = None

            # Get "truth" image
            s = self.get_data(*self.ref_loc, desired_name)
            if s is not None:
                desired = s
                if desired_extn is not None:
                    desired = "{}[{}]".format(desired, desired_extn)
            print("\nComparing:\n {} \nto\n {}".format(actual, desired))
            if actual.endswith('fits'):
                # Build HDULists for comparison based on user-specified extensions
                if extn_list is not None:
                    actual = build_hdulist(actual, extn_list)
                    desired = build_hdulist(desired, extn_list)

                # Working with FITS files...
                fdiff = FITSDiff(actual, desired, rtol=rtol, atol=atol,
                                 ignore_hdus=ignore_hdus,
                                 ignore_keywords=ignore_keywords)
                creature_report += fdiff.report()
                if not fdiff.identical:
                    # Only keep track of failed results which need to
                    # be used to replace the truth files (if OK).
                    updated_outputs.append((actual, desired))
                if not fdiff.identical and all_okay:
                    all_okay = False
            elif desired_extn is not None:
                # Specific element of FITS file specified
                actual_hdu = get_hdu(actual)
                desired_hdu = get_hdu(desired)

                # Working with FITS Binary table with header...
                fdiff = HDUDiff(actual_hdu, desired_hdu, rtol=rtol, atol=atol,
                                 ignore_keywords=ignore_keywords,
                                 ignore_fields=ignore_fields)
                creature_report += fdiff.report()
                if not fdiff.identical:
                    # Only keep track of failed results which need to
                    # be used to replace the truth files (if OK).
                    updated_outputs.append((actual, desired))
                if not fdiff.identical and all_okay:
                    all_okay = False
            else:
                # ASCII-based diff
                with open(actual) as afile:
                    actual_lines = afile.readlines()
                with open(desired) as dfile:
                    desired_lines = dfile.readlines()
                udiff = unified_diff(actual_lines, desired_lines,
                                     fromfile=actual, tofile=desired)

                old_stdout = sys.stdout
                udiffIO = StringIO()
                sys.stdout = udiffIO
                sys.stdout.writelines(udiff)
                sys.stdout = old_stdout
                udiff_report = udiffIO.getvalue()
                creature_report += udiff_report
                if len(udiff_report) > 2 and all_okay:
                    all_okay = False
                if len(udiff_report) > 2:
                    # Only keep track of failed results which need to
                    # be used to replace the truth files (if OK).
                    updated_outputs.append((actual, desired))

        if not all_okay and self.results_root is not None:
            tree = os.path.join(self.results_root, self.input_loc,
                            dt, testdir) + os.sep
            test_cwd = os.getcwd()
            # Write out JSON file to enable retention of different results
            new_truths = [os.path.join(test_cwd, os.path.basename(i[1])) for i in updated_outputs]

            for actual_files,new_truth in zip(updated_outputs, new_truths):
                print("Renaming {} as new 'truth' file: {}".format(
                      actual_files[0], new_truth))
                shutil.move(actual_files[0], new_truth)
            log_pattern = [os.path.join(os.path.dirname(x), '*.log') for x in new_truths]
            generate_upload_schema(pattern=new_truths + log_pattern,
                           testname=testname,
                           target= tree)

        if not all_okay and raise_error:
            raise AssertionError(os.linesep + creature_report)

        return creature_report

def get_hdu(filename):
    """Return the HDU for the file and extension specified in the filename.

       This routine expects the filename to be of the format:
           <filename>.fits[extn]

        For example, "jw99999-a3001_t1_nircam_f140m-maskbar_i2d.fits[hdrtab]"
    """
    froot, fextn = filename.split('[')
    fextn = fextn.replace(']','')
    fits_file = fits.open(froot)
    return fits_file[fextn]

def build_hdulist(filename, extn_list):
    """Create a new HDUList object based on extensions specified in extn_list"""
    f = fits.open(filename)
    fhdu = [f[extn] for extn in extn_list]

    return fhdu
