"""Regression test helpers."""
import os

import pytest
from astropy.io import fits
from astropy.utils.data import conf

from ci_watson.artifactory_helpers import get_bigdata
from ci_watson.artifactory_helpers import compare_outputs

# Base classes for actual tests.
# NOTE: Named in a way so pytest will not pick them up here.
class BaseTest(object):
    prevdir = os.getcwd()
    use_ftp_crds = False
    timeout = 30  # seconds
    tree = ''
    docopy = False

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

    def get_input_path(self):
        """
        Return path within repository to remote source of input data
        """
        path = [self.tree, self.input_loc]

        return path

    def get_data(self, *args, **kwargs):
        """
        Download `filename` into working directory using
        `artifactory_helpers/get_bigdata()`.
        This will then return the full path to the local copy of the file.
        """
        # If user has specified action for no_copy, apply it with
        # default behavior being whatever was defined in the base class.
        docopy = kwargs.get('docopy', self.docopy)
        """
        local_file = get_bigdata(self.tree,
                                 self.input_loc,
                                 *args)
        """
        local_file = get_bigdata(*self.get_input_path(),
                                 *args,
                                 repo=self.input_repo,
                                 docopy=docopy)

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
        # If user has specified action for no_copy, apply it with
        # default behavior being whatever was defined in the base class.
        no_copy = kwargs.get('no_copy', self.no_copy)

        self.get_data(*args, no_copy=no_copy)

    def compare_outputs(self, outputs, raise_error=True, **kwargs):

        # Parse any user-specified kwargs
        ignore_keywords = kwargs.get('ignore_keywords', self.ignore_keywords)
        ignore_hdus = kwargs.get('ignore_hdus', self.ignore_hdus)
        ignore_fields = kwargs.get('ignore_fields', self.ignore_fields)
        rtol = kwargs.get('rtol', self.rtol)
        atol = kwargs.get('atol', self.atol)

        compare_kws = dict(ignore_fields=ignore_fields, ignore_hdus=ignore_hdus,
                        ignore_keywords=ignore_keywords,
                        rtol=rtol, atol=atol)

        input_path = [self.tree, self.input_loc, *self.ref_loc]
        return compare_outputs(outputs, raise_error=True,
                               input_path=input_path,
                               input_loc=self.input_loc,
                               docopy = self.docopy,
                               input_repo = self.input_repo, 
                               results_root = self.results_root,
                               **compare_kws)

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
