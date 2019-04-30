from glob import glob as _sys_glob
import os.path as op
import os
import sys
import pytest
import requests

from ci_watson.artifactory_helpers import (
    BigdataError,
    check_url,
    get_bigdata,
    get_bigdata_root,
)
from .compare_outputs import compare_outputs

from jwst.associations import load_asn

__all__ = [
    'BaseJWSTTest',
]

# Define location of default Artifactory API key, for Jenkins use only
ARTIFACTORY_API_KEY_FILE = '/eng/ssb2/keys/svc_rodata.key'


@pytest.mark.usefixtures('_jail')
@pytest.mark.bigdata
class BaseJWSTTest:
    '''
    Base test class from which to derive JWST regression tests
    '''
    rtol = 0.00001
    atol = 0

    input_loc = ''  # root directory for 'input' files
    ref_loc = []    # root path for 'truth' files: ['test1','truth'] or ['test3']

    ignore_table_keywords = []
    ignore_fields = []
    ignore_hdus = ['ASDF']
    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX', 'FILENAME']

    env = 'dev'

    @pytest.fixture(autouse=True)
    def config_access(self, pytestconfig):
        self.inputs_root = pytestconfig.getini('inputs_root')[0]
        self.results_root = pytestconfig.getini('results_root')[0]

    @property
    def repo_path(self):
        return [self.inputs_root, self.env, self.input_loc]

    def get_data(self, *pathargs, docopy=True):
        """
        Download `filename` into working directory using
        `artifactory_helpers/get_bigdata()`.
        This will then return the full path to the local copy of the file.
        """
        local_file = get_bigdata(*self.repo_path, *pathargs, docopy=docopy)
        return local_file

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

        input_path = [self.inputs_root, self.env, self.input_loc, *self.ref_loc]

        return compare_outputs(outputs,
                               input_path=input_path,
                               docopy=True,
                               results_root=self.results_root,
                               **compare_kws)

    def data_glob(self, *pathargs, glob='*'):
        """Retrieve file list matching glob

        Parameters
        ----------
        pathargs: (str[, ...])
            Path components

        glob: str
            The file name match criterion

        Returns
        -------
        file_paths: [str[, ...]]
            File paths that match the glob criterion.
            Note that the TEST_BIGDATA and `repo_path`
            roots are removed so these results can be fed
            back into `get_data()`
        """

        # Get full path and proceed depending on whether
        # is a local path or URL.
        root = op.join(get_bigdata_root(), *self.repo_path)
        path = op.join(root, *pathargs)
        if op.exists(path):
            file_paths = _data_glob_local(path, glob)
        elif check_url(path):
            file_paths = _data_glob_url(path, glob)
        else:
            raise BigdataError('Path cannot be found: {}'.format(path))

        # Remove the root from the paths
        root_len = len(root) + 1  # +1 to account for the folder delimiter.
        file_paths = [
            file_path[root_len:]
            for file_path in file_paths
        ]
        return file_paths


# Pytest function to support the parameterization of BaseJWSTTestSteps
def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    idlist = [funcargs['id'] for funcargs in funcarglist]
    del argnames[argnames.index('id')]
    metafunc.parametrize(argnames, [[funcargs[name] for name in argnames]
            for funcargs in funcarglist], ids=idlist)


class BaseJWSTTestSteps(BaseJWSTTest):

    params = {'test_steps':[dict(input="",
                                 test_dir=None,
                                 step_class=None,
                                 step_pars=dict(),
                                 output_truth="",
                                 output_hdus=[])
                            ]
             }

    def test_steps(self, input, test_dir, step_class, step_pars,
                   output_truth, output_hdus):
        """
        Template method for parameterizing all the tests of JWST pipeline
        processing steps.
        """

        if test_dir is None:
            return

        self.test_dir = test_dir
        self.ref_loc = [self.test_dir, 'truth']

        # can be removed once all truth files have been updated
        self.ignore_keywords += ['FILENAME']

        input_file = self.get_data(self.test_dir, input)
        result = step_class.call(input_file, save_results=True, **step_pars)

        output_file = result.meta.filename
        result.close()

        output_pars = None
        if isinstance(output_truth, tuple):
            output_pars = output_truth[1]
            output_truth = output_truth[0]

        if not output_pars:
            if output_hdus:
                output_spec = (output_file, output_truth, output_hdus)
            else:
                output_spec = (output_file, output_truth)
        else:
            output_spec = {'files':(output_file, output_truth),
                           'pars':output_pars}
        outputs = [output_spec]
        self.compare_outputs(outputs)


def raw_from_asn(asn_file):
    """
    Return a list of all input files from a given association.

    Parameters
    ----------
    asn_file : str
        Filename for the ASN file.

    Returns
    -------
    members : list of str
        A list of all input files in the association

    """

    members = []
    with open(asn_file) as f:
        asn = load_asn(f)

    for product in asn['products']:
        for member in product['members']:
            members.append(member['expname'])

    return members


def _data_glob_local(path, glob='*'):
    """Perform a glob on the local path

    Parameters
    ----------
    path: File path-like object
        The path to check.

    glob: str
        The file name match criterion

    Returns
    -------
    file_paths: [str[, ...]]
        Full file paths that match the glob criterion
    """
    full_glob = op.join(path, glob)
    return _sys_glob(full_glob)


def _data_glob_url(url, glob='*'):
    """
    Parameters
    ----------
    url: str
        The URL to check

    glob: str
        The file name match criterion

    Returns
    -------
    url_paths: [str[, ...]]
        Full URLS that match the glob criterion
    """
    try:
        envkey = os.environ['API_KEY_FILE']
    except KeyError:
        envkey = ARTIFACTORY_API_KEY_FILE

    try:
        with open(envkey) as fp:
            headers = {'X-JFrog-Art-Api': fp.readline().strip()}
    except (PermissionError, FileNotFoundError):
        print("Warning: Anonymous Artifactory search requests are limited to "
            "1000 results. Use an API key and define API_KEY_FILE environment "
            "variable to get full search results.", file=sys.stderr)
        headers = None

    search_url = op.join(get_bigdata_root(), 'api/search/artifact')
    # Pick out "jwst-pipeline", the repo name
    repo = url.split('/')[4]
    params = {'name': glob, 'repos': repo}
    with requests.get(search_url, params=params, headers=headers) as r:
        url_paths = [a['uri'].replace('api/storage/', '') for a in r.json()['results']]

    return url_paths
