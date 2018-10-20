import pytest

from ci_watson.artifactory_helpers import (
    check_url,
    get_bigdata_root,
    get_bigdata,
    compare_outputs,
)


@pytest.mark.usefixtures('_jail')
@pytest.mark.bigdata
class BaseJWSTTest:
    '''Base test class from which to derive JWST regression tests
    '''
    rtol = 0.00001
    atol = 0

    input_loc = ''  # root directory for 'input' files
    ref_loc = []    # root path for 'truth' files: ['test1','truth'] or ['test3']

    ignore_table_keywords = []
    ignore_fields = []
    ignore_hdus = ['ASDF']
    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX', 'FILENAME']

    input_repo = 'jwst-pipeline'
    results_root = 'jwst-pipeline-results'
    env = 'dev'

    bigdata_root = get_bigdata_root()
    if bigdata_root and check_url(bigdata_root):
        docopy = True
    else:
        docopy = False

    @property
    def repo_path(self):
        return [self.input_repo, self.env, self.input_loc]

    def get_data(self, *pathargs, docopy=True):
        """
        Download `filename` into working directory using
        `artifactory_helpers/get_bigdata()`.
        This will then return the full path to the local copy of the file.
        """
        # If user has specified action for no_copy, apply it with
        # default behavior being whatever was defined in the base class.
        repo_path = [
            self.input_repo,
            self.env,
            self.input_loc
        ]
        local_file = get_bigdata(*self.repo_path, *pathargs, docopy=self.docopy)

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

        input_path = [self.input_repo, self.env, self.input_loc, *self.ref_loc]

        return compare_outputs(outputs, raise_error=True,
                               input_path=input_path,
                               docopy=self.docopy,
                               results_root=self.results_root,
                               **compare_kws)


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

        result = step_class.call(input_file, **step_pars)

        output_file = result.meta.filename
        result.save(output_file)
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
    Return a list of all MEMBER input files from a given ASN.

    Parameters
    ----------
    asn_file : str
        Filename for the ASN file.

    Returns
    -------
    raw_files : list of str
        A list of input files to process.

    """
    import json

    raw_files = []
    tab = json.load(open(asn_file))

    for product in tab['products']:
        for member in product['members']:
            raw_files.append(member['expname'])

    return raw_files
