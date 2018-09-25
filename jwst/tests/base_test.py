import pytest
#from ci_watson.base_classes import BaseTest
from .base_classes import BaseTest

@pytest.mark.usefixtures('_jail')
class BaseJWSTTest(BaseTest):
    ignore_hdus = ['ASDF']
    ignore_keywords = ['DATE', 'CAL_VER', 'CAL_VCS', 'CRDS_VER', 'CRDS_CTX', 'FILENAME']
    input_repo = 'jwst-pipeline'
    results_root = 'jwst-pipeline-results'

    copy_local = False  # Do not make additional copy by default
    rtol = 0.00001

    def set_environ(self):
        self.tree = ''

    def raw_from_asn(self, asn_file):
        return raw_from_asn(asn_file)


# Pytest function to support the parameterization of these classes
def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    idlist = [funcargs['id'] for funcargs in funcarglist]
    del argnames[argnames.index('id')]
    metafunc.parametrize(argnames, [[funcargs[name] for name in argnames]
            for funcargs in funcarglist], ids=idlist)


@pytest.mark.usefixtures('_jail')
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
        Template method for parameterizing all the tests of MIRI pipeline
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


# function(s) for tests written using functional form, not as BaseTest class tests
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
