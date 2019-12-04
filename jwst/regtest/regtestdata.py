import os
import pprint

import asdf
from asdf import util
from ci_watson.artifactory_helpers import (
    get_bigdata_root,
    get_bigdata,
    BigdataError,
)


class RegtestData:
    """Defines data paths on Artifactory and data retrieval methods"""

    def __init__(self, env="dev", inputs_root="jwst-pipeline",
        results_root="jwst-pipeline-results", docopy=True,
        input=None, input_remote=None, output=None, truth=None,
        truth_remote=None):
        self._env = env
        self._inputs_root = inputs_root
        self._results_root = results_root

        self.docopy = docopy

        self._input = input
        self._input_remote = input_remote
        self._output = output
        self._truth = truth
        self._truth_remote = truth_remote
        self._bigdata_root = get_bigdata_root()

    def __repr__(self):
        return pprint.pformat(
            dict(input=self.input, output=self.output, truth=self.truth,
            input_remote=self.input_remote, truth_remote=self.truth_remote),
            indent=2
        )

    @property
    def input_remote(self):
        if self._input_remote is not None:
            return os.path.join(*self._input_remote)
        else:
            return None

    @input_remote.setter
    def input_remote(self, value):
        self._input_remote = value.split(os.sep)

    @property
    def truth_remote(self):
        if self._truth_remote is not None:
            return os.path.join(*self._truth_remote)
        else:
            return None

    @truth_remote.setter
    def truth_remote(self, value):
        self._truth_remote = value.split(os.sep)

    @property
    def input(self):
        return self._input

    @input.setter
    def input(self, value):
        self._input = os.path.abspath(value)

    @property
    def truth(self):
        return self._truth

    @truth.setter
    def truth(self, value):
        self._truth = os.path.abspath(value)

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, value):
        self._output = os.path.abspath(value)

    @property
    def bigdata_root(self):
        return self._bigdata_root

    @bigdata_root.setter
    def bigdata_root(self, value):
        return NotImplementedError("Set TEST_BIGDATA environment variable "
            "to change this value.")

    # The methods
    def get_data(self, path=None):
        """Copy data from Artifactory remote resource to the CWD

        Updates self.input and self.input_remote upon completion
        """
        if path is None:
            path = self.input_remote
        else:
            self.input_remote = path
        self.input = get_bigdata(self._inputs_root, self._env,
            os.path.dirname(path), os.path.basename(path), docopy=self.docopy)

        return self.input

    def get_truth(self, path=None):
        """Copy truth data from Artifactory remote resource to the CWD/truth

        Updates self.truth and self.truth_remote on completion
        """
        if path is None:
            path = self.truth_remote
        else:
            self.truth_remote = path
        os.makedirs('truth', exist_ok=True)
        os.chdir('truth')
        try:
            self.truth = get_bigdata(self._inputs_root, self._env,
                os.path.dirname(path), os.path.basename(path), docopy=self.docopy)
            self.truth_remote = os.path.join(self._bigdata_root, self._inputs_root, self._env, path)
        except BigdataError:
            os.chdir('..')
            raise
        os.chdir('..')

        return self.truth

    def get_association(self, asn):
        return NotImplemented

    def to_asdf(self):
        tree = eval(self.__repr__())
        af = asdf.AsdfFile(tree=tree)

    def read(self, filename):
        return NotImplemented
