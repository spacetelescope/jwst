# #############################################################
# Setup a base class and instantiate it in order to provide the
# SDP pool lists for test parametrization.
# #############################################################
from jwst.regtest.regtestdata import RegtestData


class SDPPoolsSource(RegtestData):
    """Retreive the SDP pools from the test data source

    These tests are very much tied to using `pytest` and the `ci-watson` plugin.
    In particular, there are references to `pytest.config` that only exist when
    running under pytest. Such references are stubbed out with best defaults used.
    """

    # The following will be set during both plugin instantiation
    # of `BaseJWSTTest`, or explicitly during dynamic parametrization
    # The empty string defaults indicate that data source roots have
    # been defined.
    test_dir = 'associations/sdp'
    ref_loc = [test_dir, 'truth']

    _inputs_root = ''
    _pool_paths = None
    _results_root = ''
    _truth_paths = None

    @property
    def pool_paths(self):
        """Get the association pools"""
        if self._pool_paths is None:
            self._pool_paths = self.data_glob(self.test_dir + '/pools', glob='*.csv')
        return self._pool_paths

    def truth_paths(self, pool):
        """Get the truth associations"""
        truth_pool_path = '/'.join(self.ref_loc) + '/' + pool
        self._truth_paths = self.data_glob(truth_pool_path, glob='*.json')
        return self._truth_paths
