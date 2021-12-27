# #############################################################
# Setup a base class and instantiate it in order to provide the
# SDP pool lists for test parametrization.
# #############################################################
import os

from jwst.regtest.regtestdata import RegtestData


class SDPPoolsSource(RegtestData):
    """Retrieve the SDP pools from the test data source

    This class defines a common relationship between test collection, via
    `jwst.regtest.conftest.pool_path_fixture` and the test using said fixture,
    `jwst.regtest.test_associations_sdp_pools`.
    """

    # Setup some defaults.
    test_dir = 'associations/sdp'
    ref_loc = [test_dir, 'truth']

    _inputs_root = ''
    _pool_paths = None
    _results_root = ''

    def __init__(self, okify_op='folder_copy', **kwargs):
        super().__init__(okify_op=okify_op, **kwargs)

    @property
    def pool_paths(self):
        """Get the association pools"""
        if self._pool_paths is None:
            self._pool_paths = self.data_glob(self.test_dir + '/pools', glob='*.csv')
        return self._pool_paths

    def truth_paths(self, pool):
        """Get the truth associations"""
        paths = []
        truth_pool_path = '/'.join(self.ref_loc) + '/' + pool
        for path in self.data_glob(truth_pool_path, glob='*.json'):
            paths.append(self.get_truth(path))
        self.truth_remote = os.path.join(self._inputs_root, self.env, truth_pool_path)

        return paths
