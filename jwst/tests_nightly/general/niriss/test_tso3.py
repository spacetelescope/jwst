import pytest
from jwst.pipeline import Tso3Pipeline

from jwst.tests.base_classes import BaseJWSTTest, raw_from_asn


@pytest.mark.bigdata
class TestTso3Pipeline(BaseJWSTTest):
    input_loc = 'niriss'
    ref_loc = ['test_caltso3', 'truth']
    test_dir = 'test_caltso3'

    def test_tso3_pipeline_nis(self):
        """Regression test of calwebb_tso3 on NIRISS SOSS simulated data.
        """
        asn_file = self.get_data(self.test_dir,
                                "jw87600-a3001_20170527t111213_tso3_001_asn.json")
        for file in raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        Tso3Pipeline.call(asn_file)

        outputs = [
            # Compare level-2c product
            ('jw87600024001_02101_00001_nis_a3001_crfints.fits',
             'jw87600-a3001_t1_niriss_clear-gr700xd_crfints_ref.fits',
             ['primary', 'sci', 'dq', 'err']),

            # Compare level-3 product
            ('jw87600-a3001_t1_niriss_clear-gr700xd_x1dints.fits',
             'jw87600-a3001_t1_niriss_clear-gr700xd_x1dints_ref.fits',
             ['primary', 'extract1d']),
            ('jw87600-a3001_t1_niriss_clear-gr700xd_x1dints.ecsv',
             'jw87600-a3001_t1_niriss_clear-gr700xd_x1dints_ref.ecsv'),
        ]
        self.compare_outputs(outputs)
