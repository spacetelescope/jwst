"""Test module pointing_summary"""

import sys

from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.diff import report_diff_values
from stdatamodels.jwst.datamodels import ImageModel

import jwst.lib.pointing_summary as ps


class TestCalcPointing:
    def setup_class(self):
        """Create model with needed header parameters"""
        model = ImageModel()
        model.meta.target.ra = 90.75541667
        model.meta.target.dec = -66.56055556
        model.meta.pointing.ra_v1 = 91.08142005
        model.meta.pointing.dec_v1 = -66.60547869
        model.meta.wcsinfo.ra_ref = 90.70377653
        model.meta.wcsinfo.dec_ref = -66.59540224
        self.model = model

    def test_calc_pointing_deltas(self):
        """Test calc_pointing_deltas basic running"""
        truth = (
            "Delta(target=<SkyCoord (ICRS): (ra, dec) in deg"
            "\n    (90.75541667, -66.56055556)>, v1=<SkyCoord (ICRS): (ra, dec) in deg"
            "\n    (91.08142005, -66.60547869)>, refpoint=<SkyCoord (ICRS): (ra, dec) in deg"
            "\n    (90.70377653, -66.59540224)>, delta_v1=<Angle 0.13712727"
            " deg>, delta_refpoint=<Angle 0.04044315 deg>)"
        )
        deltas = ps.calc_pointing_deltas(self.model)
        assert truth == str(deltas)

    def test_calc_deltas(self):
        """Test calc_deltas basic running"""
        deltas = ps.calc_deltas([self.model])

        truth = Table.read(
            get_pkg_data_filename("data/calc_deltas_truth.ecsv", package="jwst.lib.tests")
        )

        # exposure value does not matter
        truth[0][0] = "<ImageModel>"

        # round the delta values to a reasonable level
        deltas[0][4] = round(deltas[0][4], 8)
        deltas[0][5] = round(deltas[0][5], 8)
        truth[0][4] = round(truth[0][4], 8)
        truth[0][5] = round(truth[0][5], 8)

        assert report_diff_values(truth, deltas, fileobj=sys.stderr)
