from pathlib import Path

from stdatamodels.jwst import datamodels

from jwst.ami.tests.helpers import example_model
from jwst.pipeline import Ami3Pipeline


def test_niriss_ami3(tmp_cwd):
    """Smoke test for ami3 for NIRISS."""
    input_model = example_model()
    input_name = "test_sci_cal.fits"
    input_model.save(input_name)
    input_model.close()

    psf_model = example_model()
    psf_name = "test_psf_cal.fits"
    psf_model.save(psf_name)
    psf_model.close()

    ami3_asn = {
        "asn_id": "c1000",
        "asn_pool": "test_pool.csv",
        "products": [
            {
                "name": "test_ami3",
                "members": [
                    {"expname": input_name, "exptype": "science"},
                    {"expname": psf_name, "exptype": "psf"},
                ],
            }
        ],
    }

    Ami3Pipeline.call(ami3_asn)

    # Check for expected output files
    expected = [
        "test_sci_c1000_ami-oi.fits",
        "test_psf_c1000_psf-ami-oi.fits",
        "test_ami3_aminorm-oi.fits",
    ]
    for filename in expected:
        assert Path(filename).exists()
        with datamodels.open(filename) as model:
            assert model.meta.cal_step.ami_analyze == "COMPLETE"
            if "norm" in filename:
                assert model.meta.cal_step.ami_normalize == "COMPLETE"
