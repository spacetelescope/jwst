import json

import pytest
from stdatamodels.jwst import datamodels

from jwst.ami.tests.helpers import example_model
from jwst.pipeline import Ami3Pipeline


@pytest.fixture()
def ami3_asn(tmp_path):
    input_model = example_model()
    input_name = "test_sci_cal.fits"
    input_model.save(str(tmp_path / input_name))
    input_model.close()

    psf_model = example_model()
    psf_name = "test_psf_cal.fits"
    psf_model.save(str(tmp_path / psf_name))
    psf_model.close()

    asn_path = tmp_path / "test_ami3_asn.json"
    asn = {
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

    with asn_path.open("w") as asn_file:
        json.dump(asn, asn_file)

    return str(asn_path)


def test_niriss_ami3(tmp_path, ami3_asn):
    """Smoke test for ami3 for NIRISS."""
    Ami3Pipeline.call(ami3_asn, output_dir=str(tmp_path))

    # Check for expected output files
    expected = [
        "test_sci_c1000_ami-oi.fits",
        "test_psf_c1000_psf-ami-oi.fits",
        "test_ami3_aminorm-oi.fits",
    ]
    for filename in expected:
        assert (tmp_path / filename).exists()
        with datamodels.open(str(tmp_path / filename)) as model:
            assert model.meta.cal_step.ami_analyze == "COMPLETE"
            if "norm" in filename:
                assert model.meta.cal_step.ami_normalize == "COMPLETE"
