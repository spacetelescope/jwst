from glob import glob
import os

import pytest


@pytest.mark.bigdata
def test_regtestdata_get_data(rtdata, _jail):
    rtdata.get_data("infrastructure/test_regtestdata/file1_rate.fits")
    rtdata.output = "file1_cal.fits"

    assert rtdata.input == os.path.join(os.getcwd(), "file1_rate.fits")


@pytest.mark.bigdata
def test_regtestdata_get_truth(rtdata, _jail):
    rtdata.get_truth("infrastructure/test_regtestdata/file1_rate.fits")
    rtdata.output = "file1_rate.fits"

    assert rtdata.truth == os.path.join(os.getcwd(), "truth", "file1_rate.fits")


@pytest.mark.bigdata
def test_regtestdata_get_asn(_jail, rtdata):
    rtdata.get_asn("infrastructure/test_regtestdata/my_asn.json")
    files = glob("*.fits")
    rtdata.output = "file1_rate.fits"

    assert os.path.isfile("my_asn.json")
    assert len(files) == 3


def test_fitsdiff_defaults(fitsdiff_default_kwargs):
    assert 'ASDF' in fitsdiff_default_kwargs['ignore_hdus']
