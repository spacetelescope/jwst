"""Regression tests for MIRI MRS PFPC corrections."""

import os

import pytest

# from jwst.regtest.regtestdata import RTData
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]

# Define artifactory source and truth
INPUT_DATA_PATH = "miri/mrs"
TRUTH_PATH = "truth/test_miri_mrs_pfpc"
ASN_FILE = "jw01206-o001_mod_spec3_asn.json"
ASN_DESCRIPTOR = {
    "file_name": ASN_FILE,
    "path": INPUT_DATA_PATH,
    "from_mast": True,
    "asn_files_from_mast": True,
    "mod_code": "_trim_asn",
    "asn_files": [
        "jw01206001001_02103_00001_mirifushort_cal.fits",
        "jw01206001001_02103_00002_mirifushort_cal.fits",
        "jw01206001001_02103_00003_mirifushort_cal.fits",
        "jw01206001001_02103_00004_mirifushort_cal.fits",
        "jw01206001001_02105_00001_mirifulong_cal.fits",
        "jw01206001001_02105_00002_mirifulong_cal.fits",
        "jw01206001001_02105_00003_mirifulong_cal.fits",
        "jw01206001001_02105_00004_mirifulong_cal.fits",
    ],
}

# When RTData changes are merged (#10729), set INPUT_DATA:
# INPUT_DATA = {ASN_FILE: RTData(**ASN_DESCRIPTOR)}


@pytest.fixture(scope="module")
def run_spec3_pfpc(rtdata_module):
    """Run the Spec3Pipeline on association with PFPC."""
    rtdata = rtdata_module

    # todo: should this override reference file be noted in input data?
    rtdata.get_data(f"{INPUT_DATA_PATH}/jwst_miri_pfpc_0000.fits")

    # When RTData changes are merged (#10729), get_asn can use INPUT_DATA directly:
    # rtdata.get_asn(INPUT_DATA[ASN_FILE])

    # For now, get ASN from the filename
    rtdata.get_asn(f"{INPUT_DATA_PATH}/{ASN_FILE}")

    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.pfpc.skip=false",
        "--steps.pfpc.save_results=true",
        "--steps.pfpc.override_pfpc=jwst_miri_pfpc_0000.fits",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.parametrize("band", ["ch1-medium", "ch2-medium", "ch3-short", "ch4-short"])
@pytest.mark.parametrize("suffix", ["x1d", "pfpc"])
def test_miri_mrs_spec3_pfpc(run_spec3_pfpc, fitsdiff_default_kwargs, band, suffix):
    """Regression test for spec3 with PFPC correction."""

    rtdata = run_spec3_pfpc

    output = f"jw01206-o001_t001_miri_{band}_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    # ignore the pfpc file keyword since it will contain a
    # full path for the override file
    fitsdiff_default_kwargs["ignore_keywords"].append("R_PFPC")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def _trim_asn():
    """
    Trim the PID 1206 obs1 ASN file from MAST for faster processing.

    Keeps only the dithers that include ch1 medium and ch3 short,
    so spectral leak will be applied.

    Assumes jw01206-o001_mast_spec3_asn.json exists in the current directory.
    Writes to jw01206-o001_mod_spec3_asn.json.
    """
    import json
    from pathlib import Path

    with Path(ASN_FILE.replace("mod", "mast")).open() as fh:
        asn = json.load(fh)

    # keep only ch1/2 medium and ch3/4 short so spectral leak will work
    members = []
    for member in asn["products"][0]["members"]:
        expname = member["expname"]
        if "2103" in expname and "short" in expname:
            members.append(member)
        elif "2105" in expname and "long" in expname:
            members.append(member)
    asn["products"][0]["members"] = members

    with Path(ASN_FILE).open("w") as fh:
        json.dump(asn, fh, indent=4)
