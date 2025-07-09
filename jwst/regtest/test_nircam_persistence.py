import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_persistence_step(rtdata_module):
    """
    Run persistence step on NIRCAM data.

    Returns
    -------
    rtdata : RegtestData
        Updated RegtestData object with inputs set
    """
    rtdata = rtdata_module

    # Run persistence step on the _linearity file of the previous exposure to get the
    # _trapsfilled file
    # Need custom trap_density file as the default reference file is all zeros
    rtdata.get_data("nircam/persistence/trap_density.fits")
    trapdensity = rtdata.input
    rtdata.get_data("nircam/persistence/jw01076101001_02101_00001_nrca1_linearity.fits")
    args = [
        "jwst.persistence.PersistenceStep",
        rtdata.input,
        "--override_trapdensity=" + trapdensity,
    ]
    Step.from_cmdline(args)

    # Now run the step on the following exposure using the _trapsfilled file
    # created by the above run
    rtdata.get_data("nircam/persistence/jw01076101001_02101_00002_nrca1_linearity.fits")
    args = [
        "jwst.persistence.PersistenceStep",
        rtdata.input,
        "--override_trapdensity=" + trapdensity,
        "--input_trapsfilled=jw01076101001_02101_00001_nrca1_trapsfilled.fits",
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
def test_persistence_step(run_persistence_step, fitsdiff_default_kwargs):
    """Test the output from the persistence step."""
    rtdata = run_persistence_step
    output = "jw01076101001_02101_00002_nrca1_persistencestep.fits"
    rtdata.output = output
    rtdata.get_truth(
        "truth/test_nircam_persistence/jw01076101001_02101_00002_nrca1_persistencestep.fits"
    )

    # Ignore the custom trap density file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_TRPDEN")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
