import logging
import pytest
import subprocess

from jwst.stpipe import Step
from jwst.assign_wcs.util import NoDataOnDetectorError
from jwst.pipeline import Spec2Pipeline
from jwst.tests.helpers import LogWatcher


@pytest.mark.bigdata
def test_nirspec_missing_msa_fail(rtdata, fitsdiff_default_kwargs, monkeypatch):
    """
    Test of calwebb_spec2 pipeline performed on NIRSpec MSA exposure
    that's missing an MSAMETFL. Exception should be raised.
    """

    # Get the input file, don't get the MSA file
    rtdata.get_data("nirspec/mos/jw01180025001_05101_00001_nrs2_rate.fits")

    # Run the calwebb_spec2 pipeline
    args = ["calwebb_spec2", rtdata.input]

    # watch for an error log message, we don't use caplog here because
    # something in the test suite messes up the logging during some runs
    # (probably stpipe or the association generator) and causes caplog
    # to sometimes miss the message
    watcher = LogWatcher("Missing MSA meta (MSAMETFL) file")
    monkeypatch.setattr(logging.getLogger("jwst.assign_wcs.nirspec"), "error", watcher)

    with pytest.raises(Exception):
        Step.from_cmdline(args)

    watcher.assert_seen()


@pytest.mark.bigdata
def test_nirspec_missing_msa_nofail(rtdata, fitsdiff_default_kwargs, monkeypatch):
    """
    Test of calwebb_spec2 pipeline performed on NIRSpec MSA exposure
    that's missing an MSAMETFL. Exception should NOT be raised.
    """

    # Get the input file, don't get the MSA file
    rtdata.get_data("nirspec/mos/jw01180025001_05101_00001_nrs2_rate.fits")

    # Run the calwebb_spec2 pipeline
    args = ["calwebb_spec2", rtdata.input, "--fail_on_exception=False"]

    # watch for an error log message, we don't use caplog here because
    # something in the test suite messes up the logging during some runs
    # (probably stpipe or the association generator) and causes caplog
    # to sometimes miss the message
    watcher = LogWatcher("Missing MSA meta (MSAMETFL) file")
    monkeypatch.setattr(logging.getLogger("jwst.assign_wcs.nirspec"), "error", watcher)

    Step.from_cmdline(args)

    watcher.assert_seen()


@pytest.mark.bigdata
def test_nirspec_assignwcs_skip(rtdata, fitsdiff_default_kwargs, caplog):
    """
    Test of calwebb_spec2 pipeline performed on NIRSpec MSA exposure
    with the AssignWcs step skipped. The pipeline should abort.
    """

    # Get the input file
    rtdata.get_data("nirspec/mos/jw01180025001_05101_00001_nrs2_rate.fits")

    # Run the calwebb_spec2 pipeline
    args = ["calwebb_spec2", rtdata.input, "--steps.assign_wcs.skip=True"]

    Step.from_cmdline(args)

    assert "Aborting remaining processing for this exposure." in caplog.text


@pytest.mark.bigdata
def test_nirspec_nrs2_nodata_api(rtdata, fitsdiff_default_kwargs):
    """
    Test of calwebb_spec2 pipeline performed on NIRSpec IFU exposure
    that has a filter/grating combination that produces no data on
    the NRS2 detector. Pipeline should raise an exception.
    """

    # Get the input file
    rtdata.get_data("nirspec/ifu/jw01128009001_0310c_00004_nrs2_rate.fits")

    # Call the Spec2Pipeline
    step = Spec2Pipeline()
    step.assign_wcs.skip = False

    with pytest.raises(NoDataOnDetectorError):
        step.run(rtdata.input)


@pytest.mark.bigdata
def test_nirspec_nrs2_nodata_strun(rtdata, fitsdiff_default_kwargs, caplog):
    """
    Test of calwebb_spec2 pipeline performed on NIRSpec IFU exposure
    that has a filter/grating combination that produces no data on
    the NRS2 detector. Pipeline should return with non-zero exit status.
    """

    # Get the input file
    rtdata.get_data("nirspec/ifu/jw01128009001_0310c_00004_nrs2_rate.fits")

    # Call the Spec2Pipeline
    cmd = ["strun", "jwst.pipeline.Spec2Pipeline", rtdata.input]

    status = subprocess.run(cmd)

    assert status.returncode == 64
