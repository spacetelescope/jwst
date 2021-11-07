import pytest
import subprocess

from jwst.stpipe import Step
from jwst.assign_wcs.util import NoDataOnDetectorError
from jwst.pipeline import Spec2Pipeline


@pytest.mark.bigdata
def test_nirspec_missing_msa_fail(_jail, rtdata, fitsdiff_default_kwargs, caplog):
    """
        Test of calwebb_spec2 pipeline performed on NIRSpec MSA exposure
        that's missing an MSAMETFL. Exception should be raised.
    """

    # Get the input file
    rtdata.get_data('nirspec/mos/f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod.fits')

    # Run the calwebb_spec2 pipeline
    args = ["calwebb_spec2", rtdata.input]

    with pytest.raises(Exception):
        Step.from_cmdline(args)

    assert 'Missing MSA meta (MSAMETFL) file' in caplog.text


@pytest.mark.bigdata
def test_nirspec_missing_msa_nofail(_jail, rtdata, fitsdiff_default_kwargs, caplog):
    """
        Test of calwebb_spec2 pipeline performed on NIRSpec MSA exposure
        that's missing an MSAMETFL. Exception should NOT be raised.
    """

    # Get the input file
    rtdata.get_data('nirspec/mos/f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod.fits')

    # Run the calwebb_spec2 pipeline
    args = ["calwebb_spec2",
            rtdata.input,
            '--fail_on_exception=False']

    Step.from_cmdline(args)

    assert 'Missing MSA meta (MSAMETFL) file' in caplog.text


@pytest.mark.bigdata
def test_nirspec_assignwcs_skip(_jail, rtdata, fitsdiff_default_kwargs, caplog):
    """
        Test of calwebb_spec2 pipeline performed on NIRSpec MSA exposure
        with the AssignWcs step skipped. The pipeline should abort.
    """

    # Get the input file
    rtdata.get_data('nirspec/mos/f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod.fits')

    # Run the calwebb_spec2 pipeline
    args = ["calwebb_spec2",
            rtdata.input,
            '--steps.assign_wcs.skip=True']

    Step.from_cmdline(args)

    assert 'Aborting remaining processing for this exposure.' in caplog.text


@pytest.mark.bigdata
def test_nirspec_nrs2_nodata_api(_jail, rtdata, fitsdiff_default_kwargs):
    """
        Test of calwebb_spec2 pipeline performed on NIRSpec IFU exposure
        that has a filter/grating combination that produces no data on
        the NRS2 detector. Pipeline should raise an exception.
    """

    # Get the input file
    rtdata.get_data('nirspec/ifu/jw84700006001_02101_00001_nrs2_rate.fits')

    # Call the Spec2Pipeline
    step = Spec2Pipeline()
    step.assign_wcs.skip = False

    with pytest.raises(NoDataOnDetectorError):
        step.run(rtdata.input)


@pytest.mark.bigdata
def test_nirspec_nrs2_nodata_strun(_jail, rtdata, fitsdiff_default_kwargs, caplog):
    """
        Test of calwebb_spec2 pipeline performed on NIRSpec IFU exposure
        that has a filter/grating combination that produces no data on
        the NRS2 detector. Pipeline should return with non-zero exit status.
    """

    # Get the input file
    rtdata.get_data('nirspec/ifu/jw84700006001_02101_00001_nrs2_rate.fits')

    # Call the Spec2Pipeline
    cmd = [
        'strun',
        'jwst.pipeline.Spec2Pipeline',
        rtdata.input]

    status = subprocess.run(cmd)

    assert status.returncode == 64
