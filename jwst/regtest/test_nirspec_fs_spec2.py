import os

import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

import jwst.datamodels as dm
from jwst.lib.suffix import replace_suffix
from jwst.pathloss import PathLossStep
from jwst.stpipe import Step

file_roots = [
    'jw02072-o002_20221206t143745_spec2_00001_asn.json',
    'jw01245-o002_20230107t223023_spec2_00001_asn.json',
    'jw01309-o022_20230113t025924_spec2_00001_asn.json'
]

asn_memberdict = {
    'jw01309-o022_20230113t025924_spec2_00001_asn.json':
        ['jw01309022001_04102_00004_nrs2_rate.fits',
         'jw01309022001_04102_00002_nrs2_rate.fits',
         'jw01309022001_04102_00001_nrs2_rate.fits'],
    'jw01245-o002_20230107t223023_spec2_00001_asn.json':
        ['jw01245002001_04102_00002_nrs1_rate.fits',
         'jw01245002001_04102_00001_nrs1_rate.fits'],
    'jw02072-o002_20221206t143745_spec2_00001_asn.json':
        ['jw02072002001_05101_00001_nrs1_rate.fits',
         'jw02072002001_05101_00002_nrs1_rate.fits',
         'jw02072002001_05101_00003_nrs1_rate.fits']
}
# ids = ["fullframe", "S400A1-subarray", "ALLSLITS-subarray"]


@pytest.fixture(scope="module", params=file_roots)# ids=ids)
def run_pipeline(jail, rtdata_module, request):
    """Run the calwebb_spec2 pipeline on NIRSpec Fixed-Slit exposures.
       We currently test the following types of inputs:
         1) Full-frame exposure (all slits will be extracted)
         2) ALLSLITS subarray exposure (all slits will be extracted)
         3) S400A1 subarray exposure (1 slit extracted)"""

    rtdata = rtdata_module

    for fle in asn_memberdict[request.param]:
        rtdata.get_data('nirspec/fs/' + fle)
    # Get the input exposure
    rtdata.get_data('nirspec/fs/' + request.param)

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = ["calwebb_spec2", rtdata.input,
            "--steps.assign_wcs.save_results=true",
            "--steps.extract_2d.save_results=true",
            "--steps.wavecorr.save_results=true",
            "--steps.srctype.save_results=true",
            "--steps.flat_field.save_results=true",
            "--steps.pathloss.save_results=true"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", [
    "assign_wcs", "extract_2d", "wavecorr", "flat_field", "pathloss", "srctype",
    "cal", "s2d", "x1d"])
def test_nirspec_fs_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on a
       NIRSpec FS exposures."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = replace_suffix(
        os.path.splitext(asn_memberdict[os.path.basename(rtdata.input)][0])[0], suffix) + '.fits'
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_fs_spec2", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_pathloss_corrpars(rtdata):
    """Test PathLossStep using correction_pars"""
    with dm.open(rtdata.get_data('nirspec/fs/nrs1_flat_field.fits')) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.use_correction_pars = True
        corrected_corrpars = pls.run(data)

    bad_slits = []
    for idx, slits in enumerate(zip(corrected.slits, corrected_corrpars.slits)):
        corrected_slit, corrected_corrpars_slit = slits
        if not np.allclose(corrected_slit.data, corrected_corrpars_slit.data, equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f'correction_pars failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_pathloss_inverse(rtdata):
    """Test PathLossStep using inversion"""
    with dm.open(rtdata.get_data('nirspec/fs/nrs1_flat_field.fits')) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.inverse = True
        corrected_inverse = pls.run(corrected)

        bad_slits = []
        for idx, slits in enumerate(zip(data.slits, corrected_inverse.slits)):
            data_slit, corrected_inverse_slit = slits
            non_nan = ~np.isnan(corrected_inverse_slit.data)
            if not np.allclose(data_slit.data[non_nan], corrected_inverse_slit.data[non_nan]):
                bad_slits.append(idx)

    assert not bad_slits, f'Inversion failed for slits {bad_slits}'


@pytest.mark.bigdata
def test_pathloss_source_type(rtdata):
    """Test PathLossStep forcing source type"""
    with dm.open(rtdata.get_data('nirspec/fs/nrs1_flat_field.fits')) as data:
        pls = PathLossStep()
        pls.source_type = 'extended'
        pls.run(data)

    bad_slits = []
    for idx, slit in enumerate(pls.correction_pars.slits):
        if slit:
            if not np.allclose(slit.data, slit.pathloss_uniform, equal_nan=True):
                bad_slits.append(idx)
    assert not bad_slits, f'Force to uniform failed for slits {bad_slits}'
