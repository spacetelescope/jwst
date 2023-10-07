"""
Unit test for Cube Build testing reading in MIRI cubepars ref file and using it
"""

import numpy as np
import pytest
from astropy.io import fits

from gwcs import WCS
from stdatamodels.jwst.datamodels import IFUImageModel

from jwst import assign_wcs
from jwst.cube_build import CubeBuildStep
from jwst.cube_build.file_table import ErrorNoAssignWCS
from jwst.cube_build.cube_build import ErrorNoChannels


@pytest.fixture(scope='module')
def miri_cube_pars(tmpdir_factory):
    """ Set up the miri cube pars reference file  """

    filename = tmpdir_factory.mktemp('cube_pars')
    filename = str(filename.join('miri_cube_pars.fits'))
    hdu0 = fits.PrimaryHDU()
    hdu0.header['REFTYPE'] = 'CUBEPAR'
    hdu0.header['INSTRUME'] = 'MIRI'
    hdu0.header['MODELNAM'] = 'FM'
    hdu0.header['DETECTOR'] = 'N/A'
    hdu0.header['EXP_TYPE'] = 'MIR_MRS'

    # make the first extension
    channel = np.array(['1', '1', '1', '2', '2', '2', '3', '3', '3', '4', '4', '4'])
    subchannel = np.array(['SHORT', 'MEDIUM', 'LONG', 'SHORT', 'MEDIUM', 'LONG',
                           'SHORT', 'MEDIUM', 'LONG', 'SHORT', 'MEDIUM', 'LONG'])

    spsize = np.array([0.13, 0.13, 0.13, 0.17, 0.17, 0.17, 0.2, 0.2, 0.2, 0.35, 0.35, 0.35])
    wsamp = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.006, 0.006, 0.006])

    wmin = np.array([4.89, 5.65, 6.52, 7.49, 8.65, 10.00, 11.53, 13.37, 15.44, 17.66, 20.54, 23.95])
    wmax = np.array([5.75, 6.64, 7.66, 8.78, 10.14, 11.7, 13.48, 15.63, 18.05, 20.92, 24.40, 28.45])

    col1 = fits.Column(name='CHANNEL', format='1A', array=channel)
    col2 = fits.Column(name='BAND', format='6A', array=subchannel)
    col3 = fits.Column(name='WAVEMIN', format='E', array=wmin, unit='micron')
    col4 = fits.Column(name='WAVEMAX', format='E', array=wmax, unit='micron')
    col5 = fits.Column(name='SPAXELSIZE', format='E', array=spsize, unit='arcsec')
    col6 = fits.Column(name='SPECTRALSTEP', format='D', array=wsamp, unit='micron')

    hdu1 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    hdu1.header['EXTNAME'] = 'CUBEPAR'

    # make the second extension
    roispat = np.array([0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.20, 0.20, 0.20, 0.40, 0.40, 0.40])
    roispec = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.006, 0.006, 0.006])

    power = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
    softrad = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])

    col1 = fits.Column(name='CHANNEL', format='1A', array=channel)
    col2 = fits.Column(name='BAND', format='6A', array=subchannel)
    col3 = fits.Column(name='ROISPATIAL', format='E', array=roispat, unit='arcsec')
    col4 = fits.Column(name='ROISPECTRAL', format='E', array=roispec, unit='micron')
    col5 = fits.Column(name='POWER', format='I', array=power)
    col6 = fits.Column(name='SOFTRAD', format='E', array=softrad, unit='arcsec')

    hdu2 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    hdu2.header['EXTNAME'] = 'CUBEPAR_MSM'

    # make the third extension
    # Define the multiextension wavelength solution - only use a few number for testing
    finalwave = np.array([5, 10, 15, 20, 25])
    roispat = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    roispec = np.array([0.001, 0.002, 0.003, 0.004, 0.005])
    power = np.array([1, 2, 3, 4, 5])
    softrad = np.array([0.01, 0.02, 0.03, 0.04, 0.05])

    col1 = fits.Column(name='WAVELENGTH', format='D', array=finalwave, unit='micron')
    col2 = fits.Column(name='ROISPATIAL', format='E', array=roispat, unit='arcsec')
    col3 = fits.Column(name='ROISPECTRAL', format='E', array=roispec, unit='micron')
    col4 = fits.Column(name='POWER', format='I', array=power)
    col5 = fits.Column(name='SOFTRAD', format='E', array=softrad, unit='arcsec')

    hdu3 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])
    hdu3.header['EXTNAME'] = 'MULTICHANNEL_MSM'

    hdu = fits.HDUList([hdu0, hdu1, hdu2, hdu3])
    hdu.writeto(filename, overwrite=True)
    return filename


@pytest.fixture(scope='function')
def miri_image():

    image = IFUImageModel((20, 20))
    image.data = np.random.random((20, 20))
    image.meta.instrument.name = 'MIRI'
    image.meta.instrument.detector = 'MIRIFULONG'
    image.meta.exposure.type = 'MIR_MRS'
    image.meta.instrument.channel = '12'
    image.meta.instrument.band = 'SHORT'
    image.meta.filename = 'test_miri.fits'
    return image


@pytest.mark.parametrize("as_filename", [True, False])
def test_call_cube_build(_jail, miri_cube_pars, miri_image, tmp_path, as_filename):
    """ test defaults of step are set up and user input are defined correctly """
    if as_filename:
        fn = tmp_path / 'miri.fits'
        miri_image.save(fn)
        step_input = fn
    else:
        step_input = miri_image

    # we do not want to run the CubeBuild through to completion because
    # the image needs to be a full image and this take too much time
    # in a unit test

    # Test ErrorNoAssignWCS is raised
    with pytest.raises(ErrorNoAssignWCS):
        step = CubeBuildStep()
        step.override_cubepar = miri_cube_pars
        step.channel = '3'
        step.run(step_input)

    # Test some defaults to step are setup correctly and
    # is user specifies channel is set up correctly
    step = CubeBuildStep()
    step.override_cubepar = miri_cube_pars
    step.channel = '1'

    try:
        step.run(step_input)
    except ErrorNoAssignWCS:
        pass

    assert step.pars_input['channel'] == ['1']
    assert step.interpolation == 'drizzle'
    assert step.weighting == 'drizzle'
    assert step.coord_system == 'skyalign'

    # Set Assign WCS has been run but the user input to channels is wrong
    miri_image.meta.cal_step.assign_wcs = 'COMPLETE'
    # save file with modifications
    if as_filename:
        miri_image.save(step_input)
    with pytest.raises(ErrorNoChannels):
        step = CubeBuildStep()
        step.override_cubepar = miri_cube_pars
        step.channel = '3'
        step.run(step_input)


@pytest.fixture(scope='function')
def nirspec_data():
    image = IFUImageModel((2048, 2048))
    image.data = np.random.random((2048, 2048))
    image.meta.instrument.name = 'NIRSPEC'
    image.meta.instrument.detector = 'NRS1'
    image.meta.exposure.type = 'NRS_IFU'
    image.meta.filename = 'test_nirspec.fits'
    image.meta.observation.date = '2023-10-06'
    image.meta.observation.time = '00:00:00.000'
    # below values taken from regtest using file
    # jw01249005001_03101_00004_nrs1_cal.fits
    image.meta.instrument.filter = 'F290LP'
    image.meta.instrument.grating = 'G395H'
    image.meta.wcsinfo.v2_ref = 299.83548
    image.meta.wcsinfo.v3_ref = -498.256805
    image.meta.wcsinfo.ra_ref = 358.0647567841019
    image.meta.wcsinfo.dec_ref = -2.167207258876695
    image.meta.cal_step.assign_wcs = 'COMPLETE'
    step = assign_wcs.assign_wcs_step.AssignWcsStep()
    refs = {}
    for reftype in assign_wcs.assign_wcs_step.AssignWcsStep.reference_file_types:
        refs[reftype] = step.get_reference_file(image, reftype)
    pipe = assign_wcs.nirspec.create_pipeline(image, refs, slit_y_range=[-.5, .5])
    image.meta.wcs = WCS(pipe)
    return image


@pytest.mark.parametrize("as_filename", [True, False])
def test_call_cube_build_nirspec(_jail, nirspec_data, tmp_path, as_filename):
    if as_filename:
        fn = tmp_path / 'test_nirspec.fits'
        nirspec_data.save(fn)
        step_input = fn
    else:
        step_input = nirspec_data
    step = CubeBuildStep()
    step.channel = '1'
    step.coord_system = 'internal_cal'
    step.run(step_input)
