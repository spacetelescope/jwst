"""
Unit test for Cube Build testing reading in MIRI cubepars ref file and using it
"""

import numpy as np
import pytest
import math
from astropy.io import fits
from jwst.cube_build import ifu_cube
from jwst.cube_build import cube_build_io_util
from jwst.cube_build import instrument_defaults


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
    hdu.close()

    return filename


def test_miri_use_cubepars(_jail, miri_cube_pars):
    """ Test reading in the miri cube pars file """

    instrument_info = instrument_defaults.InstrumentInfo()
    all_channel = []
    all_subchannel = []
    all_grating = []
    all_filter = []
    all_channel.append('1')
    all_subchannel.append('medium')

    cube_build_io_util.read_cubepars(miri_cube_pars,
                                     'MIRI',
                                     'msm',
                                     all_channel,
                                     all_subchannel,
                                     all_grating,
                                     all_filter,
                                     instrument_info)

    par1 = '1'
    par2 = 'medium'

    ascale, bscale, wscale = instrument_info.GetScale(par1, par2)

    # check that the values are read in correctly
    assert math.isclose(ascale, 0.13, abs_tol=0.00001)
    assert math.isclose(bscale, 0.13, abs_tol=0.00001)
    assert math.isclose(wscale, 0.001, abs_tol=0.00001)

    roiw = instrument_info.GetWaveRoi(par1, par2)
    rois = instrument_info.GetSpatialRoi(par1, par2)
    power = instrument_info.GetMSMPower(par1, par2)
    wavemin = instrument_info.GetWaveMin(par1, par2)
    wavemax = instrument_info.GetWaveMax(par1, par2)

    assert math.isclose(roiw, 0.001, abs_tol=0.00001)
    assert math.isclose(rois, 0.1, abs_tol=0.00001)
    assert math.isclose(power, 2, abs_tol=0.00001)
    assert math.isclose(wavemin, 5.65, abs_tol=0.00001)
    assert math.isclose(wavemax, 6.64, abs_tol=0.00001)

    # set up the ifucube class
    pars_cube = {
        'scale1': 0.0,
        'scale2': 0.0,
        'scalew': 0.0,
        'interpolation': 'pointcloud',
        'weighting': 'msm',
        'weight_power': 2,
        'coord_system': 'world',
        'rois': 0.0,
        'roiw': 0.0,
        'wavemin': None,
        'wavemax': None,
        'skip_dqflagging': False,
        'xdebug': None,
        'ydebug': None,
        'zdebug': None,
        'debug_pixel': 0,
        'spaxel_debug': None}

    pipeline = 3
    filename = None
    input_model = None
    output_name_base = None
    output_type = 'band'
    instrument = 'MIRI'
    list_par1 = all_channel
    list_par2 = all_subchannel
    master_table = None
    instrument_info = instrument_info
    this_cube = ifu_cube.IFUCubeData(
        pipeline,
        filename,
        input_model,
        output_name_base,
        output_type,
        instrument,
        list_par1,
        list_par2,
        instrument_info,
        master_table,
        **pars_cube)

    this_cube.num_files = 1  # set in ifu cube
    # test that the correct values read from the table are filled
    # in in the this_cube class. Also check that for this configuration
    # linear_wavelength = True
    this_cube.determine_cube_parameters()

    assert math.isclose(this_cube.wavemin, wavemin, abs_tol=0.00001)
    assert math.isclose(this_cube.wavemax, wavemax, abs_tol=0.00001)
    assert this_cube.linear_wavelength

    assert math.isclose(this_cube.weight_power, 2, abs_tol=0.00001)
    assert math.isclose(this_cube.roiw, 0.001, abs_tol=0.00001)
    rois = 0.1 * 1.5  # increase for single file
    assert math.isclose(this_cube.rois, rois, abs_tol=0.00001)
    assert math.isclose(this_cube.spatial_size, 0.13, abs_tol=0.00001)


def test_miri_cubepars_user_defaults(_jail, miri_cube_pars):
    """ Read in the miri cube pars file and override some defaults """

    instrument_info = instrument_defaults.InstrumentInfo()
    all_channel = []
    all_subchannel = []
    all_grating = []
    all_filter = []
    all_channel.append('4')
    all_subchannel.append('long')

    cube_build_io_util.read_cubepars(miri_cube_pars,
                                     'MIRI',
                                     'msm',
                                     all_channel,
                                     all_subchannel,
                                     all_grating,
                                     all_filter,
                                     instrument_info)

    # test another band
    par1 = '4'
    par2 = 'long'

    # first check that it reads in correct values for this band
    # from the reference file
    ascale, bscale, wscale = instrument_info.GetScale(par1, par2)

    assert math.isclose(ascale, 0.35, abs_tol=0.00001)
    assert math.isclose(bscale, 0.35, abs_tol=0.00001)
    assert math.isclose(wscale, 0.006, abs_tol=0.00001)

    roiw = instrument_info.GetWaveRoi(par1, par2)
    rois = instrument_info.GetSpatialRoi(par1, par2)
    power = instrument_info.GetMSMPower(par1, par2)
    wavemin = instrument_info.GetWaveMin(par1, par2)
    wavemax = instrument_info.GetWaveMax(par1, par2)

    assert math.isclose(roiw, 0.006, abs_tol=0.00001)
    assert math.isclose(rois, 0.4, abs_tol=0.00001)
    assert math.isclose(power, 2, abs_tol=0.00001)
    assert math.isclose(wavemin, 23.95, abs_tol=0.00001)
    assert math.isclose(wavemax, 28.45, abs_tol=0.00001)

    # set up the ifucube class
    pars_cube = {
        'scale1': 0.0,
        'scale2': 0.0,
        'scalew': 0.0,
        'interpolation': 'pointcloud',
        'weighting': 'msm',
        'weight_power': 2,
        'coord_system': 'world',
        'rois': 0.0,
        'roiw': 0.0,
        'wavemin': None,
        'wavemax': None,
        'skip_dqflagging': False,
        'xdebug': None,
        'ydebug': None,
        'zdebug': None,
        'debug_pixel': 0,
        'spaxel_debug': None}

    pipeline = 3
    filename = None
    input_model = None
    output_name_base = None
    output_type = 'band'
    instrument = 'MIRI'
    list_par1 = all_channel
    list_par2 = all_subchannel
    master_table = None
    instrument_info = instrument_info
    this_cube = ifu_cube.IFUCubeData(
        pipeline,
        filename,
        input_model,
        output_name_base,
        output_type,
        instrument,
        list_par1,
        list_par2,
        instrument_info,
        master_table,
        **pars_cube)

    this_cube.num_files = 1  # set in ifu cube
    # test that the correct values read from the table are filled
    # in in the this_cube class. Also check that for this configuration
    # linear_wavelength = True
    this_cube.determine_cube_parameters()
    # now test if the user has provided input to build cube

    assert math.isclose(this_cube.wavemin, wavemin, abs_tol=0.00001)
    assert math.isclose(this_cube.wavemax, wavemax, abs_tol=0.00001)
    assert this_cube.linear_wavelength

    assert math.isclose(this_cube.weight_power, 2, abs_tol=0.00001)
    assert math.isclose(this_cube.roiw, roiw, abs_tol=0.00001)
    rois = rois * 1.5  # increase for single file
    assert math.isclose(this_cube.rois, rois, abs_tol=0.00001)
    assert math.isclose(this_cube.spatial_size, 0.35, abs_tol=0.00001)

    user_ascale = 0.2
    user_wscale = 0.05
    user_power = 1
    user_wave_min = 24.5
    user_wave_max = 27.5
    user_rois = 0.6
    user_roiw = 0.8
    pars_cube = {
        'scale1': user_ascale,
        'scale2': user_ascale,
        'scalew': user_wscale,
        'interpolation': 'pointcloud',
        'weighting': 'msm',
        'weight_power': user_power,
        'coord_system': 'world',
        'rois': user_rois,
        'roiw': user_roiw,
        'wavemin': user_wave_min,
        'wavemax': user_wave_max,
        'skip_dqflagging': False,
        'xdebug': None,
        'ydebug': None,
        'zdebug': None,
        'debug_pixel': 0,
        'spaxel_debug': None}

    this_cube = ifu_cube.IFUCubeData(
        pipeline,
        filename,
        input_model,
        output_name_base,
        output_type,
        instrument,
        list_par1,
        list_par2,
        instrument_info,
        master_table,
        **pars_cube)

    this_cube.num_files = 1  # set in check_ifucube
    this_cube.determine_cube_parameters()
    # do they match the user provided ones
    assert math.isclose(this_cube.wavemin, user_wave_min, abs_tol=0.00001)
    assert math.isclose(this_cube.wavemax, user_wave_max, abs_tol=0.00001)
    assert this_cube.linear_wavelength
    assert math.isclose(this_cube.spatial_size, user_ascale, abs_tol=0.00001)
    assert math.isclose(this_cube.spectral_size, user_wscale, abs_tol=0.00001)
    assert math.isclose(this_cube.weight_power, user_power, abs_tol=0.00001)
    assert math.isclose(this_cube.roiw, user_roiw, abs_tol=0.00001)
    assert math.isclose(this_cube.rois, user_rois, abs_tol=0.00001)


def test_miri_cubepars_multiple_bands(_jail, miri_cube_pars):
    """Read in the miri cube pars file. Test cube has correct values when
    multiple bands are used
    """

    instrument_info = instrument_defaults.InstrumentInfo()
    all_channel = []
    all_subchannel = []
    all_grating = []
    all_filter = []
    # set up all_channel and all_subchannel - 1 to 1 matching between the two
    all_channel.append('1')
    all_channel.append('1')
    all_channel.append('1')
    all_channel.append('2')
    all_channel.append('2')
    all_channel.append('2')
    all_channel.append('3')
    all_channel.append('3')
    all_channel.append('3')

    all_subchannel.append('short')
    all_subchannel.append('medium')
    all_subchannel.append('long')
    all_subchannel.append('short')
    all_subchannel.append('medium')
    all_subchannel.append('long')
    all_subchannel.append('short')
    all_subchannel.append('medium')
    all_subchannel.append('long')

    cube_build_io_util.read_cubepars(miri_cube_pars,
                                     'MIRI',
                                     'msm',
                                     all_channel,
                                     all_subchannel,
                                     all_grating,
                                     all_filter,
                                     instrument_info)

    # test reading in another band we have not checked before
    par1 = '3'
    par2 = 'medium'

    # first check that it reads in correct values for this band
    # from the reference file
    ascale, bscale, wscale = instrument_info.GetScale(par1, par2)

    assert math.isclose(ascale, 0.2, abs_tol=0.00001)
    assert math.isclose(bscale, 0.2, abs_tol=0.00001)
    assert math.isclose(wscale, 0.003, abs_tol=0.00001)

    roiw = instrument_info.GetWaveRoi(par1, par2)
    rois = instrument_info.GetSpatialRoi(par1, par2)
    power = instrument_info.GetMSMPower(par1, par2)
    wavemin = instrument_info.GetWaveMin(par1, par2)
    wavemax = instrument_info.GetWaveMax(par1, par2)

    assert math.isclose(roiw, 0.003, abs_tol=0.00001)
    assert math.isclose(rois, 0.2, abs_tol=0.00001)
    assert math.isclose(power, 2, abs_tol=0.00001)
    assert math.isclose(wavemin, 13.37, abs_tol=0.00001)
    assert math.isclose(wavemax, 15.63, abs_tol=0.00001)

    # set up the ifucube class
    pars_cube = {
        'scale1': 0.0,
        'scale2': 0.0,
        'scalew': 0.0,
        'interpolation': 'pointcloud',
        'weighting': 'msm',
        'weight_power': 2,
        'coord_system': 'world',
        'rois': 0.0,
        'roiw': 0.0,
        'wavemin': None,
        'wavemax': None,
        'skip_dqflagging': False,
        'xdebug': None,
        'ydebug': None,
        'zdebug': None,
        'debug_pixel': 0,
        'spaxel_debug': None}

    pipeline = 3
    filename = None
    input_model = None
    output_name_base = None
    output_type = 'multi'
    instrument = 'MIRI'
    list_par1 = all_channel
    list_par2 = all_subchannel
    master_table = None
    instrument_info = instrument_info
    this_cube = ifu_cube.IFUCubeData(
        pipeline,
        filename,
        input_model,
        output_name_base,
        output_type,
        instrument,
        list_par1,
        list_par2,
        instrument_info,
        master_table,
        **pars_cube)

    this_cube.num_files = 12  # set in ifu cube
    # test that the correct values read from the table are filled
    # in in the this_cube class when multiple bands and output_type = multi
    # are set. Also check that for this configuration linear_wavelength = False
    this_cube.determine_cube_parameters()

    # for multiple bands the smallest spatial scale is chosen
    assert math.isclose(this_cube.spatial_size, 0.13, abs_tol=0.00001)
    assert not this_cube.linear_wavelength

    # wavemin - min for channels 1-3 (min of channel 1 short)
    assert math.isclose(this_cube.wavemin, 4.89, abs_tol=0.00001)
    # wavemas = max for channels 1-3 (max of channel 3 long)
    assert math.isclose(this_cube.wavemax, 18.05, abs_tol=0.00001)

    weight_test = np.array([1, 2, 3, 4])
    roiw_test = np.array([0.001, 0.002, 0.003, 0.004])
    rois_test = np.array([0.1, 0.2, 0.3, 0.4])
    wave_test = np.array([5, 10, 15, 20])
    assert np.allclose(this_cube.weight_power_table, weight_test, rtol=0.00001)
    assert np.allclose(this_cube.roiw_table, roiw_test, rtol=0.00001)
    assert np.allclose(this_cube.rois_table, rois_test, rtol=0.00001)
    assert np.allclose(this_cube.wavelength_table, wave_test, rtol=0.00001)
