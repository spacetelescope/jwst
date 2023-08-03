"""
Unit test for Cube Build for reading nirspec cube pars ref file and using it
"""

import numpy as np
import pytest
import math
from astropy.io import fits
from jwst.cube_build import ifu_cube
from jwst.cube_build import cube_build_io_util
from jwst.cube_build import instrument_defaults


@pytest.fixture(scope='module')
def nirspec_cube_pars(tmpdir_factory):
    """ Set up the nirspec cube pars reference file  """

    filename = tmpdir_factory.mktemp('cube_pars')
    filename = str(filename.join('nirspec_cube_pars.fits'))
    hdu0 = fits.PrimaryHDU()
    hdu0.header['REFTYPE'] = 'CUBEPAR'
    hdu0.header['INSTRUME'] = 'NIRSPEC'
    hdu0.header['MODELNAM'] = 'FM'
    hdu0.header['DETECTOR'] = 'N/A'
    hdu0.header['EXP_TYPE'] = 'NRS_IFU'
    hdu0.header['BAND'] = 'N/A'
    hdu0.header['CHANNEL'] = 'N/A'

    # make the first extension
    disp = np.array(['PRISM', 'G140M', 'G140M', 'G140H', 'G140H', 'G235M', 'G235H', 'G395M', 'G395H'])
    filt = np.array(['CLEAR', 'F070LP', 'F100LP', 'F070LP', 'F100LP', 'F170LP', 'F170LP', 'F290LP', 'F290LP'])
    spsize = np.array([0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10])
    wsamp = np.array([0.005, 0.001, 0.0006, 0.0002, 0.0002, 0.001, 0.0004, 0.0017, 0.0007])

    wmin = np.array([0.6, 0.7, 0.97, 0.7, 0.97, 1.66, 1.66, 2.87, 2.87])
    wmax = np.array([5.3, 1.27, 1.89, 1.27, 1.89, 3.17, 3.17, 5.27, 5.27])

    col1 = fits.Column(name='DISPERSER', format='5A', array=disp)
    col2 = fits.Column(name='FILTER', format='6A', array=filt)
    col3 = fits.Column(name='WAVEMIN', format='E', array=wmin, unit='micron')
    col4 = fits.Column(name='WAVEMAX', format='E', array=wmax, unit='micron')
    col5 = fits.Column(name='SPAXELSIZE', format='E', array=spsize, unit='arcsec')
    col6 = fits.Column(name='SPECTRALSTEP', format='D', array=wsamp, unit='micron')

    hdu1 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    hdu1.header['EXTNAME'] = 'CUBEPAR'

    # make the second extension
    disp = np.array(['PRISM', 'G140M', 'G140M', 'G140H', 'G140H', 'G235M', 'G235H', 'G395M', 'G395H'])
    filt = np.array(['CLEAR', 'F070LP', 'F100LP', 'F070LP', 'F100LP', 'F170LP', 'F170LP', 'F290LP', 'F290LP'])

    roispat = np.array([0.201, 0.202, 0.203, 0.204, 0.205, 0.206, 0.207, 0.208, 0.209])
    roispec = np.array([0.011, 0.0012, 0.0012, 0.0004, 0.0004, 0.002, 0.0008, 0.003, 0.0012])

    power = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2])
    softrad = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])

    col1 = fits.Column(name='DISPERSER', format='5A', array=disp)
    col2 = fits.Column(name='FILTER', format='6A', array=filt)
    col3 = fits.Column(name='ROISPATIAL', format='E', array=roispat, unit='arcsec')
    col4 = fits.Column(name='ROISPECTRAL', format='E', array=roispec, unit='micron')
    col5 = fits.Column(name='POWER', format='I', array=power)
    col6 = fits.Column(name='SOFTRAD', format='E', array=softrad, unit='arcsec')

    hdu2 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    hdu2.header['EXTNAME'] = 'CUBEPAR_MSM'

    # make the third extension
    # Define the multiextension wavelength solution
    finalwave = np.arange(0.6, 5.3, 0.1)
    nelem = len(finalwave)

    # Linear relation of spatial roi with wavelength
    roispat = np.ones(nelem) * 0.2
    # Linear relation of spectral roi with wavelength
    roispec = np.ones(nelem) * 0.01
    # Power is 2 at all wavelengths
    power = np.ones(nelem, dtype=int) * 2
    # Softening radius is 0.01 at all wavelengths
    softrad = np.ones(nelem) * 0.01

    col1 = fits.Column(name='WAVELENGTH', format='D', array=finalwave, unit='micron')
    col2 = fits.Column(name='ROISPATIAL', format='E', array=roispat, unit='arcsec')
    col3 = fits.Column(name='ROISPECTRAL', format='E', array=roispec, unit='micron')
    col4 = fits.Column(name='POWER', format='I', array=power)
    col5 = fits.Column(name='SOFTRAD', format='E', array=softrad, unit='arcsec')

    hdu3 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])
    hdu3.header['EXTNAME'] = 'MULTICHAN_PRISM_MSM'

    # make the 4th extension
    # Define the multiextension wavelength solution
    finalwave = np.arange(0.7, 7.7, 0.1)
    nelem = len(finalwave)
    # Linear relation of spatial roi with wavelength
    roispat = np.ones(nelem) * 0.2
    # Linear relation of spectral roi with wavelength
    roispec = np.ones(nelem) * 0.01
    # Power is 2 at all wavelengths
    power = np.ones(nelem, dtype=int) * 2
    # Softening radius is 0.01 at all wavelengths
    softrad = np.ones(nelem) * 0.01

    col1 = fits.Column(name='WAVELENGTH', format='D', array=finalwave, unit='micron')
    col2 = fits.Column(name='ROISPATIAL', format='E', array=roispat, unit='arcsec')
    col3 = fits.Column(name='ROISPECTRAL', format='E', array=roispec, unit='micron')
    col4 = fits.Column(name='POWER', format='I', array=power)
    col5 = fits.Column(name='SOFTRAD', format='E', array=softrad, unit='arcsec')

    hdu4 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])
    hdu4.header['EXTNAME'] = 'MULTICHAN_MED_MSM'

    hdu = fits.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4])
    hdu.writeto(filename, overwrite=True)
    return filename


def test_nirspec_cubepars(_jail, nirspec_cube_pars):
    """ Read in the nirspec cube pars file """

    instrument_info = instrument_defaults.InstrumentInfo()
    all_channel = []
    all_subchannel = []
    all_grating = []
    all_filter = []
    all_grating.append('prism')
    all_filter.append('clear')

    cube_build_io_util.read_cubepars(nirspec_cube_pars,
                                     'NIRSPEC',
                                     'msm',
                                     all_channel,
                                     all_subchannel,
                                     all_grating,
                                     all_filter,
                                     instrument_info)

    par1 = 'prism'
    par2 = 'clear'
    ascale, bscale, wscale = instrument_info.GetScale(par1, par2)
    assert math.isclose(ascale, 0.1, abs_tol=0.00001)
    assert math.isclose(bscale, 0.1, abs_tol=0.00001)
    assert math.isclose(wscale, 0.005, abs_tol=0.00001)

    roiw = instrument_info.GetWaveRoi(par1, par2)
    rois = instrument_info.GetSpatialRoi(par1, par2)
    power = instrument_info.GetMSMPower(par1, par2)
    wavemin = instrument_info.GetWaveMin(par1, par2)
    wavemax = instrument_info.GetWaveMax(par1, par2)

    assert math.isclose(roiw, 0.011, abs_tol=0.00001)
    assert math.isclose(rois, 0.201, abs_tol=0.00001)
    assert math.isclose(power, 2, abs_tol=0.00001)
    assert math.isclose(wavemin, 0.6, abs_tol=0.00001)
    assert math.isclose(wavemax, 5.3, abs_tol=0.00001)

    # set up the ifucube class

    pars_cube = {
        'scalexy': 0.0,
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
        'debug_spaxel': '0 0 0'}

    pipeline = 3
    input_model = None
    output_name_base = None
    output_type = 'band'
    instrument = 'NIRSPEC'
    list_par1 = all_grating
    list_par2 = all_filter
    master_table = None
    instrument_info = instrument_info
    this_cube = ifu_cube.IFUCubeData(
        pipeline,
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
    # linear_wavelength = 'True'
    this_cube.determine_cube_parameters()

    assert math.isclose(this_cube.wavemin, wavemin, abs_tol=0.00001)
    assert math.isclose(this_cube.wavemax, wavemax, abs_tol=0.00001)
    assert this_cube.linear_wavelength is True

    assert math.isclose(this_cube.weight_power, 2, abs_tol=0.00001)
    assert math.isclose(this_cube.roiw, 0.011, abs_tol=0.00001)
    rois = 0.201 * 1.5  # increase for single file
    assert math.isclose(this_cube.rois, rois, abs_tol=0.00001)
    assert math.isclose(this_cube.spatial_size, 0.1, abs_tol=0.00001)

    # now test if the user has provided input to build cube

    user_ascale = 0.2
    user_wscale = 0.05
    user_power = 1
    user_wave_min = 0.8
    user_wave_max = 4.5
    user_rois = 0.6
    user_roiw = 0.8
    pars_cube = {
        'scalexy': user_ascale,
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
        'debug_spaxel': '0 0 0'}

    this_cube = ifu_cube.IFUCubeData(
        pipeline,
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
    assert this_cube.linear_wavelength is True
    assert math.isclose(this_cube.spatial_size, user_ascale, abs_tol=0.00001)
    assert math.isclose(this_cube.spectral_size, user_wscale, abs_tol=0.00001)
    assert math.isclose(this_cube.weight_power, user_power, abs_tol=0.00001)
    assert math.isclose(this_cube.roiw, user_roiw, abs_tol=0.00001)
    assert math.isclose(this_cube.rois, user_rois, abs_tol=0.00001)
