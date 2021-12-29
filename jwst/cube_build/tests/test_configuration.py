"""
Unit test for Cube Build testing setting up configuration
"""

import pytest

from jwst import datamodels
from jwst.cube_build import cube_build
from jwst.cube_build import file_table

wcsinfo = {
    'dec_ref': -0.00244536159612126,
    'ra_ref': -0.00205553321270217,
    'roll_ref': -0.0,
    'v2_ref': -503.65447,
    'v3_ref': -318.74246,
    'v3yangle': 0.0,
    'vparity': -1
}

mirifushort_short = {
    'detector': 'MIRIFUSHORT',
    'channel': '12',
    'band': 'SHORT',
    'name': 'MIRI'
}

mirifushort_medium = {
    'detector': 'MIRIFUSHORT',
    'channel': '12',
    'band': 'MEDIUM',
    'name': 'MIRI'
}

mirifushort_long = {
    'detector': 'MIRIFUSHORT',
    'channel': '12',
    'band': 'LONG',
    'name': 'MIRI'
}

mirifulong_short = {
    'detector': 'MIRIFULONG',
    'channel': '34',
    'band': 'SHORT',
    'name': 'MIRI'
}

mirifulong_medium = {
    'detector': 'MIRIFULONG',
    'channel': '34',
    'band': 'MEDIUM',
    'name': 'MIRI'
}

mirifulong_long = {
    'detector': 'MIRIFULONG',
    'channel': '34',
    'band': 'LONG',
    'name': 'MIRI'
}

nirspec_G140M = {
    'detector': 'NRS1',
    'grating': 'G140M',
    'filter': 'F100LP',
    'name': 'NIRSPEC'
}

nirspec_G235M = {
    'detector': 'NRS1',
    'grating': 'G235M',
    'filter': 'F170LP',
    'name': 'NIRSPEC'
}

observation = {
    'date': '2019-01-01',
    'time': '17:00:00'}

subarray = {
    'fastaxis': 1,
    'name': 'FULL',
    'slowaxis': 2,
    'xsize': 1032,
    'xstart': 1,
    'ysize': 1024,
    'ystart': 1
}

subarray_nirspec = {
    'fastaxis': 2,
    'name': 'FULL',
    'slowaxis': 1,
    'xsize': 2048,
    'xstart': 1,
    'ysize': 2048,
    'ystart': 1
}


@pytest.fixture(scope='function')
def miri_ifushort_short():
    """ Generate a IFU image """

    input_model = datamodels.IFUImageModel()
    input_model.meta.wcsinfo._instance.update(wcsinfo)
    input_model.meta.instrument._instance.update(mirifushort_short)
    input_model.meta.observation._instance.update(observation)
    input_model.meta.subarray._instance.update(subarray)
    input_model.meta.cal_step.assign_wcs = 'COMPLETE'
    return input_model


@pytest.fixture(scope='function')
def miri_full_coverage():
    """ Generate a IFU images SHORT, LONG for all three bands """

    input_model1 = datamodels.IFUImageModel()
    input_model1.meta.wcsinfo._instance.update(wcsinfo)
    input_model1.meta.instrument._instance.update(mirifushort_short)
    input_model1.meta.observation._instance.update(observation)
    input_model1.meta.subarray._instance.update(subarray)
    input_model1.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model2 = datamodels.IFUImageModel()
    input_model2.meta.wcsinfo._instance.update(wcsinfo)
    input_model2.meta.instrument._instance.update(mirifushort_medium)
    input_model2.meta.observation._instance.update(observation)
    input_model2.meta.subarray._instance.update(subarray)
    input_model2.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model3 = datamodels.IFUImageModel()
    input_model3.meta.wcsinfo._instance.update(wcsinfo)
    input_model3.meta.instrument._instance.update(mirifushort_long)
    input_model3.meta.observation._instance.update(observation)
    input_model3.meta.subarray._instance.update(subarray)
    input_model3.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model4 = datamodels.IFUImageModel()
    input_model4.meta.wcsinfo._instance.update(wcsinfo)
    input_model4.meta.instrument._instance.update(mirifulong_short)
    input_model4.meta.observation._instance.update(observation)
    input_model4.meta.subarray._instance.update(subarray)
    input_model4.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model5 = datamodels.IFUImageModel()
    input_model5.meta.wcsinfo._instance.update(wcsinfo)
    input_model5.meta.instrument._instance.update(mirifulong_medium)
    input_model5.meta.observation._instance.update(observation)
    input_model5.meta.subarray._instance.update(subarray)
    input_model5.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model6 = datamodels.IFUImageModel()
    input_model6.meta.wcsinfo._instance.update(wcsinfo)
    input_model6.meta.instrument._instance.update(mirifulong_long)
    input_model6.meta.observation._instance.update(observation)
    input_model6.meta.subarray._instance.update(subarray)
    input_model6.meta.cal_step.assign_wcs = 'COMPLETE'

    # full range and 2 dithers (12 files, 2 dithers of each band)
    input_models = []
    input_models.append(input_model1)
    input_models.append(input_model1)
    input_models.append(input_model2)
    input_models.append(input_model2)
    input_models.append(input_model3)
    input_models.append(input_model3)
    input_models.append(input_model4)
    input_models.append(input_model4)
    input_models.append(input_model5)
    input_models.append(input_model5)
    input_models.append(input_model6)
    input_models.append(input_model6)
    return input_models


@pytest.fixture(scope='function')
def nirspec_medium_coverage():
    """ Generate a IFU images NIRSpec G140M, G235M """

    input_model1 = datamodels.IFUImageModel()
    input_model1.meta.wcsinfo._instance.update(wcsinfo)
    input_model1.meta.instrument._instance.update(nirspec_G140M)
    input_model1.meta.observation._instance.update(observation)
    input_model1.meta.subarray._instance.update(subarray_nirspec)
    input_model1.meta.cal_step.assign_wcs = 'COMPLETE'

    input_model2 = datamodels.IFUImageModel()
    input_model2.meta.wcsinfo._instance.update(wcsinfo)
    input_model2.meta.instrument._instance.update(nirspec_G235M)
    input_model2.meta.observation._instance.update(observation)
    input_model2.meta.subarray._instance.update(subarray_nirspec)
    input_model2.meta.cal_step.assign_wcs = 'COMPLETE'

    input_models = []
    input_models.append(input_model1)
    input_models.append(input_model1)
    input_models.append(input_model2)
    input_models.append(input_model2)
    return input_models


def test_calspec2_config(_jail, miri_ifushort_short):
    """ Determine cube based on calspec2 setup """

    pars_input = {}
    pars_input['channel'] = []
    pars_input['subchannel'] = []
    pars_input['filter'] = []
    pars_input['grating'] = []
    weighting = 'msm'
    output_type = 'multi'  # calspec 2 setup. Only 1 cube create from 2 channels
    single = False
    par_filename = 'None'

    input_file = 'test.fits'
    input_models = []
    input_filenames = []
    input_models.append(miri_ifushort_short)
    input_filenames.append(input_file)

    pars = {
        'channel': pars_input['channel'],
        'subchannel': pars_input['subchannel'],
        'grating': pars_input['grating'],
        'filter': pars_input['filter'],
        'weighting': weighting,
        'single': single,
        'output_type': output_type}

    cubeinfo = cube_build.CubeData(
        input_models,
        input_filenames,
        par_filename,
        **pars)

    master_table = file_table.FileTable()
    this_instrument = master_table.set_file_table(
        cubeinfo.input_models, cubeinfo.input_filenames)

    assert this_instrument == 'MIRI'

    cubeinfo.instrument = this_instrument
    cubeinfo.determine_band_coverage(master_table)
    assert cubeinfo.all_channel == ['1', '2']
    assert cubeinfo.all_subchannel == ['short', 'short']

    num_cubes, cube_pars = cubeinfo.number_cubes()
    assert num_cubes == 1
    assert cube_pars['1']['par1'] == ['1', '2']
    assert cube_pars['1']['par2'] == ['short', 'short']


def test_calspec3_config_miri(_jail, miri_full_coverage):
    """ Test CalSpec3 MIRI configuration default band cubes"""

    pars_input = {}
    pars_input['channel'] = []
    pars_input['subchannel'] = []
    pars_input['filter'] = []
    pars_input['grating'] = []
    weighting = 'msm'
    output_type = 'band'
    single = False
    par_filename = 'None'

    input_file = 'test.fits'
    num_files = len(miri_full_coverage)
    input_filenames = []
    for i in range(num_files):
        input_filenames.append(input_file)

    pars = {
        'channel': pars_input['channel'],
        'subchannel': pars_input['subchannel'],
        'grating': pars_input['grating'],
        'filter': pars_input['filter'],
        'weighting': weighting,
        'single': single,
        'output_type': output_type}

    cubeinfo = cube_build.CubeData(
        miri_full_coverage,
        input_filenames,
        par_filename,
        **pars)

    master_table = file_table.FileTable()
    this_instrument = master_table.set_file_table(
        cubeinfo.input_models, cubeinfo.input_filenames)

    assert this_instrument == 'MIRI'

    cubeinfo.instrument = this_instrument
    cubeinfo.determine_band_coverage(master_table)
    num_cubes, cube_pars = cubeinfo.number_cubes()
    assert num_cubes == 12
    assert cubeinfo.all_channel == ['1', '1', '1', '2', '2', '2',
                                    '3', '3', '3', '4', '4', '4']
    assert cubeinfo.all_subchannel == ['short', 'medium', 'long',
                                       'short', 'medium', 'long',
                                       'short', 'medium', 'long',
                                       'short', 'medium', 'long']

    assert cube_pars['1']['par1'] == ['1']
    assert cube_pars['1']['par2'] == ['short']
    assert cube_pars['2']['par1'] == ['1']
    assert cube_pars['2']['par2'] == ['medium']
    assert cube_pars['3']['par1'] == ['1']
    assert cube_pars['3']['par2'] == ['long']

    assert cube_pars['4']['par1'] == ['2']
    assert cube_pars['4']['par2'] == ['short']
    assert cube_pars['5']['par1'] == ['2']
    assert cube_pars['5']['par2'] == ['medium']
    assert cube_pars['6']['par1'] == ['2']
    assert cube_pars['6']['par2'] == ['long']

    assert cube_pars['7']['par1'] == ['3']
    assert cube_pars['7']['par2'] == ['short']
    assert cube_pars['8']['par1'] == ['3']
    assert cube_pars['8']['par2'] == ['medium']
    assert cube_pars['9']['par1'] == ['3']
    assert cube_pars['9']['par2'] == ['long']

    assert cube_pars['10']['par1'] == ['4']
    assert cube_pars['10']['par2'] == ['short']
    assert cube_pars['11']['par1'] == ['4']
    assert cube_pars['11']['par2'] == ['medium']
    assert cube_pars['12']['par1'] == ['4']
    assert cube_pars['12']['par2'] == ['long']


def test_calspec3_config_miri_multi(_jail, miri_full_coverage):
    """ Test CalSpec3 MIRI configuration default band cubes"""

    pars_input = {}
    pars_input['channel'] = []
    pars_input['subchannel'] = []
    pars_input['filter'] = []
    pars_input['grating'] = []
    weighting = 'msm'
    output_type = 'multi'
    single = False
    par_filename = 'None'

    input_file = 'test.fits'
    num_files = len(miri_full_coverage)
    input_filenames = []
    for i in range(num_files):
        input_filenames.append(input_file)

    pars = {
        'channel': pars_input['channel'],
        'subchannel': pars_input['subchannel'],
        'grating': pars_input['grating'],
        'filter': pars_input['filter'],
        'weighting': weighting,
        'single': single,
        'output_type': output_type}

    cubeinfo = cube_build.CubeData(
        miri_full_coverage,
        input_filenames,
        par_filename,
        **pars)

    master_table = file_table.FileTable()
    this_instrument = master_table.set_file_table(
        cubeinfo.input_models, cubeinfo.input_filenames)

    assert this_instrument == 'MIRI'

    cubeinfo.instrument = this_instrument
    cubeinfo.determine_band_coverage(master_table)
    num_cubes, cube_pars = cubeinfo.number_cubes()
    assert num_cubes == 1
    assert cubeinfo.all_channel == ['1', '1', '1', '2', '2', '2',
                                    '3', '3', '3', '4', '4', '4']
    assert cubeinfo.all_subchannel == ['short', 'medium', 'long',
                                       'short', 'medium', 'long',
                                       'short', 'medium', 'long',
                                       'short', 'medium', 'long']

    assert cube_pars['1']['par1'] == ['1', '1', '1',
                                      '2', '2', '2',
                                      '3', '3', '3',
                                      '4', '4', '4']
    assert cube_pars['1']['par2'] == ['short', 'medium', 'long',
                                      'short', 'medium', 'long',
                                      'short', 'medium', 'long',
                                      'short', 'medium', 'long']


def test_calspec3_config_nirspec(_jail, nirspec_medium_coverage):
    """ Test CalSpec3 configuration for NIRSpec - default band cubes"""

    pars_input = {}
    pars_input['channel'] = []
    pars_input['subchannel'] = []
    pars_input['filter'] = []
    pars_input['grating'] = []
    weighting = 'msm'
    output_type = 'band'
    single = False
    par_filename = 'None'

    input_file = 'test.fits'
    num_files = len(nirspec_medium_coverage)
    input_filenames = []
    for i in range(num_files):
        input_filenames.append(input_file)

    pars = {
        'channel': pars_input['channel'],
        'subchannel': pars_input['subchannel'],
        'grating': pars_input['grating'],
        'filter': pars_input['filter'],
        'weighting': weighting,
        'single': single,
        'output_type': output_type}

    cubeinfo = cube_build.CubeData(
        nirspec_medium_coverage,
        input_filenames,
        par_filename,
        **pars)

    master_table = file_table.FileTable()
    this_instrument = master_table.set_file_table(
        cubeinfo.input_models, cubeinfo.input_filenames)

    assert this_instrument == 'NIRSPEC'

    cubeinfo.instrument = this_instrument
    cubeinfo.determine_band_coverage(master_table)
    num_cubes, cube_pars = cubeinfo.number_cubes()

    assert num_cubes == 2

    assert cubeinfo.all_grating == ['g140m', 'g235m']
    assert cubeinfo.all_filter == ['f100lp', 'f170lp']

    assert cube_pars['1']['par1'] == ['g140m']
    assert cube_pars['1']['par2'] == ['f100lp']
    assert cube_pars['2']['par1'] == ['g235m']
    assert cube_pars['2']['par2'] == ['f170lp']


def test_calspec3_config_nirspec_multi(_jail, nirspec_medium_coverage):
    """ Test CalSpec3 configuration for NIRSpec - Multiband cubes"""

    pars_input = {}
    pars_input['channel'] = []
    pars_input['subchannel'] = []
    pars_input['filter'] = []
    pars_input['grating'] = []
    weighting = 'msm'
    output_type = 'multi'
    single = False
    par_filename = 'None'

    input_file = 'test.fits'
    num_files = len(nirspec_medium_coverage)
    input_filenames = []
    for i in range(num_files):
        input_filenames.append(input_file)

    pars = {
        'channel': pars_input['channel'],
        'subchannel': pars_input['subchannel'],
        'grating': pars_input['grating'],
        'filter': pars_input['filter'],
        'weighting': weighting,
        'single': single,
        'output_type': output_type}

    cubeinfo = cube_build.CubeData(
        nirspec_medium_coverage,
        input_filenames,
        par_filename,
        **pars)

    master_table = file_table.FileTable()
    this_instrument = master_table.set_file_table(
        cubeinfo.input_models, cubeinfo.input_filenames)

    assert this_instrument == 'NIRSPEC'

    cubeinfo.instrument = this_instrument
    cubeinfo.determine_band_coverage(master_table)
    num_cubes, cube_pars = cubeinfo.number_cubes()

    assert num_cubes == 1
    assert cubeinfo.all_grating == ['g140m', 'g235m']
    assert cubeinfo.all_filter == ['f100lp', 'f170lp']

    assert cube_pars['1']['par1'] == ['g140m', 'g235m']
    assert cube_pars['1']['par2'] == ['f100lp', 'f170lp']
