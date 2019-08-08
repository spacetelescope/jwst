""" 
Unit tests for straylight correction
"""

from jwst.datamodels import IFUImageModel
from jwst.straylight import StraylightStep 
from jwst.straylight.straylight import correct_mrs_modshepard, shepard_2d_kernel
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
import numpy as np
import pytest


def test_correct_mrs_modshepard():
    
    image = IFUImageModel((50,50))
    image.data = np.ones((50,50)) + 30.0

    slice_map = np.ones((50,50))
    # create 2 slice gaps

    slice_map[:,10:15] = 0
    slice_map[:,30:35] = 0
    image.data[:,10:15] = 0.5
    image.data[:,30:35] = 0.5
    roi = 4
    power = 1

    result = correct_mrs_modshepard(image, slice_map, roi, power)
    assert type(image) is type(result)

def test_correct_detector():
    image = IFUImageModel((20,20))
    image.data = np.random.random((20,20))
    image.meta.instrument.name  = 'MIRI'
    image.meta.instrument.detector = 'MIRIFULONG'
    image.data = np.random.random((50,50))
    collect_pipeline_cfgs('./config')
    result = StraylightStep.call(image,
                                 config_file='config/straylight.cfg')
    assert result.meta.cal_step.straylight == 'SKIPPED'
    assert type(image) is type(result)


def test_not_nirspec():
    image = IFUImageModel((20,20))
    image.data = np.random.random((20,20))
    image.meta.instrument.name  = 'NIRSPEC'
    image.meta.instrument.detector = 'NRS1'
    collect_pipeline_cfgs('./config')
    result = StraylightStep.call(image,
                                 config_file='config/straylight.cfg')
    assert result.meta.cal_step.straylight == 'SKIPPED'
    assert type(image) is type(result)


def test_roi_parameter1():
    image = IFUImageModel((1032,1024))
    collect_pipeline_cfgs('./config')
    image.meta.instrument.name  = 'MIRI'
    image.meta.observation.date = '2018-09-07'
    image.meta.observation.time = '10:32:20.181'
    image.meta.instrument.detector = 'MIRIFUSHORT'
    image.meta.instrument.channel = '12'
    image.meta.instrument.band = 'SHORT'
    image.meta.exposure.type = 'MIR_MRS'
    result = StraylightStep.call(image, 
                                 config_file = 'config/straylight.cfg',
                                 method='ModShepard',
                                 roi=5, # should be even value
                                 power=2)
    assert result.meta.cal_step.straylight == 'SKIPPED'
    assert type(image) is type(result)
