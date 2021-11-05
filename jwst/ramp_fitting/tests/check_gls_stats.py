# program to check GLS cpu and memory usage

import time
import numpy as np

from jwst.datamodels import dqflags
from jwst.datamodels import RampModel
from jwst.datamodels import GainModel, ReadnoiseModel

from stcal.ramp_fitting import ramp_fit

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)


def setup_inputs(ngroups=10, readnoise=10, nints=1,
                 nrows=1032, ncols=1024, nframes=1, grouptime=1.0, gain=1, deltatime=1):
    times = np.array(list(range(ngroups)), dtype=np.float64) * deltatime
    gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
    err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    pixdq = np.zeros(shape=(nrows, ncols), dtype=np.float64)
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
    model1 = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
    model1.meta.instrument.name = 'MIRI'
    model1.meta.instrument.detector = 'MIRIMAGE'
    model1.meta.instrument.filter = 'F480M'
    model1.meta.observation.date = '2015-10-13'
    model1.meta.exposure.type = 'MIR_IMAGE'
    model1.meta.exposure.group_time = deltatime
    model1.meta.subarray.name = 'FULL'
    model1.meta.subarray.xstart = 1
    model1.meta.subarray.ystart = 1
    model1.meta.subarray.xsize = 1024
    model1.meta.subarray.ysize = 1032
    model1.meta.exposure.frame_time = deltatime
    model1.meta.exposure.ngroups = ngroups
    model1.meta.exposure.group_time = deltatime
    model1.meta.exposure.nframes = 1
    model1.meta.exposure.groupgap = 0
    gain = GainModel(data=gain)
    gain.meta.instrument.name = 'MIRI'
    gain.meta.subarray.xstart = 1
    gain.meta.subarray.ystart = 1
    gain.meta.subarray.xsize = 1024
    gain.meta.subarray.ysize = 1032
    rnModel = ReadnoiseModel(data=read_noise)
    rnModel.meta.instrument.name = 'MIRI'
    rnModel.meta.subarray.xstart = 1
    rnModel.meta.subarray.ystart = 1
    rnModel.meta.subarray.xsize = 1024
    rnModel.meta.subarray.ysize = 1032
    return model1, gdq, rnModel, pixdq, err, gain


def simple_ramp(fit_method='OLS', ngroups=10, cr_group=None):
    # Here given a 10 group ramp with an exact slope of 20/group.
    # The output slope should be 20.
    slope_per_group = 20.
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups)
    for k in range(ngroups):
        model1.data[0, k, :, :] = 10.0 + k * slope_per_group

    if cr_group:
        model1.groupdq[0, cr_group, :, :] = dqflags.group['JUMP_DET']

    slopes = ramp_fit(model1, 1024 * 1., False, rnModel, gain, fit_method,
                      'optimal')
    return slopes


if __name__ == '__main__':

    solve_method = ['GLS', 'OLS']
    ngroups = [10, 20, 40]
    crgroup = [None, 4]
    for cmethod in solve_method:
        for cgroups in ngroups:
            for cur_cr in crgroup:
                start_time = time.process_time()
                slopes = simple_ramp(fit_method=cmethod,
                                     ngroups=cgroups,
                                     cr_group=cur_cr)
                end_time = time.process_time()
                print('timing: ', cmethod, cgroups, cur_cr, end_time - start_time)
