import time
import logging
import numpy as np
import warnings
#import multiprocessing
from multiprocessing.pool import Pool as Pool
import dill
#from pathos.multiprocessing import ProcessingPool as Pool
import psutil
import astropy
import pickle


from .. import datamodels
from ..datamodels import dqflags
from ..lib import pipe_utils
from . import ramp_fit as rf

from . import gls_fit           # used only if algorithm is "GLS"
from . import utils
log = logging.getLogger(__name__)

def ols_ramp_fit_multi(input_model, buffsize, save_opt, readnoise_2d, gain_2d,
                 weighting):
    total_cores = psutil.cpu_count(logical=False)
    #   nrows = int(cubeshape[1]/numslice)
    #   numslice = 1
    model_file = open('file.tst', 'ab')
#    pickle._dump(input_model,model_file)
    number_slices = min(8, total_cores)
    log.info("number of processes being used is %d" % number_slices)
    print("total cores ",total_cores)
    print("number_slices ",number_slices)

    total_rows = input_model.data.shape[2]
    total_cols = input_model.data.shape[3]
    number_of_integrations = input_model.data.shape[0]
    number_of_groups = input_model.data.shape[1]
    data = np.zeros_like(input_model.data)
    err = np.zeros_like(input_model.err)
    groupdq = np.zeros_like(input_model.groupdq)
    pixeldq = np.zeros_like(input_model.pixeldq)
    rows_per_slice = round(total_rows / number_slices)
    pool = Pool(processes=number_slices)
    slices = []
    for i in range(number_slices - 1):

        readnoise_slice = readnoise_2d[i * rows_per_slice: (i + 1) * rows_per_slice, :]
        gain_slice = gain_2d[i * rows_per_slice: (i + 1) * rows_per_slice, :]
        data = input_model.data[:,:,i * rows_per_slice: (i + 1) * rows_per_slice, :].copy()
        err = input_model.err[:, :, i * rows_per_slice: (i + 1) * rows_per_slice, :].copy()
        groupdq = input_model.groupdq[:, :, i * rows_per_slice: (i + 1) * rows_per_slice, :].copy()
        pixeldq = input_model.pixeldq[ i * rows_per_slice: (i + 1) * rows_per_slice, :].copy()
        if pipe_utils.is_tso(input_model) and hasattr(input_model, 'int_times'):
            int_times = input_model.int_times
        else:
            int_times = None
        slices.insert(i, (data, err, groupdq, pixeldq, buffsize, save_opt, readnoise_slice, gain_slice, weighting,
                      input_model.meta.instrument.name, input_model.meta.exposure.frame_time,
                      input_model.meta.exposure.ngroups, input_model.meta.exposure.group_time,
                      input_model.meta.exposure.groupgap, input_model.meta.exposure.nframes,
                      input_model.meta.exposure.drop_frames1, int_times))

    # last slice gets the rest
    readnoise_slice = readnoise_2d[(number_slices - 1) * rows_per_slice: total_rows, :]
    gain_slice = gain_2d[(number_slices - 1) * rows_per_slice: total_rows, :]
    data = input_model.data[:, :, (number_slices - 1) * rows_per_slice: total_rows, :].copy()
    err = input_model.err[:, :, (number_slices - 1) * rows_per_slice: total_rows, :].copy()
    groupdq = input_model.groupdq[:, :, (number_slices - 1) * rows_per_slice: total_rows, :].copy()
    pixeldq = input_model.pixeldq[(number_slices - 1) * rows_per_slice: total_rows, :].copy()
    slices.insert(number_slices - 1, (data, err, groupdq, pixeldq, buffsize, save_opt, readnoise_slice, gain_slice, weighting,
                                      input_model.meta.instrument.name, input_model.meta.exposure.frame_time,
                                      input_model.meta.exposure.ngroups, input_model.meta.exposure.group_time,
                                      input_model.meta.exposure.groupgap, input_model.meta.exposure.nframes,
                                      input_model.meta.exposure.drop_frames1, int_times))

    log.info("Creating %d processes for ramp fitting " % number_slices)
  #  pool = multiprocessing.Pool(processes=number_slices)
  #  real_results = pool.map(rf.ols_ramp_fit, slices[0], slices[1], slices[2], slices[3])
    real_results = pool.starmap(rf.ols_ramp_fit, slices)
    k = 0
    log.info("All processes complete")
    # Create new model for the primary output.
    imshape = (total_rows, total_cols)
    out_model = datamodels.ImageModel(data=np.zeros(imshape, dtype=np.float32),
                                      dq=np.zeros(imshape, dtype=np.uint32),
                                      var_poisson=np.zeros(imshape, dtype=np.float32),
                                      var_rnoise=np.zeros(imshape, dtype=np.float32),
                                      err=np.zeros(imshape, dtype=np.float32))
    out_model.update(input_model)  # ... and add all keys from input
    #create per integrations model, if this is a multi-integration exposure
    if number_of_integrations > 1:
        int_model = datamodels.CubeModel(
            data=np.zeros((number_of_integrations,) + imshape, dtype=np.float32),
            dq=np.zeros(imshape, dtype=np.unit32),
            var_poisson=np.zeros((number_of_integrations,) + imshape, dtype=np.float32),
            var_rnoise=np.zeros((number_of_integrations,) + imshape, dtype=np.float32),
            int_times=np.zeros((number_of_integrations,) + imshape, dtype=np.float32),
            err=np.zeros((number_of_integrations,) + imshape, dtype=np.float32))
        int_model.update(input_model)  # ... and add all keys from input
    else:
        int_model = None

    # Create model for the optional output
    if save_opt:
        opt_model = datamodels.RampFitOutputModel(
                slope=np.zeros((number_of_integrations,) + (1,) + imshape, dtype=np.float32),
                sigslope=np.zeros((number_of_integrations,) + (1,) + imshape, dtype=np.float32),
                var_poisson=np.zeros((number_of_integrations,) + (1,) + imshape, dtype=np.float32),
                var_rnoise=np.zeros((number_of_integrations,) + (1,) + imshape, dtype=np.float32),
                yint=np.zeros((number_of_integrations,) + (1,) + imshape, dtype=np.float32),
                sigyint=np.zeros((number_of_integrations,) + (1,) + imshape, dtype=np.float32),
                pedestal=np.zeros((number_of_integrations,) + imshape, dtype=np.float32),
                weights=np.zeros((number_of_integrations,) + (1,) + imshape, dtype=np.float32),
                crmag=np.zeros((number_of_integrations,) + (1,) + imshape, dtype=np.float32))

        opt_model.meta.filename = input_model.meta.filename
        opt_model.update(input_model)  # ... and add all keys from input
    else:
        opt_model = None

    for resultslice in real_results:
        if len(real_results) == k + 1:  # last result
            out_model.data[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[0]
            out_model.dq[k * rows_per_slice:total_rows, :] = resultslice[1]
            out_model.var_poisson[:, :, k * rows_per_slice:total_rows, :] = resultslice[2]
            out_model.var_rnoise[:, :, k * rows_per_slice:total_rows, :] = resultslice[3]
            out_model.err[:, :, k * rows_per_slice:total_rows, :] = resultslice[4]
            if resultslice[5] is not None:
                int_model.data[:, :, k * rows_per_slice:total_rows, :] = resultslice[5]
                int_model.dq[k * rows_per_slice:total_rows, :] = resultslice[6]
                int_model.var_poisson[:, :, k * rows_per_slice:total_rows, :] = resultslice[7]
                int_model.var_rnoise[:, :, k * rows_per_slice:total_rows, :] = resultslice[8]
                int_model.err[:, :, k * rows_per_slice:total_rows, :] = resultslice[9]
                int_model.int_times[:, :, k * rows_per_slice:total_rows, :] = resultslice[10]
            if resultslice[11] is not None:
                opt_model.slope[:, :, k * rows_per_slice:total_rows, :] = resultslice[11]
                opt_model.sigslope[:, :, k * rows_per_slice:total_rows, :] = resultslice[12]
                opt_model.var_poisson[:, :, k * rows_per_slice:total_rows, :] = resultslice[13]
                opt_model.var_rnoise[:, :, k * rows_per_slice:total_rows, :] = resultslice[14]
                opt_model.yint[:, :, k * rows_per_slice:total_rows, :] = resultslice[15]
                opt_model.sigyint[:, :, k * rows_per_slice:total_rows, :] = resultslice[16]
                opt_model.pedestal[:, :, k * rows_per_slice:total_rows, :] = resultslice[17]
                opt_model.weights[:, :, k * rows_per_slice:total_rows, :] = resultslice[18]
                opt_model.crmag[:, :, k * rows_per_slice:total_rows, :] = resultslice[19]
        else:
            out_model.data[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[0]
            out_model.dq[k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[1]
            out_model.var_poisson[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[2]
            out_model.var_rnoise[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[3]
            out_model.err[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[4]
            if resultslice[5] is not None:
                int_model.data[:, :, k * rows_per_slice: (k + 1) * rows_per_slice, :] = resultslice[5]
                int_model.dq[k * rows_per_slice: (k + 1) * rows_per_slice, :] = resultslice[6]
                int_model.var_poisson[:, :, k * rows_per_slice: (k + 1) * rows_per_slice, :] = resultslice[7]
                int_model.var_rnoise[:, :, k * rows_per_slice: (k + 1) * rows_per_slice, :] = resultslice[8]
                int_model.err[:, :, k * rows_per_slice: (k + 1) * rows_per_slice, :] = resultslice[9]
                int_model.int_times[:, :, k * rows_per_slice: (k + 1) * rows_per_slice, :] = resultslice[10]
            if resultslice[11] is not None:
                opt_model.slope[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[11]
                opt_model.sigslope[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[12]
                opt_model.var_poisson[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[13]
                opt_model.var_rnoise[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[14]
                opt_model.yint[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[15]
                opt_model.sigyint[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[16]
                opt_model.pedestal[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[17]
                opt_model.weights[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[18]
                opt_model.crmag[:, :, k * rows_per_slice: (k + 1) *rows_per_slice, :] = resultslice[19]
        k = k + 1
    return out_model, int_model, opt_model
