import pytest
import numpy as np
from astropy.io import fits
from jwst.jump.jump import detect_jumps
from jwst.datamodels import dqflags
from jwst.datamodels import MIRIRampModel
from jwst.datamodels import GainModel, ReadnoiseModel
from jwst.jump import JumpStep
from itertools import cycle

def test_one_core():
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 100
    CR_fraction = 2
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1,89) if x % 5 == 0]
    CR_locs = [x for x in range(1032*1024) if x % CR_fraction == 0]
    CR_x_locs = [x % 1032 for x in CR_locs]
    CR_y_locs = [np.int(x / 1032) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0,CR_group:,CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0,CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    print("number of CRs "+ len(CR_x_locs).__str__())
    gain = np.ones(shape=(1024, 1032), dtype=np.float64) * ingain
    hdr = fits.Header()
    hdr['INSTRUME'] = 'MIRI'
    hdr['SUBARRAY'] = 'FULL'
    hdr['SUBSTRT1'] = 1
    hdr['SUBSIZE1'] = 1032
    hdr['SUBSTRT2'] = 1
    hdr['SUBSIZE2'] = 1024
    hdu = fits.PrimaryHDU(gain,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=gain))
    hdul.writeto('gain.fits', overwrite=True)

    rnoise = np.ones(shape=(1024, 1032), dtype=np.float64) * inreadnoise
    hdu = fits.PrimaryHDU(rnoise,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=rnoise))
    hdul.writeto('readnoise.fits', overwrite=True)


    out_model = JumpStep.call(model1, override_gain='gain.fits', override_readnoise = 'readnoise.fits',
                              maximum_cores='one')
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        print,CR_group
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))

def test_one_core_NIRCAM():
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 100
    CR_fraction = 2
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,nrows=2048, ncols=2048,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1,89) if x % 5 == 0]
    CR_locs = [x for x in range(2048*2048) if x % CR_fraction == 0]
    CR_x_locs = [x % 2048 for x in CR_locs]
    CR_y_locs = [np.int(x / 2048) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0,CR_group:,CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0,CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    print("number of CRs "+ len(CR_x_locs).__str__())
    gain = np.ones(shape=(2048, 2048), dtype=np.float64) * ingain
    hdr = fits.Header()
    hdr['INSTRUME'] = 'MIRI'
    hdr['SUBARRAY'] = 'FULL'
    hdr['SUBSTRT1'] = 1
    hdr['SUBSIZE1'] = 2048
    hdr['SUBSTRT2'] = 1
    hdr['SUBSIZE2'] = 2048
    hdu = fits.PrimaryHDU(gain,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=gain))
    hdul.writeto('gain.fits', overwrite=True)

    rnoise = np.ones(shape=(2048, 2048), dtype=np.float64) * inreadnoise
    hdu = fits.PrimaryHDU(rnoise,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=rnoise))
    hdul.writeto('readnoise.fits', overwrite=True)


    out_model = JumpStep.call(model1, override_gain='gain.fits', override_readnoise = 'readnoise.fits',
                              maximum_cores='one')
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        print,CR_group
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))

def test_all_cores_NIRCAM():
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 100
    CR_fraction = 2
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,nrows=2048, ncols=2048,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1,89) if x % 5 == 0]
    CR_locs = [x for x in range(2048*2048) if x % CR_fraction == 0]
    CR_x_locs = [x % 2048 for x in CR_locs]
    CR_y_locs = [np.int(x / 2048) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0,CR_group:,CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0,CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    print("number of CRs "+ len(CR_x_locs).__str__())
    gain = np.ones(shape=(2048, 2048), dtype=np.float64) * ingain
    hdr = fits.Header()
    hdr['INSTRUME'] = 'MIRI'
    hdr['SUBARRAY'] = 'FULL'
    hdr['SUBSTRT1'] = 1
    hdr['SUBSIZE1'] = 2048
    hdr['SUBSTRT2'] = 1
    hdr['SUBSIZE2'] = 2048
    hdu = fits.PrimaryHDU(gain,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=gain))
    hdul.writeto('gain.fits', overwrite=True)

    rnoise = np.ones(shape=(2048, 2048), dtype=np.float64) * inreadnoise
    hdu = fits.PrimaryHDU(rnoise,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=rnoise))
    hdul.writeto('readnoise.fits', overwrite=True)


    out_model = JumpStep.call(model1, override_gain='gain.fits', override_readnoise = 'readnoise.fits',
                              maximum_cores='all')
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        print,CR_group
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))

def test_half_cores():
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 100
    CR_fraction = 2
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1,89) if x % 5 == 0]
    CR_locs = [x for x in range(1032*1024) if x % CR_fraction == 0]
    CR_x_locs = [x % 1032 for x in CR_locs]
    CR_y_locs = [np.int(x / 1032) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0,CR_group:,CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0,CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    print("number of CRs "+ len(CR_x_locs).__str__())
    gain = np.ones(shape=(1024, 1032), dtype=np.float64) * ingain
    hdr = fits.Header()
    hdr['INSTRUME'] = 'MIRI'
    hdr['SUBARRAY'] = 'FULL'
    hdr['SUBSTRT1'] = 1
    hdr['SUBSIZE1'] = 1032
    hdr['SUBSTRT2'] = 1
    hdr['SUBSIZE2'] = 1024
    hdu = fits.PrimaryHDU(gain,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=gain))
    hdul.writeto('gain.fits', overwrite=True)

    rnoise = np.ones(shape=(1024, 1032), dtype=np.float64) * inreadnoise
    hdu = fits.PrimaryHDU(rnoise,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=rnoise))
    hdul.writeto('readnoise.fits', overwrite=True)


    out_model = JumpStep.call(model1, override_gain='gain.fits', override_readnoise = 'readnoise.fits',
                              maximum_cores='all')
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        print,CR_group
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))

def test_one_core_two_CRs():
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 100
    CR_fraction = 3
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1,89) if x % 5 == 0]
    CR_locs = [x for x in range(1032*1024) if x % CR_fraction == 0]
    CR_x_locs = [x % 1032 for x in CR_locs]
    CR_y_locs = [np.int(x / 1032) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0,CR_group:,CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0,CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500
        model1.data[0, CR_group+8:, CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0, CR_group+8:, CR_y_locs[i], CR_x_locs[i]] + 700


    print("number of CRs "+ len(CR_x_locs).__str__())
    gain = np.ones(shape=(1024, 1032), dtype=np.float64) * ingain
    hdr = fits.Header()
    hdr['INSTRUME'] = 'MIRI'
    hdr['SUBARRAY'] = 'FULL'
    hdr['SUBSTRT1'] = 1
    hdr['SUBSIZE1'] = 1032
    hdr['SUBSTRT2'] = 1
    hdr['SUBSIZE2'] = 1024
    hdu = fits.PrimaryHDU(gain,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=gain))
    hdul.writeto('gain.fits', overwrite=True)

    rnoise = np.ones(shape=(1024, 1032), dtype=np.float64) * inreadnoise
    hdu = fits.PrimaryHDU(rnoise,header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(name="SCI", data=rnoise))
    hdul.writeto('readnoise.fits', overwrite=True)


    out_model = JumpStep.call(model1, override_gain='gain.fits', override_readnoise = 'readnoise.fits',
                              maximum_cores='one')
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))
        assert (4 == np.max(out_model.groupdq[0, CR_group+8, CR_y_locs[i], CR_x_locs[i]]))



def setup_inputs(ngroups=10, readnoise=10, nints=1,
                 nrows=1024, ncols=1032, nframes=1, grouptime=1.0, gain=1, deltatime=1):
    print('readnoise', readnoise)
    print('gain', gain)
    times = np.array(list(range(ngroups)), dtype=np.float64) * deltatime
    gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
    err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    pixdq = np.zeros(shape=(nrows, ncols), dtype=np.float64)
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
    model1 = MIRIRampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
    model1.meta.instrument.name = 'MIRI'
    model1.meta.instrument.detector = 'MIRIMAGE'
    model1.meta.instrument.filter = 'F480M'
    model1.meta.observation.date = '2015-10-13'
    model1.meta.exposure.type = 'MIR_IMAGE'
    model1.meta.exposure.group_time = deltatime
    model1.meta.subarray.name = 'FULL'
    model1.meta.subarray.xstart = 1
    model1.meta.subarray.ystart = 1
    model1.meta.subarray.xsize = 1032
    model1.meta.subarray.ysize = 1024
    model1.meta.exposure.frame_time = deltatime
    model1.meta.exposure.ngroups = ngroups
    model1.meta.exposure.group_time = deltatime
    model1.meta.exposure.nframes = 1
    model1.meta.exposure.groupgap = 0
    gain = GainModel(data=gain)
    gain.meta.instrument.name = 'MIRI'
    gain.meta.subarray.xstart = 1
    gain.meta.subarray.ystart = 1
    gain.meta.subarray.xsize = 1032
    gain.meta.subarray.ysize = 1024
    rnModel = ReadnoiseModel(data=read_noise)
    rnModel.meta.instrument.name = 'MIRI'
    rnModel.meta.subarray.xstart = 1
    rnModel.meta.subarray.ystart = 1
    rnModel.meta.subarray.xsize = 1032
    rnModel.meta.subarray.ysize = 1024
    return model1, gdq, rnModel, pixdq, err, gain