import os
import pytest
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_nrs_fs_multi_spec2_2():
    """

    Regression test of calwebb_spec2 pipeline performed on NIRSpec fixed-slit data.

    """
    step = Spec2Pipeline()
    step.save_bsub = True
    step.save_results = True
    step.resample_spec.save_results = True
    step.cube_build.save_results = True
    step.extract_1d.save_results = True
    step.run(BIGDATA+'/pipelines/jwtest1013001_01101_00001_NRS1_rate.fits')

    na = 'jwtest1013001_01101_00001_NRS1_cal.fits'
    nb = BIGDATA+'/pipelines/jwtest1013001_01101_00001_NRS1_cal_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci',1],h['err',1],h['dq',1],h['relsens',1],h['wavelength',1],
                                    h['pathloss_pointsource',1],h['wavelength_pointsource',1],
                                    h['pathloss_uniformsource',1],h['wavelength_uniformsource',1],
                                    h['sci',2],h['err',2],h['dq',2],h['relsens',2],h['wavelength',2],
                                    h['pathloss_pointsource',2],h['wavelength_pointsource',2],
                                    h['pathloss_uniformsource',2],h['wavelength_uniformsource',2],
                                    h['sci',3],h['err',3],h['dq',3],h['relsens',3],h['wavelength',3],
                                    h['sci',4],h['err',4],h['dq',4],h['relsens',4],h['wavelength',4],
                                    h['pathloss_pointsource',4],h['wavelength_pointsource',4],
                                    h['pathloss_uniformsource',4],h['wavelength_uniformsource',4]])
    newhref = pf.HDUList([href['primary'],href['sci',1],href['err',1],href['dq',1],href['relsens',1],href['wavelength',1],
                                          href['pathloss_pointsource',1],href['wavelength_pointsource',1],
                                          href['pathloss_uniformsource',1],href['wavelength_uniformsource',1],
                                          href['sci',2],href['err',2],href['dq',2],href['relsens',2],href['wavelength',2],
                                          href['pathloss_pointsource',2],href['wavelength_pointsource',2],
                                          href['pathloss_uniformsource',2],href['wavelength_uniformsource',2],
                                          href['sci',3],href['err',3],href['dq',3],href['relsens',3],href['wavelength',3],
                                          href['sci',4],href['err',4],href['dq',4],href['relsens',4],href['wavelength',4],
                                          href['pathloss_pointsource',4],href['wavelength_pointsource',4],
                                          href['pathloss_uniformsource',4],href['wavelength_uniformsource',4]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    na = 'jwtest1013001_01101_00001_NRS1_s2d.fits'
    nb = BIGDATA+'/pipelines/jwtest1013001_01101_00001_NRS1_s2d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci',1],h['wht',1],h['con',1],h['relsens',1],
                                    h['sci',2],h['wht',2],h['con',2],h['relsens',2],
                                    h['sci',3],h['wht',3],h['con',3],h['relsens',3],
                                    h['sci',4],h['wht',4],h['con',4],h['relsens',4]])
    newhref = pf.HDUList([href['primary'],href['sci',1],href['wht',1],href['con',1],href['relsens',1],
                                          href['sci',2],href['wht',2],href['con',2],href['relsens',2],
                                          href['sci',3],href['wht',3],href['con',3],href['relsens',3],
                                          href['sci',4],href['wht',4],href['con',4],href['relsens',4]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    na = 'jwtest1013001_01101_00001_NRS1_x1d.fits'
    nb = BIGDATA+'/pipelines/jwtest1013001_01101_00001_NRS1_x1d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['extract1d',1],h['extract1d',2],h['extract1d',3],h['extract1d',4]])
    newhref = pf.HDUList([href['primary'],href['extract1d',1],href['extract1d',2],href['extract1d',3],href['extract1d',4]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

