# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test jwst.transforms
"""
import pytest
from numpy.testing.utils import assert_allclose
from ..import models


#_RANDOM_SEED = 0x1337

"""
def test_logical():

    with NumpyRNGContext(0x1337):
        compareto = np.random.randn(10)
    with NumpyRNGContext(0x1338):
        val = np.random.randn(10)
    with NumpyRNGContext(0x1339):
        x = np.random.randn(10)
    l = models.Logical('GT', .5, 10)

    res = l(x)
    y = x.copy()
    y[np.greater(x, .5)] = 10
    assert_allclose(res, npres)
    l = models.Logical('lt', compareto, val)
    cond = np.less(x, compareto)
    y = x.copy()
    y[cond] = val[cond]
    assert_allclose(res, npres)
"""

def test_ideal_to_v23_roundtrip():
    """
    Test roundtripping of the transforms.
    """
    v2i = models.V2V3ToIdeal(.4, 450, 730, 1)
    x, y = 450, 730
    assert_allclose(v2i.inverse(*v2i(x, y)), (x, y))


@pytest.mark.parametrize(('wavelength', 'n'),
                         [(1e-6, 1.43079543),
                          (2e-6,  1.42575377),
                          (5e-6, 1.40061966)
                          ])
def test_refraction_index(wavelength, n):
    """
    Tests the computation of the refraction index.
    True values are from the ESA pipeline.
    Reference values are from the PRISM reference file from CV3.
    """
    temp_sys = 37.06107795068881 # in K
    tref = 35 # in K
    pref = 0 # in atm
    pressure_sys = 0 # in atm
    kcoef =  [0.58339748, 0.46085267, 3.8915394]
    lcoef = [0.00252643, 0.010078333, 1200.556]
    tcoef = [-2.66e-05, 0.0, 0.0, 0.0, 0.0, 0.0]
    n_pipeline = models.Snell.compute_refraction_index(wavelength, temp_sys,
                                                       tref, pref, pressure_sys,
                                                       kcoef, lcoef, tcoef)
    assert_allclose(n_pipeline, n)


# slicer
def slicer(x, y):
    import numpy as np
    xc = 0.00039999998989515007
    yc = 0
    xsize = 0.0007999999797903001
    ysize = 0.012000000104308128
    xref = 0
    yref = 0
    rot = 0
    xs = xref + (xc + x * xsize) * np.cos(rot) - (yc + y * ysize)*np.sin(rot)
    ys = yref + (xc + x*xsize)*np.sin(rot) + (yc + y*ysize)*np.cos(rot)
    return xs, ys


def ifupost(lam, x, y):
    pass

import numpy as np
def pcf(x, y):

    xc_in = 0.000355277216
    yc_in = -3.0089012e-05
    xc_out = 0.0487158154447
    yc_out = 0.00856211956976
    xfactor = 0.100989874454
    yfactor = 0.100405184145
    rot = np.deg2rad(-0.129043957046)
    x_forw = {'c0_0': 1117.59960155,
              'c0_1':-239359.848453,
              'c0_2':-907984.325649,
              'c0_3':2253859.06297,
              'c0_4':-288766.031491,
              'c0_5': -10242425.3773,
              'c1_0': -71116.2325003,
              'c1_1': 19985196.6339,
              'c1_2':54765121.9595,
              'c1_3':-92432337.4712,
              'c1_4':16703492.9955,
              'c2_0':1099943.11386,
              'c2_1':-625376936.455,
              'c2_2':-1100537360.67,
              'c2_3':943264210.148,
              'c3_0':15384471.6978,
              'c3_1':8692554365.56,
              'c3_2':7369320811.29,
              'c4_0':-553676297.745,
              'c4_1':-45284426015.5,
              'c5_0':3922639947.11
              }
    x_forw_dist = {'c0_0': -539180414.74,
              'c0_1': 116349945628.0,
              'c0_2':222581438832.0,
              'c0_3':-770393744399.0,
              'c0_4':88796457467.3,
              'c0_5': 729159099666.0,
              'c1_0': 34174541751.1,
              'c1_1': -9.63784782164e+12,
              'c1_2': -1.33095385033e+13,
              'c1_3':3.15895383852e+13,
              'c1_4':-2.50610837913e+12,
              'c2_0':-526120142427.0,
              'c2_1':2.99290404545e+14,
              'c2_2':2.65054129632e+14,
              'c2_3':-3.23571712442e+14,
              'c3_0':-7.36049125988e+12,
              'c3_1':-4.1294633963e+15,
              'c3_2':-1.75789040172e+15,
              'c4_0':2.63537506555e+14,
              'c4_1':2.13599040186e+16,
              'c5_0':-1.86014367616e+15
              }
    y_forw = {'c0_0':-82.3492267824,
              'c0_1':29234.6982762,
              'c0_2':-540260.780853,
              'c0_3':771881.305018,
              'c0_4':-2563462.26848,
              'c0_5': 29914272.1164,
              'c1_0':4513.04082605,
              'c1_1':-2212869.44311,
              'c1_2':32875633.0303,
              'c1_3':-29923698.5288,
              'c1_4':27293902.5636,
              'c2_0':-39820.4434726,
              'c2_1':62431493.9962,
              'c2_2':-667197265.033,
              'c2_3':297253538.182,
              'c3_0':-1838860.86305,
              'c3_1':-777169857.2,
              'c3_2':4514693865.7,
              'c4_0':42790637.764,
              'c4_1':3596423850.94,
              'c5_0':-260274017.448,
              }
    y_forw_dist = {'c0_0':188531839.97,
              'c0_1':-43453434864.0,
              'c0_2':70807756765.8,
              'c0_3':-308272809909.0,
              'c0_4':159768473071.0,
              'c0_5':9.71263334459e+12,
              'c1_0':-11762923852.9,
              'c1_1':3.54593887319e+12,
              'c1_2':-4.19864365542e+12,
              'c1_3':1.25456429831e+13,
              'c1_4':-1.17070515916e+13,
              'c2_0':173091230285.0,
              'c2_1':-1.08534069056e+14,
              'c2_2':8.28933480976e+13,
              'c2_3':-1.24708740989e+14,
              'c3_0':2.77438975799e+12,
              'c3_1':1.4767797203e+15,
              'c3_2':-5.45358301961e+14,
              'c4_0':-9.31015579941e+13,
              'c4_1':-7.53689063943e+15,
              'c5_0':6.46310545048e+14,
              }
    xres = xc_out + xfactor * (np.cos(rot) * (x-xc_in) + np.sin(rot)*(y - yc_in))
    yres = yc_out + yfactor * (-np.sin(rot) * (x-xc_in) + np.cos(rot)*(y-yc_in))
    return xres, yres