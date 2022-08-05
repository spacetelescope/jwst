"""
Unit tests for straylight correction
"""

from jwst.straylight.straylight import makemodel_ccode, makemodel_composite
import numpy as np


def test_correct_mrs_xartifact():
    """ Test Correct Straylight routine gives expected results for small region """

    image = np.zeros((1024, 1032))
    image[:, 300] = 100
    istart, istop = 0, 516
    xvec = np.arange(1032)
    lorfwhm = np.zeros(1024) + 100
    lorscale = np.zeros(1024) + 0.001
    gauxoff = np.zeros(1024) + 10
    gaufwhm = np.zeros(1024) + 7
    gauscale1 = np.zeros(1024) + 0.001
    gauscale2 = np.zeros(1024) + 0.001

    # Test C version
    result_c = makemodel_ccode(image, xvec, istart, istop, lorfwhm, lorscale,
                               gaufwhm, gauxoff, gauscale1, gauscale2)
    cutout_c = result_c[500, 270:330]

    # Test python version
    result_py = makemodel_composite(image, xvec, istart, istop, lorfwhm, lorscale,
                                    gaufwhm, gauxoff, gauscale1, gauscale2)
    cutout_py = result_py[500, 270:330]

    compare = np.array([0.09783203, 0.10662437, 0.11656731, 0.12742336, 0.13880997,
                        0.15021261, 0.1610211, 0.17058835, 0.17830785, 0.18370678,
                        0.18655573, 0.18699956, 0.18570036, 0.18393369, 0.18349825,
                        0.18625805, 0.19326517, 0.20376261, 0.21473958, 0.2216788,
                        0.22045676, 0.2094179, 0.1902925, 0.16733257, 0.14527816,
                        0.12747414, 0.11511109, 0.10763159, 0.10367233, 0.10188942,
                        0.10139532, 0.10188942, 0.10367233, 0.10763159, 0.11511109,
                        0.12747414, 0.14527816, 0.16733257, 0.1902925, 0.2094179,
                        0.22045676, 0.2216788, 0.21473958, 0.20376261, 0.19326517,
                        0.18625805, 0.18349825, 0.18393369, 0.18570036, 0.18699956,
                        0.18655573, 0.18370678, 0.17830785, 0.17058835, 0.1610211,
                        0.15021261, 0.13880997, 0.12742336, 0.11656731, 0.10662437])

    assert np.allclose(compare, cutout_c, rtol=1e-6)
    assert np.allclose(compare, cutout_py, rtol=1e-6)
