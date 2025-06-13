"""Unit tests for matrix_dft module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from jwst.ami import matrix_dft


def test_matrix_dft(circular_pupil):
    # roundtrip
    nlam_d = 100
    npix = circular_pupil.shape[0]
    dft = matrix_dft.matrix_dft(circular_pupil, nlam_d, npix)
    idft = matrix_dft.matrix_idft(dft, nlam_d, npix)

    # test that the IDFT is close to the original pupil
    # can't just use assert_allclose because there are always a few pixels that are not
    # very close given that the input has sharp edges
    isclose = np.isclose(np.abs(idft), circular_pupil, atol=0.1)
    fraction_close = np.sum(isclose) / circular_pupil.size
    assert fraction_close > 0.99

    # test nonzero offset
    offset = (10.9, -7.1)
    dft_offset = matrix_dft.matrix_dft(
        circular_pupil, nlam_d, npix, centering="ADJUSTABLE", offset=offset
    )

    # ensure max value of dft is at the offset from center of array
    argmax_2d = np.unravel_index(np.argmax(np.abs(dft_offset)), dft.shape)
    center = np.array(dft_offset.shape) / 2
    offset_idx = center + np.floor(np.array(offset))
    assert_allclose(argmax_2d, offset_idx)

    # test centering options
    dft_symmetric = matrix_dft.matrix_dft(circular_pupil, nlam_d, npix, centering="SYMMETRIC")
    dft_adjustable = matrix_dft.matrix_dft(circular_pupil, nlam_d, npix, centering="ADJUSTABLE")
    assert_allclose(np.abs(dft_symmetric), np.abs(dft_adjustable))

    # test expected failure with offset non-tuple
    with pytest.raises(ValueError):
        dft = matrix_dft.matrix_dft(circular_pupil, nlam_d, npix, centering="ADJUSTABLE", offset=1)

    # test nlam_d and npix as tuples
    dft = matrix_dft.matrix_dft(circular_pupil, (nlam_d, nlam_d), (npix, npix))
    assert dft.shape == (npix, npix)

    # test expected failure with nlam_d and npix as non-tuple
    npix_bad = (npix, npix, npix)
    nlam_d_bad = (nlam_d, nlam_d, nlam_d)
    with pytest.raises(ValueError):
        dft = matrix_dft.matrix_dft(circular_pupil, nlam_d_bad, npix)
    with pytest.raises(ValueError):
        dft = matrix_dft.matrix_dft(circular_pupil, nlam_d, npix_bad)

    # test expected failure invalid centering style
    with pytest.raises(ValueError):
        dft = matrix_dft.matrix_dft(circular_pupil, nlam_d, npix, centering="INVALID")
