from jwst.lib.winclip import get_clipped_pixels
import numpy as np


def test_clip_to_large_window():
    xi, yi, ai, ki = get_clipped_pixels(
        [1.5, 2, 2.5],
        [3.5, 4, 4],
        1, 10, 10, 1, 1
    )

    assert np.all(xi == [1, 1, 1, 2, 2, 2, 2])
    assert np.all(yi == [3, 3, 4, 3, 4, 3, 4])
    assert np.all(ki == [0, 1, 1, 1, 1, 2, 2])
    assert np.allclose(ai, [1., 0.25, 0.25, 0.25, 0.25, 0.5, 0.5])


def test_clip_large_rect_to_large_window():
    xi, yi, ai, ki = get_clipped_pixels(
        [0.1, 2, 2.5],
        [3.9, 4, 4],
        1, 10, 10, 1.2, 1.8
    )

    assert np.all(xi == [0, 0, 1, 1, 2, 2, 1, 1, 2, 2, 3, 3])
    assert np.all(yi == [3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4])
    assert np.all(ki == [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2])
    assert np.allclose(
        ai,
        [0.7, 0.56, 0.54, 0.54, 0.54, 0.54,
         0.09, 0.09, 0.90, 0.90, 0.09, 0.09]
    )


def test_clip_small_rect_at_edges():
    xi, yi, ai, ki = get_clipped_pixels(
        [1.5, 1.9, -0.2, -0.15, -0.15, -0.15],
        [3.5, 0.9, 0.9, 0.9, -2, -0.15],
        0, 10, 10, 0.4, 0.4
    )
    assert np.all(xi == [1, 1, 0, 0])
    assert np.all(yi == [3, 0, 0, 0])
    assert np.all(ki == [0, 1, 3, 5])
    assert np.allclose(ai, [0.16, 0.09, 0.015, 0.0025])
