"""
Test for flat_field.interpolate_flat
"""

import numpy as np

from jwst.flatfield.flat_field import interpolate_flat

nz = 6
ny = 5
nx = 7

image_wl = np.arange(nz, dtype=np.float32) + 4.5
wl = np.arange(nx * ny, dtype=np.float32).reshape(ny, nx) - 5.0


def test_interpolate_flat_1():

    image_flat = np.arange(nx * ny, dtype=np.float32).reshape(ny, nx) - 5.0
    image_err = np.arange(nx * ny, dtype=np.float32).reshape(ny, nx) - 7.0
    image_dq = np.arange(nx * ny, dtype=np.int32).reshape(ny, nx)

    # Since image_flat is 2-D, the inputs will be returned unchanged.
    output = interpolate_flat(image_flat, image_dq, image_err, image_wl, wl)
    assert np.allclose(output[0], image_flat, atol=1.0e-6)
    assert np.allclose(output[1], image_dq, atol=0)
    assert np.allclose(output[2], image_err, atol=1.0e-6)


def test_interpolate_flat_2():

    # 2-D, but reshaped to 3-D
    image_flat = np.arange(nx * ny, dtype=np.float32).reshape(1, ny, nx) - 5.0
    image_err = np.arange(nx * ny, dtype=np.float32).reshape(1, ny, nx) - 7.0
    # 2-D
    image_dq = np.arange(nx * ny, dtype=np.int32).reshape(ny, nx)
    output = interpolate_flat(image_flat, image_dq, image_err, image_wl, wl)
    assert np.allclose(output[0], image_flat[0, :, :], atol=1.0e-6)
    assert np.allclose(output[1], image_dq, atol=0)
    assert np.allclose(output[2], image_err[0, :, :], atol=1.0e-6)


def test_interpolate_flat_3():

    # 2-D, but reshaped to 3-D
    image_flat = np.arange(nx * ny, dtype=np.float32).reshape(1, ny, nx) - 5.0
    image_err = np.arange(nx * ny, dtype=np.float32).reshape(1, ny, nx) - 7.0
    # 2-D, but reshaped to 3-D
    image_dq = np.arange(nx * ny, dtype=np.int32).reshape(1, ny, nx)
    output = interpolate_flat(image_flat, image_dq, image_err, image_wl, wl)
    assert np.allclose(output[0], image_flat[0, :, :], atol=1.0e-6)
    assert np.allclose(output[1], image_dq[0, :, :], atol=0)
    assert np.allclose(output[2], image_err[0, :, :], atol=1.0e-6)


def test_interpolate_flat_4():

    # 3-D
    image_flat = np.arange(nx * ny * nz, dtype=np.float32).reshape(nz, ny, nx) - 5.0
    image_err = np.arange(nx * ny * nz, dtype=np.float32).reshape(nz, ny, nx) - 7.0
    image_dq = np.zeros((nz, ny, nx), dtype=np.uint32)
    image_dq[:, 1, 3] = 2
    output = interpolate_flat(image_flat, image_dq, image_err, image_wl, wl)

    expected_value_0 = np.array(
        [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 22.5, 58.5, 94.5, 130.5],
            [166.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ],
        dtype=np.float32,
    )

    expected_value_1 = np.array(
        [
            [1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 2, 0, 0, 0],
            [0, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1],
        ],
        dtype=np.uint32,
    )

    expected_value_2 = np.array(
        [
            [-339.5, -303.5, -267.5, -231.5, -195.5, -159.5, -123.5],
            [-87.5, -51.5, -15.5, 20.5, 56.5, 92.5, 128.5],
            [164.5, 200.5, 236.5, 272.5, 308.5, 344.5, 380.5],
            [416.5, 452.5, 488.5, 524.5, 560.5, 596.5, 632.5],
            [668.5, 704.5, 740.5, 776.5, 812.5, 848.5, 884.5],
        ],
        dtype=np.float32,
    )

    assert np.allclose(output[0], expected_value_0, atol=1.0e-6)
    assert np.allclose(output[1], expected_value_1, atol=0)
    assert np.allclose(output[2], expected_value_2, atol=1.0e-6)
