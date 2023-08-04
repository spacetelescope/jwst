import numpy as np

from jwst.refpix.irs2_subtract_reference import decode_mask, clobber_ref


def test_clobber_ref():
    data = np.ones((2, 3, 5, 3200))

    output = np.array([1, 1, 2, 2, 3, 3, 4, 4], dtype=np.int16)
    odd_even = np.array([1, 2, 1, 2, 1, 2, 1, 2], dtype=np.int16)
    mask = np.array([1 + 2**1,
                     2**2 + 2**3,
                     2**30 + 2**31,
                     2**26 + 2**27,
                     2**5 + 2**7,
                     2**11 + 2**13,
                     0,
                     0],
                    dtype=np.uint32)

    clobber_ref(data, output, odd_even, mask)

    compare = np.ones((2, 3, 5, 3200))
    compare[..., 648: 648+4] = 0.
    compare[..., 668: 668+4] = 0.
    compare[..., 690: 690+4] = 0.
    compare[..., 710: 710+4] = 0.
    compare[..., 1890: 1890+4] = 0.
    compare[..., 1910: 1910+4] = 0.
    compare[..., 1808: 1808+4] = 0.
    compare[..., 1828: 1828+4] = 0.
    compare[..., 2028: 2028+4] = 0.
    compare[..., 2068: 2068+4] = 0.
    compare[..., 2150: 2150+4] = 0.
    compare[..., 2190: 2190+4] = 0.

    assert np.allclose(data, compare)


def test_decode_mask():

    output = np.array([1, 1, 2, 2, 3, 3, 4, 4], dtype=np.int16)

    mask = np.array([1048608, 0, 8464, 8, 16448, 33554944, 9, 32897],
                    dtype=np.uint32)

    nrows = len(output)
    check = np.zeros(nrows, dtype=bool)
    compare = [[5, 20],
               [],
               [4, 8, 13],
               [3],
               [6, 14],
               [9, 25],
               [0, 3],
               [0, 7, 15]]
    for row in range(nrows):
        bits = decode_mask(output[row], mask[row])
        check[row] = (bits == compare[row])

    assert np.all(check)
