import numpy as np

from jwst.refpix.irs2_subtract_reference import decode_mask, clobber_ref


def test_clobber_ref():
    data = np.ones((2, 3, 5, 3200))
    data[:, :, :] = np.arange(3200, dtype=float)

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
    is_irs2 = np.full(3200, True)

    # mark a couple pixels already bad
    ref_flags = np.full(3200, False)
    ref_flags[1908:1910] = True
    ref_flags[1912:1914] = True

    clobber_ref(data, output, odd_even, mask, ref_flags, is_irs2)

    compare = np.ones((2, 3, 5, 3200))
    compare[:, :, :] = np.arange(3200, dtype=float)

    # next pixel is bad, no lower value, neighbor is okay:
    # replace with neighbor
    compare[..., 648: 648+2] = [646, 647]
    # lower pixel is bad, replace with upper pixel
    compare[..., 668: 668+2] = [688, 689]
    # upper pixel is bad, replace with lower pixel
    compare[..., 690: 690+2] = [670, 671]
    # lower pixel is bad, replace with upper pixel
    compare[..., 710: 710+2] = [730, 731]
    # upper pixel is bad, replace with lower pixel
    compare[..., 1890: 1890+2] = [1870, 1871]
    # lower bad, no upper, no good neighbors, replace with 0.0
    compare[..., 1910: 1910+2] = 0.0
    # upper is bad, replace with lower
    compare[..., 1808: 1808+2] = [1788, 1789]
    # lower is bad, replace with upper
    compare[..., 1828: 1828+2] = [1848, 1849]
    # both good, replace with average
    compare[..., 2028: 2028+2] = [2028, 2029]
    compare[..., 2068: 2068+2] = [2068, 2069]
    compare[..., 2150: 2150+2] = [2150, 2151]
    compare[..., 2190: 2190+2] = [2190, 2191]

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
