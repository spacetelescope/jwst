import numpy as np

from jwst.refpix.irs2_subtract_reference import decode_mask, clobber_ref


def test_clobber_ref():
    data = np.ones((2, 3, 5, 3200), dtype=np.float32)

    output = np.array([1, 1, 2, 2, 3, 3, 4, 4], dtype=np.int16)
    odd_even = np.array([1, 2, 1, 2, 1, 2, 1, 2], dtype=np.int16)
    mask = np.array([1 + 2**1,
                     2**2 + 2**3,
                     2**30 + 2**31,
                     2**26 + 2**27,
                     2**5 + 2**7,
                     2**11 + 2**13,
                     0,
                     2**4],
                    dtype=np.uint32)

    clobber_ref(data, output, odd_even, mask)

    compare = np.ones((2, 3, 5, 3200), dtype=np.float32)
    compare[:, :, :, 648] = 0.
    compare[:, :, :, 649] = 0.
    compare[:, :, :, 668] = 0.
    compare[:, :, :, 669] = 0.
    compare[:, :, :, 690] = 0.
    compare[:, :, :, 691] = 0.
    compare[:, :, :, 710] = 0.
    compare[:, :, :, 711] = 0.
    compare[:, :, :, 1290] = 0.
    compare[:, :, :, 1291] = 0.
    compare[:, :, :, 1310] = 0.
    compare[:, :, :, 1311] = 0.
    compare[:, :, :, 1368] = 0.
    compare[:, :, :, 1369] = 0.
    compare[:, :, :, 1388] = 0.
    compare[:, :, :, 1389] = 0.
    compare[:, :, :, 2028] = 0.
    compare[:, :, :, 2029] = 0.
    compare[:, :, :, 2068] = 0.
    compare[:, :, :, 2069] = 0.
    compare[:, :, :, 2150] = 0.
    compare[:, :, :, 2151] = 0.
    compare[:, :, :, 2190] = 0.
    compare[:, :, :, 2191] = 0.
    compare[:, :, :, 3108] = 0.
    compare[:, :, :, 3109] = 0.

    assert np.allclose(data, compare)


def test_decode_mask():

    output = np.array([1, 1, 2, 2, 3, 3, 4, 4], dtype=np.int16)

    mask = np.array([1048608, 0, 8464, 8, 16448, 33554944, 9, 32897],
                    dtype=np.uint32)

    nrows = len(output)
    check = np.zeros(nrows, dtype=bool)
    compare = [[5, 20],
               [],
               [18, 23, 27],
               [28],
               [6, 14],
               [9, 25],
               [28, 31],
               [16, 24, 31]]
    for row in range(nrows):
        bits = decode_mask(output[row], mask[row])
        check[row] = (bits == compare[row])

    assert np.all(check)
