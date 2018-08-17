import numpy as np

class WCS:
    """Temporary interface for WCS using a 1-D array.

    nelem = wcs.get_nelem()
    dtype = wcs.get_dtype()
    (wl0, wl1) = wcs.get_range()
    pixel = wcs.to_pixel(wavelength)
    wavelength = wcs.to_world(pixel)
    wcs.close()

    Parameters
    ----------
    input: 1-D ndarray
        This is a 1-D array of world coordinate values, e.g.  wavelengths.
        The corresponding pixel coordinates are taken to be the array
        index numbers (zero indexed).
    """

    def __init__(self, input):

        self.world = None
        self.dtype = None
        self.nelem = 0
        self.monotonic = True

        self.world = np.array(input)
        if len(self.world.shape) > 1:
            raise RuntimeError("temp_wcs expected a 1-D array")

        self.dtype = self.world.dtype
        self.nelem = len(self.world)

        if self.nelem > 2:
            dw = self.world[1:] - self.world[0:-1]
            if dw.min() * dw.max() < 0.:
                print("Warning:  the specified array is not monotonic.")
                self.monotonic = False
            else:
                self.monotonic = True

    def get_nelem(self):
        """Return the number of elements in the world coordinate array."""
        return self.nelem

    def get_dtype(self):
        """Return the data type of the world coordinate array."""
        return self.dtype

    def get_range(self):
        """Return the minimum and maximum world coordinate values."""
        w0 = self.world[0]
        w1 = self.world[-1]
        if w1 < w0:
            temp = w0
            w0 = w1
            w1 = temp
        return (w0, w1)

    def to_pixel(self, w):

        if not self.nelem or self.nelem < 1:
            return None

        if not self.monotonic:
            raise RuntimeError("World coordinates are not monotonic.")

        try:
            len_w = len(w)
            w_has_len = True
        except TypeError:
            w_has_len = False

        if self.nelem == 1:
            if w_has_len:
                return np.zeros(len_w, dtype=self.dtype)
            else:
                return 0.

        if self.world[-1] > self.world[0]:
            dir = 1.
        else:
            dir = -1.

        if w_has_len:
            w = np.array(w)
        else:
            w = np.zeros(1, dtype=np.float64) + w

        pixel = w.copy()
        if dir > 0.:
            # Yes, too much looping; this is a stub!
            for i in range(len(w)):
                if w[i] < self.world[0]:
                    # Linear extrapolation.
                    pixel[i] = (w[i] - self.world[0]) / \
                               (self.world[1] - self.world[0])
                elif w[i] > self.world[-1]:
                    # Linear extrapolation.
                    pixel[i] = float(self.nelem - 1) + \
                               (w[i] - self.world[-1]) / \
                               (self.world[-1] - self.world[-2])
                else:
                    for k in range(self.nelem - 1):
                        # Linear interpolation.
                        if w[i] >= self.world[k] and \
                           w[i] <= self.world[k + 1]:
                            pixel[i] = float(k) + (w[i] - self.world[k]) / \
                                           (self.world[k + 1] - self.world[k])
                            break
        else:                                   # dir < 0
            for i in range(len(w)):
                if w[i] > self.world[0]:
                    pixel[i] = (w[i] - self.world[0]) / \
                               (self.world[1] - self.world[0])
                elif w[i] < self.world[-1]:
                    pixel[i] = float(self.nelem - 1) + \
                               (w[i] - self.world[-1]) / \
                               (self.world[-1] - self.world[-2])
                else:
                    for k in range(self.nelem - 1):
                        if w[i] <= self.world[k] and \
                           w[i] >= self.world[k + 1]:
                            pixel[i] = float(k + 1) - \
                                       (w[i] - self.world[k + 1]) / \
                                       (self.world[k] - self.world[k + 1])
                            break

        if w_has_len:
            return pixel
        else:
            return float(pixel[0])

    def to_world(self, pixel):

        if not self.nelem or self.nelem < 1:
            return None

        try:
            len_pixel = len(pixel)
            pixel_has_len = True
        except TypeError:
            pixel_has_len = False

        if self.nelem == 1:
            if pixel_has_len:
                return np.zeros(len_pixel, dtype=self.dtype) + self.world[0]
            else:
                return self.world[0]

        if pixel_has_len:
            pixel = np.array(pixel)
        else:
            pixel = np.zeros(1, dtype=np.float64) + pixel

        ix = np.floor(pixel).astype(np.intp)
        ix = np.where(ix < 0, 0, ix)
        ix = np.where(ix >= self.nelem - 2, self.nelem - 2, ix)
        v1 = self.world[ix]
        v2 = self.world[ix + 1]

        p = pixel - ix
        q = 1. - p
        value = q * v1 + p * v2

        if pixel_has_len:
            return value
        else:
            return float(value[0])

    def close(self):
        if self.world is not None:
            del self.world
            self.world = None
        self.dtype = None
        self.nelem = 0
        self.monotonic = True
