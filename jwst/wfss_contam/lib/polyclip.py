import os
from glob import glob
import numpy.ctypeslib as npct
import numpy as np
import ctypes
from ctypes import c_int

# from . import polyclip_c
this_path = os.path.split(__file__)[0]
try:
    so_file = glob(os.path.join(this_path, 'polyclip_c*.so'))[0]
    so_file = so_file[0]
except IndexError:
    print("WARNING: Cannot find polyclip_c*.so library")
    so_file = ''
polyclip = ctypes.cdll.LoadLibrary(so_file)

array_1d_int_l = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
array_1d_int_r = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
array_1d_int_b = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
array_1d_int_t = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')

array_1d_double_px = npct.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')
array_1d_double_py = npct.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')
array_1d_double_px_out = npct.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')
array_1d_double_py_out = npct.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')
array_1d_double_ri_out = npct.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')
array_1d_double_areas = npct.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')
array_1d_double_nclip_poly = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
array_1d_int_poly_inds = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
array_1d_double_inds = npct.ndpointer(dtype=np.int32, ndim=2, flags='CONTIGUOUS')
array_1d_double_x = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
array_1d_double_y = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
array_1d_double_index = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
polyclip.polyclip_multi2.restype = None

polyclip.polyclip_multi4.argtypes = [array_1d_int_l,              # l
                                     array_1d_int_r,              # r
                                     array_1d_int_b,              # b
                                     array_1d_int_t,              # t
                                     array_1d_double_px,          # px
                                     array_1d_double_py,          # py
                                     c_int,                       # n_poly
                                     array_1d_int_poly_inds,      # poly_inds
                                     array_1d_double_x,           # x
                                     array_1d_double_y,           # y
                                     array_1d_double_nclip_poly,  # nclip_poly
                                     array_1d_double_areas,       # areas
                                     array_1d_double_index        # output index
                                     ]
