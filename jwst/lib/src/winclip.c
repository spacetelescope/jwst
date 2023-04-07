/*
This module implements window clipping algorithm used for spectral extraction.
This algorithm clips (lists of) rectangles (with sides parallel to axes)
to (lists of) square windows and also one large "image rectangle".

Main function for Python: get_clipped_pixels.

Python signature: get_clipped_pixels(x, y, padding, nx, ny, w, h)

Vertices of the rectangles to be clipped are formed by adding or subtracting
w/2 to(from) x and h/2 to(from) y coordinates (so x, y are at the center of
the rectangles to be clipped). The clipping window is an intersection of two
windows: a large one representing image boundary [0, nx] and [0, ny] and a
second window formed by adding(subtracting) 'padding' to(from) coordinates
x and y converted to integer.

The output of this function is a tuple of four arrays: integer
coordinates (i, j) of clipped rectangle(s) pixels, index of the polygon
to which the pixel belongs (index in input coordinate arrays), and area of the
pixel inside the clipping window.

Parameters
----------
x : numpy.ndarray
    X-coordinate of the center rectangle to be clipped. Center of the clipping
    window is computed as int(x) and the clipping window is formed by adding or
    subtracting 'padding' from these integer coordinate values.

y : numpy.ndarray
    Y-coordinates. Same explanation as for X-coordinate above.

padding : int
    Number of pixels to add around the center pixel of the clipping window
    (int(x), int(y)).

nx : int
ny : int
    Width and height of the "image" in pixels - the large clipping window.

w : double
h : double
    Width and height of the rectangle to be clipped.

Returns
-------

x_arr : numpy.ndarray
y_arr : numpy.ndarray
    Integer coordinates of the pixels within clipped rectangles that have
    non-zeros area.

areas_arr : numpy.ndarray
    Array containing area of each pixel in x_arr and y_arr.

idx_arr : numpy.ndarray
    Indices of rectangles from input arrays. For example, input coordinate
    arrays 'x' has N elements. Each element represents the coordinate of one
    rectangle to be clipped (or a clipping window). This 'idx_arr' contains
    the index of the clipping rectangle to which each *clipped* pixel
    (in the x_arr and y_arr arrays) belongs.

*/

#include <stdlib.h>
#include <math.h>

#include <Python.h>

#define PY_ARRAY_UNIQUE_SYMBOL _jwst_winclip_numpy_api
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

typedef double dbl_type;
static const int npy_dbl = NPY_DOUBLE;
static const double full_pixel_area = 1.0;

// if float32 accuracy is enough, use the following types instead:
/*
typedef float dbl_type;
static const int npy_dbl = NPY_FLOAT;
static const float full_pixel_area = 1.0f;
*/


inline dbl_type dbl_min(dbl_type x, dbl_type y) {
    return (x < y) ? x : y;
}


int mem_alloc(int nelem, int **xv, int **yv, dbl_type **av, int **idx) {
    int *x, *y, *i;
    dbl_type *a;
    const char *msg = "Couldn't allocate memory for output arrays.";

    // x-array:
    x = (int*) realloc(*xv, nelem * sizeof(int));
    if (x) {
        *xv = x;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    // y-array:
    y = (int*) realloc(*yv, nelem * sizeof(int));
    if (y) {
        *yv = y;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    // area-array:
    a = (dbl_type*) realloc(*av, nelem * sizeof(dbl_type));
    if (a) {
        *av = a;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    // idx-array:
    i = (int*) realloc(*idx, nelem * sizeof(int));
    if (i) {
        *idx = i;
    } else {
        PyErr_SetString(PyExc_MemoryError, msg);
        return 1;
    }

    return 0;
}

/*
Caller is responsible for de-allocating memory for the
return values: x, y, areas, idx
*/
int clip_pixels(int n, dbl_type *xc, dbl_type *yc, int padding,
    int nx, int ny, double w, double h,
    int *npix, int **x, int **y, dbl_type **areas, int **idx) {
    double x1, x2, y1, y2;
    double wx1, wx2, wy1, wy2;  // temp values for fractional pixel size
    int npixe;  // estimate of number of output pixels per input coord
    int nalloc;  // number of allocated output values (not bytes)
    int l, r, b, t, ix, iy, np, chunk, tnpix = 0;
    int imin, imax, jmin, jmax, imin1, jmin1;
    dbl_type *av = NULL;  // vector for clipped pixel area
    int *xv = NULL, *yv = NULL, *iv = NULL;  // vectors of coordinates
    int i, j, k;

    npixe = (int) (dbl_min(ceil(w + 1.0), 2.0 * padding + 1.0) *
                   dbl_min(ceil(h + 1.0), 2.0 * padding + 1.0));

    // allocate enough memory to hold output pixel data (coordinates and area)
    // for about 1/10 of input coordinates (for large inputs):
    chunk =  npixe * ((n / 10) ? (n / 10) : 10);
    nalloc = n * npixe;
    if (mem_alloc(nalloc, &xv, &yv, &av, &iv)) return 1;

    nx -= 1;
    ny -= 1;

    for (k = 0; k < n; k++) {
        // compute clipping window:
        ix = (int) xc[k];
        iy = (int) yc[k];

        l = ix - padding;
        r = ix + padding;
        b = iy - padding;
        t = iy + padding;

        // clip l, r, b, t to 0...nx(ny) - image window
        if ((l > nx) || (r < 0) || (b > ny) || (t < 0)) continue;
        if (l < 0) l = 0;
        if (r > nx) r = nx;
        if (b < 0) b = 0;
        if (t > ny) t = ny;

        // compute rectangle edge positions:
        x1 = ((double) xc[k]) - w / 2.0;
        x2 = x1 + w;
        y1 = ((double) yc[k]) - h / 2.0;
        y2 = y1 + h;

        // For testing purposes, to simulate accuracy loss in wfss_contam algorithm,
        // uncomment code below:
        /*
        x1 = (float) (xc[k] - w / 2.0);
        x2 = (float) (xc[k] + w / 2.0);
        y1 = (float) (yc[k] - h / 2.0);
        y2 = (float) (yc[k] + h / 2.0);
        */

        // compute indices' ranges for the rectangle clipped to the window:
        if ((imin = (int) floor(x1)) < l) {
            imin = l;
            x1 = (double) imin;
        }
        if ((imax = (int) ceil(x2 - 1.0)) > r) {
            imax = r;
            x2 = (double) (imax + 1);
        }
        if (x1 >= x2) continue;
        if ((jmin = (int) floor(y1)) < b) {
            jmin = b;
            y1 = (double) jmin;
        }
        if ((jmax = (int) ceil(y2 - 1)) > t) {
            jmax = t;
            y2 = (double) (jmax + 1);
        }
        if (y1 >= y2) continue;
        imin1 = imin + 1;
        jmin1 = jmin + 1;

        np = (jmax - jmin + 1) * (imax - imin + 1);

        if (tnpix + np > nalloc) {
            // allocate more memory for output vectors:
            nalloc += chunk;
            if (mem_alloc(nalloc, &xv, &yv, &av, &iv)) return 1;
        }

        // pre-compute fractional pixel sizes:
        wx1 = (x2 > imin1) ? (imin1 - x1) : (x2 - x1);
        wx2 = x2 - imax;
        wy1 = (y2 > jmin1) ? (jmin1 - y1) : (y2 - y1);
        wy2 = y2 - jmax;

        // compute pixel areas

        // 1. process left column:
        //   - 1.a: bottom (small j)-left pixel
        xv[tnpix] = imin;
        yv[tnpix] = jmin;
        av[tnpix] = (dbl_type) (wx1 * wy1);
        iv[tnpix] = k;
        tnpix++;

        if (jmax > jmin) {
            //   -1.b: top (large j)-left corner pixel:
            xv[tnpix] = imin;
            yv[tnpix] = jmax;
            av[tnpix] = (dbl_type) (wx1 * wy2);
            iv[tnpix] = k;
            tnpix++;

            //   - 1.c: middle of the left column:
            for (j = jmin + 1; j < jmax; j++) {
                xv[tnpix] = imin;
                yv[tnpix] = j;
                av[tnpix] = (dbl_type) wx1;
                iv[tnpix] = k;
                tnpix++;
            }
        }

        if (imax > imin) {
            // 2. process middle columns:
            for (i = imin + 1; i < imax; i++) {
                //   -2.a: bottom pixel:
                xv[tnpix] = i;
                yv[tnpix] = jmin;
                av[tnpix] = (dbl_type) wy1;
                iv[tnpix] = k;
                tnpix++;

                if (jmax > jmin) {
                    //   -2.b: middle of each column:
                    for (j = jmin + 1; j < jmax; j++) {
                        xv[tnpix] = i;
                        yv[tnpix] = j;
                        av[tnpix] = full_pixel_area;
                        iv[tnpix] = k;
                        tnpix++;
                    }
                    //   -2.c: top pixel:
                    xv[tnpix] = i;
                    yv[tnpix] = jmax;
                    av[tnpix] = (dbl_type) wy2;
                    iv[tnpix] = k;
                    tnpix++;
                }
            }

            // 3. process right column (if distinct from left):
            //   -3.a: top-right corner pixel:
            xv[tnpix] = imax;
            yv[tnpix] = jmin;
            av[tnpix] = (dbl_type) (wx2 * wy1);
            iv[tnpix] = k;
            tnpix++;

            if (jmax > jmin) {
                //   -3.b: middle of the right column:
                for (j = jmin + 1; j < jmax; j++) {
                    xv[tnpix] = imax;
                    yv[tnpix] = j;
                    av[tnpix] = (dbl_type) wx2;
                    iv[tnpix] = k;
                    tnpix++;
                }
                //   -3.c: bottom-right corner pixel:
                xv[tnpix] = imax;
                yv[tnpix] = jmax;
                av[tnpix] = (dbl_type) (wx2 * wy2);
                iv[tnpix] = k;
                tnpix++;
            }
        }
    }

    // trim memory arrays:
    if ((tnpix < nalloc) && mem_alloc(tnpix, &xv, &yv, &av, &iv)) return 1;

    // assign output values:
    *npix = tnpix;
    *x = xv;
    *y = yv;
    *areas = av;
    *idx = iv;

    return 0;
}


PyArrayObject * ensure_array(PyObject *obj, int *is_copy) {
    if (PyArray_CheckExact(obj) &&
        PyArray_IS_C_CONTIGUOUS((PyArrayObject *) obj) &&
        PyArray_TYPE((PyArrayObject *) obj) == npy_dbl) {
        *is_copy = 0;
        return (PyArrayObject *) obj;
    } else {
        *is_copy = 1;
        return (PyArrayObject *) PyArray_FromAny(
            obj, PyArray_DescrFromType(npy_dbl), 0, 0,
            NPY_ARRAY_CARRAY | NPY_ARRAY_FORCECAST, NULL
        );
    }
}


static PyObject * get_clipped_pixels(PyObject *module, PyObject *args) {
    PyObject *result = NULL, *xco, *yco;
    double w, h;
    dbl_type *areas=NULL;
    int n, nx, ny, npix, padding;
    int *x=NULL, *y=NULL, *idx=NULL;
    int free_xc=0, free_yc=0, status=0;
    PyArrayObject *xc, *yc;
    PyArrayObject *x_arr=NULL, *y_arr=NULL, *areas_arr=NULL, *idx_arr=NULL;
    npy_intp npy_npix = 0;

    if (!PyArg_ParseTuple(args, "OOiiidd:get_clipped_pixels",
        &xco, &yco, &padding, &nx, &ny, &w, &h)) {
    return NULL;
    }

    // check that input parameters are valid:
    if (padding < 0) {
        PyErr_SetString(PyExc_ValueError,
            "'padding' must be a strictly positive number.");
        return NULL;
    }

    if ((nx < 1) || (ny < 1)) {
        PyErr_SetString(PyExc_ValueError,
            "'nx' and 'ny' must be a strictly positive integer numbers.");
        return NULL;
    }

    if ((w <= 0) || (h <= 0)) {
        PyErr_SetString(PyExc_ValueError,
            "'w' and 'h' must be a strictly positive numbers.");
        return NULL;
    }

    // ensure we are working with numpy arrays and avoid creating new ones
    // if possible:
    if ((!(xc = ensure_array(xco, &free_xc))) ||
        (!(yc = ensure_array(yco, &free_yc)))) {
        goto cleanup;
    }

    n = (int) PyArray_Size((PyObject *) xc);
    if (n != PyArray_Size((PyObject *) yc)) {
        PyErr_SetString(PyExc_ValueError,
            "Input coordinate arrays of unequal size.");
        goto cleanup;
    }
    if (!n) {
        // 0-length input arrays. Nothing to clip. Return 0-length arrays
        x_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_npix, NPY_INT, 0);
        if (!x_arr) goto fail;

        y_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_npix, NPY_INT, 0);
        if (!y_arr) goto fail;

        areas_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_npix, npy_dbl, 0);
        if (!areas_arr) goto fail;

        idx_arr = (PyArrayObject*) PyArray_EMPTY(1, &npy_npix, NPY_INT, 0);
        if (!idx_arr) goto fail;

        result = Py_BuildValue("(NNNN)", x_arr, y_arr, areas_arr, idx_arr);
        goto cleanup;
    }

    status = clip_pixels(n, (dbl_type *) PyArray_DATA(xc),
        (dbl_type *) PyArray_DATA(yc), padding, nx, ny, w, h,
        &npix, &x, &y, &areas, &idx);

    if (status) {
        goto fail;

    } else {
        // create return tuple:
        npy_npix = (npy_intp) npix;
        x_arr = (PyArrayObject*) PyArray_SimpleNewFromData(
            1, &npy_npix, NPY_INT, x
        );
        if (!x_arr) goto fail;
        x = NULL;

        y_arr = (PyArrayObject*) PyArray_SimpleNewFromData(
            1, &npy_npix, NPY_INT, y
        );
        if (!y_arr) goto fail;
        y = NULL;

        areas_arr = (PyArrayObject*) PyArray_SimpleNewFromData(
            1, &npy_npix, npy_dbl, areas
        );
        if (!areas_arr) goto fail;
        areas = NULL;

        idx_arr = (PyArrayObject*) PyArray_SimpleNewFromData(
            1, &npy_npix, NPY_INT, idx
        );
        if (!idx_arr) goto fail;
        idx = NULL;

        PyArray_ENABLEFLAGS(x_arr, NPY_ARRAY_OWNDATA);
        PyArray_ENABLEFLAGS(y_arr, NPY_ARRAY_OWNDATA);
        PyArray_ENABLEFLAGS(areas_arr, NPY_ARRAY_OWNDATA);
        PyArray_ENABLEFLAGS(idx_arr, NPY_ARRAY_OWNDATA);

        result = Py_BuildValue("(NNNN)", x_arr, y_arr, areas_arr, idx_arr);
        goto cleanup;
    }

fail:
    Py_XDECREF(x_arr);
    Py_XDECREF(y_arr);
    Py_XDECREF(areas_arr);
    Py_XDECREF(idx_arr);
    free(x);
    free(y);
    free(areas);
    free(idx);

    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_MemoryError,
            "Unable to allocate memory for output arrays.");
    }

cleanup:
    if (free_xc) Py_XDECREF(xc);
    if (free_yc) Py_XDECREF(yc);

    return result;
}


static PyMethodDef winclip_methods[] =
{
    {
        "get_clipped_pixels",  get_clipped_pixels, METH_VARARGS,
        "get_clipped_pixels(image, histogram, minValue, maxValue, binWidth)"
    },
    {0, 0}  /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "winclip",                   /* m_name */
    "Clip a rectangle to multiple rectangular windows",  /* m_doc */
    -1,                          /* m_size */
    winclip_methods,             /* m_methods */
    NULL,                        /* m_reload */
    NULL,                        /* m_traverse */
    NULL,                        /* m_clear */
    NULL,                        /* m_free */
};

PyMODINIT_FUNC PyInit_winclip(void)
{
    PyObject* m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
