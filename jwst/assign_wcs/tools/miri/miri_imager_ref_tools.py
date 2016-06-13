# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools to create reference files for MIRI Imager using IDT reference
products delivered with CDP-3:

MIRI_FM_MIRIMAGE_DISTORTION_03.02.00.fits

MIRI Imager uses 2 reference files of type:

DISTORTION
FILTEROFFSET

create_miri_imager_wcs_references() creates all reference files.
"""
from __future__ import absolute_import, division, unicode_literals, print_function


import numpy as np
from numpy.testing import assert_allclose
from astropy.io import fits
from astropy.modeling import models
from asdf import AsdfFile


def polynomial_from_coeffs_matrix(coefficients, name=None):
    n_dim = coefficients.ndim

    if n_dim == 1:
        model = models.Polynomial1D(coefficients.size - 1, name=name)
        model.parameters = coefficients
    elif n_dim == 2:
        shape = coefficients.shape
        degree = shape[0] - 1
        if shape[0] != shape[1]:
            raise TypeError("Coefficients must be an (n+1, n+1) matrix")

        coeffs = {}
        for i in range(shape[0]):
            for j in range(shape[0]):
                if i + j < degree + 1:
                    cname = 'c' + str(i) + '_' + str(j)
                    coeffs[cname] = coefficients[i, j]
        model = models.Polynomial2D(degree, name=name, **coeffs)
    return model


def create_miri_imager_filter_offset(distfile, outname):
    """
    Create an asdf reference file with the filter offsets for the MIRI imager.

    Note: The IDT supplied distortion file lists sky to pixel as the
    forward transform. Since "forward" in the JWST pipeline is from
    pixel to sky, the offsets are taken with the opposite sign.

    Parameters
    ----------
    distfile : str
        MIRI imager DISTORTION file provided by the IDT team.
    outname : str
        Name of reference file to be wriiten to disk.

    Returns
    -------
    fasdf : AsdfFile
        AsdfFile object

    Examples
    -------
    >>> create_miri_imager_filer_offset('MIRI_FM_MIRIMAGE_DISTORTION_03.02.00.fits',
                                        'jwst_miri_filter_offset_0001.asdf')
    """

    with fits.open(distfile) as f:
        data = f[9].data

    d = dict.fromkeys(data.field('FILTER'))
    for i in data:
        d[i[0]] = {'column_offset': -i[1], 'row_offset': -i[2]}
    tree = {"title": "MIRI imager filter offset - CDP4",
            "reftype": "FILTEROFFSET",
            "instrument": "MIRI",
            "detector": "MIRIMAGE",
            "pedigree": "GROUND",
            "author": "N. Dencheva",
            "exp_type": "MIR_IMAGE"
            }
    tree.update(d)
    f = AsdfFile()
    f.tree = tree
    f.write_to(outname)


def create_miri_imager_distortion(distfile, outname):
    """
    Create an asdf reference file with all distortion components for the MIRI imager.
    The filter offsets are stored in a sepaate file.

    Note: The IDT supplied distortion file lists sky to pixel as the
    forward transform. Since "forward" in the JWST pipeline is from
    pixel to sky, the meaning of forward and inverse matrices and the order
    in which they are applied is switched.

    The order of operation from pixel to sky is:
    - Apply MI matrix
    - Apply Ai and BI matrices
    - Apply the TI matrix (this gives V2/V3 coordinates)

    Parameters
    ----------
    distfile : str
        MIRI imager DISTORTION file provided by the IDT team.
    outname : str
        Name of reference file to be wriiten to disk.

    Returns
    -------
    fasdf : AsdfFile
        AsdfFile object

    Examples
    --------
    >>> create_miri_imager_distortion("MIRI_FM_MIRIMAGE_DISTORTION_03.02.00.fits", 'test.asdf')
    """
    fdist = fits.open(distfile)
    mi_matrix = fdist[8].data
    mi_col = models.Polynomial1D(1, c0=mi_matrix[0, 2], c1=mi_matrix[0, 0], name="M_column_correction")
    mi_row = models.Polynomial1D(1, c0=mi_matrix[1, 2], c1=mi_matrix[1, 1], name="M_row_correction")
    m_matrix = fdist[4].data
    m_col = models.Polynomial1D(1, c0=m_matrix[0, 2], c1=m_matrix[0, 0])
    m_row = models.Polynomial1D(1, c0=m_matrix[1, 2], c1=m_matrix[1, 1])
    mi_col.inverse = m_col
    mi_row.inverse = m_row
    m_transform = mi_col & mi_row #mi_row & mi_col

    ai_matrix = fdist[6].data
    a_matrix = fdist[2].data
    col_poly = polynomial_from_coeffs_matrix(ai_matrix, name="A_correction")
    col_poly.inverse = polynomial_from_coeffs_matrix(a_matrix)
    bi_matrix = fdist[5].data
    b_matrix = fdist[1].data
    row_poly = polynomial_from_coeffs_matrix(bi_matrix, name="B_correction")
    row_poly.inverse = polynomial_from_coeffs_matrix(b_matrix)
    poly = col_poly & row_poly

    ti_matrix = fdist[7].data
    t_matrix = fdist[3].data
    ti_col = models.Polynomial2D(1, name='TI_column_correction')
    ti_col.parameters = ti_matrix[0][::-1]
    ti_row = models.Polynomial2D(1, name='TI_row_correction')
    ti_row.parameters = ti_matrix[1][::-1]

    t_col = models.Polynomial2D(1, name='T_column_correction')
    t_col.parameters = t_matrix[0][::-1]
    t_row = models.Polynomial2D(1, name='T_row_correction')
    t_row.parameters = t_matrix[1][::-1]
    ti_col.inverse = t_col
    ti_row.inverse = t_row
    t_transform = ti_row & ti_col

    mapping = models.Mapping([0, 1, 0, 1])
    mapping.inverse = models.Identity(2)

    # ident is created here so that mapping can be assigned as inverse
    ident = models.Identity(2)
    ident.inverse = models.Mapping([0, 1, 0, 1])

    poly2t_mapping = models.Mapping([0, 1, 0, 1])
    poly2t_mapping.inverse = models.Mapping([0, 1, 0, 1])

    distortion_transform = m_transform | mapping | poly | poly2t_mapping | t_transform | ident

    fdist.close()
    f = AsdfFile()
    tree = {"title": "MIRI imager distortion - CDP4",
            "reftype": "DISTORTION",
            "instrument": "MIRI",
            "detector": "MIRIMAGE",
            "exp_type": "MIR_IMAGE",
            "pedigree": "GROUND",
            "author": "N. Dencheva",
            "model": distortion_transform
            }
    f.tree = tree
    fasdf = f.write_to(outname)
    #return fasdf


def create_miri_imager_wcs_references(filename, ref):
    """
    Create the two reference files. Writes the files in the current directory.

    Parameters
    ----------
    filename : str
        The name of the IDT file with the distortion.
        In CDP3 the file is called "MIRI_FM_MIRIMAGE_DISTORTION_03.02.00.fits"
    ref : dict
        A dictionary {reftype: refname}, e.g.
        {'DISTORTION': 'jwst_miri_distortion_0001.asdf',
         'FILTEROFFSET': 'jwst_miri_filteroffset_0001.asdf'
         }

    Examples
    --------
    >>> create_miri_imager_wcs_references('MIRI_FM_MIRIMAGE_DISTORTION_03.02.00.fits',
        {'DISTORTION': 'jwst_miri_distortion_0001.asdf',
        'FILTEROFFSET': 'jwst_miri_filter_offset_0001.asdf'})

    """
    try:
        create_miri_imager_distortion(filename, ref['DISTORTION'])
    except:
        print("Distortion file was not created.")
        raise
    try:
        create_miri_imager_filter_offset(filename, ref['FILTEROFFSET'])
    except:
        print("Filter offset file was not created.")
        raise


def test_transform(asdf_file):
    """
    Parameters
    ----------
    asdf_file: str
        reference file with distortion

    xy, v2, v3 values are from technical report with CDP-3 delivery
    """
    v2 = np.array([-2, -2, -2, -2, -1.5, -1.5, -1.5, -1.5, -1, -1, -1, -1], dtype=np.float)
    v3 = np.array([-8, -7.5, -7, -6.5, -8, -7.5, -7, -6.5, -8, -7.5, -7, -6.5], dtype=np.float)
    xy = np.array([[945.80, 728.45], [676.57, 748.63], [408.29, 768.69], [138.02, 789.09],
                   [924.30, 456.59], [655.18, 477.89], [387.05, 498.99], [116.92, 519.96],
                   [904.31, 185.02], [635.09, 207.37], [366.53, 229.45], [95.58, 250.95]],
                  dtype=np.float)

    f = AsdfFile.open(asdf_file)
    transform = f.tree['distortion']
    x, y = transform.inverse(v2, v3)
    assert_allclose(x, xy[:, 0], atol=.05)
    assert_allclose(y, xy[:, 1], atol=.05)
    s1, s2 = transform(xy[:, 0], xy[:, 1])
    assert_allclose(s1, v2, atol=0.05)
    assert_allclose(s2, v3, atol=.05)
