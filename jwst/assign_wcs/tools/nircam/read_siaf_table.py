import numpy as np
from astropy.io import ascii
from astropy.modeling import models

t = ascii.read("2015 01 27 transform info.CSV")
definitions = t.Row.as_void(t[4])

"""
Based on SIAF transforms PATHS

NIRCAM
------

NIRCAMALW_1 - NIRCAMALW
NIRCAMBLW_1 - NIRCAMBLW
"""

def get_siaf_transform(fromsys, tosys, a_degree, b_degree):
    """
    This reads in the file with transformations that the TEL team
    is using to construct the SIAF file. These transformations are
    defined as teo polynomials (A in x, and B in y) transforming coordinates
    from one coordinate system to another.

    Parameters
    ----------
    fromsys : str
        Starting system
    tosys : str
        Ending coordinate system
    a_degree : int
        Degree of A polynomial
    b_degree : int
        Degree of B polynomial

    Returns
    -------
    a_model : astropy.modeling.Model
        Correction in x
    b_model : astropy.modeling.Model
        Correction in y
    from_units : str
        Units in the starting system
    to_units : str
        Units in the ending system

    Examples
    --------
    >>> get_siaf_transform('NIRCAMALW', "OTESKY", 5, 5)

    """
    index = -1
    col2 = t.columns['col2'].tolist()
    col3 = t.columns['col3'].tolist()
    for i, s in enumerate(zip(col2, col3)):
        if s[0] == fromsys and s[1] == tosys:
            index = i
    if index < 0:
        raise ValueError("Did not find row in CSV table fomr: {0} to {1}".format(fromsys, tosys))
    row = t.Row.as_void(t[index]).tolist()
    from_units = row[3]
    to_units = row[4]
    a_indices = np.array([c.startswith('A') for c in definitions]).nonzero()[0]
    b_indices = np.array([c.startswith('B') for c in definitions]).nonzero()[0]
    a_coeff = np.asarray(row[a_indices[0] : a_indices[-1] + 1], dtype=np.float)
    b_coeff = np.asarray(row[b_indices[0] : b_indices[-1] + 1], dtype=np.float)
    a_model = to_model(a_coeff, a_degree)
    b_model = to_model(b_coeff, b_degree)
    return a_model, b_model, from_units, to_units


def to_model(coeffs, degree=5):
    """
    Creates an astropy.modeling.Model object

    Parameters
    ----------
    coeffs : array like
        Coefficients in the same order as in the ISIM transformations file.
    degree : int
        Degree of polynomial.
        Default is 5 as in the ISIM file but many of the polynomials are of
        a smaller degree.

    Returns
    -------
    poly : astropy.modeling.Polynomial2D
        Polynomial model transforming one coordinate (x or y) between two systems.
    """
    c = {}
    names = []
    for i in range(5):
        for j in range(5):
            if i+j < degree+1:
                names.append('c{0}_{1}'.format(i, j))

    for name, coe in zip(names, coeffs):
        c[name] = coe
    return models.Polynomial2D(degree, **c)

#amodel, bmodel, startunit, endunit = read_siaf_table.get_siaf_transform('NIRCAMALW_1', 'NIRCAMALW', 1, 1)
