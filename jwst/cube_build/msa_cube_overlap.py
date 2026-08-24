"""
Temporary routines for build MSA cubes using python methods.

These routines are only needed until readnoise weighting is part of the C methods.
"""

import numpy as np


def find_area_quad(minx, miny, Xcorner, Ycorner):
    """
    Find the area of an quadrilateral between clipped polygon and cube spaxel.

    Parameters
    ----------
     minx : float
        Minimum X value
     miny : float
        Minimum Y value
     Xcorner : ndarray
        X corner values of a cube spaxel
     YCorner : ndarray
        Y corner values of a cube spaxe

    Returns
    -------
     Area : double
        Area of the overlap
    """
    Area = 0

    PX = np.zeros(5)
    PY = np.zeros(5)
    PX[0] = Xcorner[0] - minx
    PX[1] = Xcorner[1] - minx
    PX[2] = Xcorner[2] - minx
    PX[3] = Xcorner[3] - minx
    PX[4] = PX[0]

    PY[0] = Ycorner[0] - miny
    PY[1] = Ycorner[1] - miny
    PY[2] = Ycorner[2] - miny
    PY[3] = Ycorner[3] - miny
    PY[4] = PY[0]

    Area = 0.5 * (
        (PX[0] * PY[1] - PX[1] * PY[0])
        + (PX[1] * PY[2] - PX[2] * PY[1])
        + (PX[2] * PY[3] - PX[3] * PY[2])
        + (PX[3] * PY[4] - PX[4] * PY[3])
    )

    return np.abs(Area)


def find_volume(area, dwave):

    vol = 0.0

    vol = area * dwave
    return vol


def find_area_poly(nvertices, xpixel, ypixel):
    """Find the area of the polygon

    Parameters
    ----------
    nvertices : int
        number of Vertices of polygon
    xpixel : numpy.ndarray
      x coordinate of vertices
    ypixel : numpy.ndarray
      y coordinate of vertices

    Returns
    -------
    area of polygon

    """
    areaPoly = 0.0
    xmin = min(xpixel)
    ymin = min(ypixel)

    for i in range(0, nvertices - 1):
        area = (xpixel[i] - xmin) * (ypixel[i + 1] - ymin) - (xpixel[i + 1] - xmin) * (
            ypixel[i] - ymin
        )
        areaPoly = areaPoly + area

    areaPoly = abs(0.5 * areaPoly)
    return areaPoly


# _____________________________________________________________________________


def calcCondition(edge, x1, y1, x2, y2, left, right, top, bottom):
    """Determine if a point is inside a polygon

    Parameters
    ----------
    edge : float
      edge of spaxel
    x1 : float
      x min coordinate of pixel
    y1 : float
      y min coordinate of pixel
    x2 : float
      x max coordinate of pixel
    y2 : float
      y max coordinate of pixel
    left : float
      left side of spaxel
    right : float
      right side of spaxel
    top : float
      top  of spaxel
    bottom : float
      bottom of spaxel

    Returns
    -------
    where the detector pixel is in relation to a side of the spaxel
    """
    stat1 = insideWindow(edge, x1, y1, left, right, top, bottom)
    stat2 = insideWindow(edge, x2, y2, left, right, top, bottom)

    if not stat1 and stat2:
        return 1
    if stat1 and stat2:
        return 2
    if stat1 and not stat2:
        return 3
    if not stat1 and not stat2:
        return 4
    return 0  # never executed


# _______________________________________________________________________


def insideWindow(edge, x, y, left, right, top, bottom):
    """Function used in determined overlap of detector pixel and spaxel

    Given the pixel edge and cener  and the left,right,top bottom of
    sides of spaxel return if detector edge is inside spaxel

    Parameters
    ----------
    edge : float
      edge of pixel
    x : float
      x center of pixel
    y : float
      y center of pixel
    left : float
      left side of spaxel
    right : float
      right side of spaxel
    top : float
      top of spaxel
    bottom : float
      bottom of spaxel

    Returns
    -------
    returns true or false values if detector point is on "correct"
    side of spaxel
    """
    CP_LEFT = 0
    CP_RIGHT = 1
    CP_BOTTOM = 2
    CP_TOP = 3

    if edge == CP_LEFT:
        return x > left
    elif edge == CP_RIGHT:
        return x < right
    elif edge == CP_BOTTOM:
        return y > bottom
    elif edge == CP_TOP:
        return y < top
    else:
        return 0


def solve_intersection(edge, x1, y1, x2, y2, left, right, top, bottom):
    """Finds the intersection of a polygon and rectangular pixel

    Parameters
    ----------
    edge : float
      one of the 4 edges of spaxel
    x1 : float
      x min of pixel
    y1 : float
      y min of pixel
    x2 : float
      x max of pixel
    y2 : float
      y max of pixel
    left : float
      left side of spaxel
    right : float
      right side of spaxel
    top : float
      top of spaxel
    bottom : float
      bottom of spaxel

    Returns
    -------
    returns x,y inside spaxel (detector detector region inside spaxel)
    """
    x = 0
    y = 0
    CP_LEFT = 0
    CP_RIGHT = 1
    CP_BOTTOM = 2
    CP_TOP = 3
    m = 0
    if x2 != x1:
        m = (y2 - y1) / (x2 - x1)
    if edge == CP_LEFT:
        x = left
        y = y1 + m * (x - x1)
    elif edge == CP_RIGHT:
        x = right
        y = y1 + m * (x - x1)
    elif edge == CP_BOTTOM:
        y = bottom
        if x1 != x2:
            x = x1 + (1.0 / m) * (y - y1)
        else:
            x = x1
    elif edge == CP_TOP:
        y = top
        if x1 != x2:
            x = x1 + (1.0 / m) * (y - y1)
        else:
            x = x1
    return x, y


# _______________________________________________________________________


def addpoint(x, y, xnew, ynew, nvertices2):
    """Adds a point to vertices of the detector pixel region inside the spaxel

    Parameters
    ----------
    x : float
     x value to add
    y : float
     y value  to add
    xnew : numpy.ndarray
      new x vertices
    ynew : numpy.ndarray
      new y vertices
    nvertices2 : int
      number of vertices

    Returns
    -------
    adds a vertice to the polygon describing the region of the detector pixel
    inside the spaxel
    """
    xnew[nvertices2] = x
    ynew[nvertices2] = y

    nvertices2 = nvertices2 + 1

    return nvertices2


# ________________________________________________________________________________


def sh_find_overlap(xcenter, ycenter, xlength, ylength, xp_corner, yp_corner):
    """Find overlap between pixel and spaxel

    Using the Sutherland_hedgeman Polygon Clipping Algorithm to solve the
    overlap region first clip the x-y detector plane by the cube's xy rectangle
    then find the overlap area

    Parameters
    ----------
    xcenter : float
      center grid point in x dimension for cube (along slice- alpha)
    ycenter : float
      center grid point in y dimension for cube (lambda)
    xlength : float
      width of spaxel in x dimension (along slice- alpha)
    ylength : float
      width of spaxel in y dimension (lambda)
    xp_corner : float
      alpha pixel corner values
    yp_Corner : float
      lambda pixel corner values

    Returns
    -------
    AreaOverlap
    """
    area_clipped = 0.0
    top = ycenter + 0.5 * ylength
    bottom = ycenter - 0.5 * ylength

    left = xcenter - 0.5 * xlength
    right = xcenter + 0.5 * xlength

    nvertices = 4  # input detector pixel vertices
    max_vertices = 9
    # initialize xPixel, yPixel to the detector pixel corners.
    # xPixel,yPixel will become the clipped polygon vertices
    # inside the cube pixel
    # xnew,ynew xpixel and ypixel of size max_vertices

    xPixel = []
    yPixel = []

    xnew = []
    ynew = []

    for j in range(0, 9):
        xnew.append(0.0)
        ynew.append(0.0)
        xPixel.append(0.0)
        yPixel.append(0.0)

    # Xpixel, YPixel closed (5 corners)
    for i in range(0, 4):
        xPixel[i] = xp_corner[i]
        yPixel[i] = yp_corner[i]
    xPixel[4] = xp_corner[0]
    yPixel[4] = yp_corner[0]

    for i in range(0, 4):  # 0:left, 1: right, 2: bottom, 3: top
        nvertices2 = 0
        for j in range(0, nvertices):
            x1 = xPixel[j]
            y1 = yPixel[j]
            x2 = xPixel[j + 1]
            y2 = yPixel[j + 1]
            condition = calcCondition(i, x1, y1, x2, y2, left, right, top, bottom)
            x = 0
            y = 0

            if condition == 1:
                x, y = solve_intersection(i, x1, y1, x2, y2, left, right, top, bottom)
                nvertices2 = addpoint(x, y, xnew, ynew, nvertices2)
                nvertices2 = addpoint(x2, y2, xnew, ynew, nvertices2)

            elif condition == 2:
                nvertices2 = addpoint(x2, y2, xnew, ynew, nvertices2)
            elif condition == 3:
                x, y = solve_intersection(i, x1, y1, x2, y2, left, right, top, bottom)
                nvertices2 = addpoint(x, y, xnew, ynew, nvertices2)

        # condition ==  4: points outside
        # Done looping over J  corners
        nvertices2 = addpoint(xnew[0], ynew[0], xnew, ynew, nvertices2)

        if nvertices2 > max_vertices:
            raise Error2DPolygon(" Failure in finding the clipped polygon, nvertices2 > 9 ")

        nvertices = nvertices2 - 1

        for k in range(0, nvertices2):
            xPixel[k] = xnew[k]
            yPixel[k] = ynew[k]

    # done loop over top,bottom,left,right
    nvertices = nvertices + 1

    if nvertices > 0:
        area_clipped = find_area_poly(nvertices, xPixel, yPixel)

    return area_clipped


# _____________________________________________________________________________


class Error2DPolygon(Exception):
    """Exception raise when no overlap between the detector pixel and
    output plane is found
    """

    pass
