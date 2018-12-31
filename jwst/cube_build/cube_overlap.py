""" Routines for creating single band, single exposure IFU Cubes with
the interoplation method = area
"""
import numpy as np
import math
from ..datamodels import dqflags


def find_area_poly(nvertices, xpixel, ypixel):
    """ Find the area of the polygon

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
        area = (xpixel[i] - xmin) * (ypixel[i + 1] - ymin) - (xpixel[i + 1] - xmin) * (ypixel[i] - ymin)
        areaPoly = areaPoly + area

    areaPoly = abs(0.5 * areaPoly)
    return areaPoly
# _____________________________________________________________________________


def find_area_quad(MinX, MinY, Xcorner, Ycorner):
    """ Find the area of an quadrilateral

    Parameters
    ----------
    MinX : float
       Minimum X value
    MinY : float
       Minimum Y value
    Xcorners : numpy.ndarray
       x corner values (use first 4 corners)
    YCorners : numpy.ndarray
       y corner values (use first 4 corners)

    Returns
    -------
    Area
    """
    PX = []
    PY = []

    PX.append(Xcorner[0] - MinX)
    PX.append(Xcorner[1] - MinX)
    PX.append(Xcorner[2] - MinX)
    PX.append(Xcorner[3] - MinX)
    PX.append(PX[0])

    PY.append(Ycorner[0] - MinY)
    PY.append(Ycorner[1] - MinY)
    PY.append(Ycorner[2] - MinY)
    PY.append(Ycorner[3] - MinY)
    PY.append(PY[0])

    Area = 0.5 * ((PX[0] * PY[1] - PX[1] * PY[0]) +
                  (PX[1] * PY[2] - PX[2] * PY[1]) +
                  (PX[2] * PY[3] - PX[3] * PY[2]) +
                  (PX[3] * PY[4] - PX[4] * PY[3]))

    return abs(Area)
# _______________________________________________________________________


def calcCondition(edge, x1, y1, x2, y2, left, right, top, bottom):
    """ Determine if a point is inside a polygon

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
    return 0   # never executed
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
        return (x > left)
    elif edge == CP_RIGHT:
        return (x < right)
    elif edge == CP_BOTTOM:
        return (y > bottom)
    elif edge == CP_TOP:
        return (y < top)
    else:
        return 0
# _______________________________________________________________________


def solve_intersection(edge, x1, y1, x2, y2,
                       left, right, top, bottom):
    """ Finds the intersection of a polygon and rectangular pixel

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
        if(x1 != x2):
            x = x1 + (1.0 / m) * (y - y1)
        else:
            x = x1
    elif edge == CP_TOP:
        y = top
        if(x1 != x2):
            x = x1 + (1.0 / m) * (y - y1)
        else:
            x = x1
    return x, y
# _______________________________________________________________________


def addpoint(x, y, xnew, ynew, nvertices2):
    """ adds a point to vertices of the detector pixel region inside the spaxel

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


def SH_FindOverlap(xcenter, ycenter, xlength, ylength, xp_corner, yp_corner):
    """ Find overlap between pixel and spaxel

    Using the Sutherland_hedgeman Polygon Clipping Algorithm to solve the
    overlap region first clip the x-y detector plane by the cube's xy rectangle
    then find the overlap area

    Parameters
    ---------
    xcenter : float
      center grid point in x dimension for cube (along slice- alpha)
    ycenter : float
      center grid point in y dimension for cube (lambda)
    xlength : float
      width of spaxel in x dimesion (along slice- alpha)
    ylength : float
      width of spaxel in y dimesion (lambda)
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
            condition = calcCondition(i, x1, y1, x2, y2,
                                      left, right, top, bottom)
            x = 0
            y = 0

            if condition == 1:
                x, y = solve_intersection(i, x1, y1, x2, y2,
                                          left, right, top, bottom)
                nvertices2 = addpoint(x, y, xnew, ynew, nvertices2)
                nvertices2 = addpoint(x2, y2, xnew, ynew, nvertices2)

            elif condition == 2:
                nvertices2 = addpoint(x2, y2, xnew, ynew, nvertices2)
            elif condition == 3:
                x, y = solve_intersection(i, x1, y1, x2, y2,
                                          left, right, top, bottom)
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


def match_det2cube(x, y, sliceno, start_slice, input_model, transform,
                   spaxel_flux,
                   spaxel_weight,
                   spaxel_iflux,
                   xcoord, zcoord,
                   crval1, crval3, cdelt1, cdelt3, naxis1, naxis2):
    """ Match detector pixels to output plane in alpha-beta coordinate system

    This routine assumes a 1-1 mapping in beta - slice no.
    This routine assumes the output coordinate systems is alpha-beta
    The user can not change scaling in beta
    Map the corners of the x,y detector values to a cube defined by alpha,beta,
    lambda. In the alpha,lambda plane find the % area of the detector pixel
    which it overlaps with in the cube. For each spaxel record the detector
    pixels that overlap with it - store flux,  % overlap, beta_distance.

    Parameters
    ----------
    x : numpy.ndarray
       x values of pixels in slice
    y : numpy.ndarray
       y values of pixels in slice
    sliceno : int
      slice number
    input_model : datamodel
      input slope model or file
    transform : transform
      wcs transform to transform x,y to alpha,beta, lambda
    spaxel : list
      list of spaxels holding information on each cube pixel.

    Returns
    -------
    spaxel filled in with needed information on overlapping detector pixels
    """
    nxc = len(xcoord)
    nzc = len(zcoord)

    sliceno_use = sliceno - start_slice + 1
    # 1-1 mapping in beta
    yy = sliceno_use - 1

    pixel_dq = input_model.dq[y, x]

    all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['DROPOUT'] +
                 dqflags.pixel['NON_SCIENCE'] +
                 dqflags.pixel['DEAD'] + dqflags.pixel['HOT'] +
                 dqflags.pixel['RC'] + dqflags.pixel['NONLINEAR'])
    # find the location of all the values to reject in cube building
    good_data = np.where((np.bitwise_and(pixel_dq, all_flags) == 0))

    # good data holds the location of pixels we want to map to cube
    x = x[good_data]
    y = y[good_data]

    # center of first pixel, x,y = 1 for Adrian's equations
    # but we want the pixel corners, x,y values passed into this
    # routine to start at 0
    pixel_flux = input_model.data[y, x]

    yy_bot = y
    yy_top = y + 1
    xx_left = x
    xx_right = x + 1

    alpha, beta, lam = transform(x, y)
    alpha1, beta1, lam1 = transform(xx_left, yy_bot)
    alpha2, beta2, lam2 = transform(xx_right, yy_bot)
    alpha3, beta3, lam3 = transform(xx_right, yy_top)
    alpha4, beta4, lam4 = transform(xx_left, yy_top)

    nn = len(x)
    # Loop over all pixels in slice
    for ipixel in range(0, nn - 1):

        # detector pixel -> 4 corners
        # In alpha,wave space
        # in beta space: beta center + width

        alpha_corner = []
        wave_corner = []

        alpha_corner.append(alpha1[ipixel])
        alpha_corner.append(alpha2[ipixel])
        alpha_corner.append(alpha3[ipixel])
        alpha_corner.append(alpha4[ipixel])

        wave_corner.append(lam1[ipixel])
        wave_corner.append(lam2[ipixel])
        wave_corner.append(lam3[ipixel])
        wave_corner.append(lam4[ipixel])

# ________________________________________________________________________________
# Now it does not matter the WCS method used
        alpha_min = min(alpha_corner)
        alpha_max = max(alpha_corner)
        wave_min = min(wave_corner)
        wave_max = max(wave_corner)

        Area = find_area_quad(alpha_min, wave_min, alpha_corner, wave_corner)

        # estimate the where the pixel overlaps in the cube
        # find the min and max values in the cube xcoord,ycoord and zcoord

        MinA = (alpha_min - crval1) / cdelt1
        MaxA = (alpha_max - crval1) / cdelt1
        ix1 = max(0, int(math.trunc(MinA)))
        ix2 = int(math.ceil(MaxA))
        if ix2 >= nxc:
            ix2 = nxc - 1

        MinW = (wave_min - crval3) / cdelt3
        MaxW = (wave_max - crval3) / cdelt3
        iz1 = int(math.trunc(MinW))
        iz2 = int(math.ceil(MaxW))
        if iz2 >= nzc:
            iz2 = nzc - 1

        # loop over possible overlapping cube pixels
        # noverlap = 0
        nplane = naxis1 * naxis2

        for zz in range(iz1, iz2 + 1):
            zcenter = zcoord[zz]
            istart = zz * nplane

            for xx in range(ix1, ix2 + 1):
                cube_index = istart + yy * naxis1 + xx  # yy = slice # -1
                xcenter = xcoord[xx]
                AreaOverlap = SH_FindOverlap(xcenter, zcenter,
                                             cdelt1, cdelt3,
                                             alpha_corner, wave_corner)

                if AreaOverlap > 0.0:
                    AreaRatio = AreaOverlap / Area
                    spaxel_flux[cube_index] = spaxel_flux[cube_index] + \
                        (AreaRatio * pixel_flux[ipixel])
                    spaxel_weight[cube_index] = spaxel_weight[cube_index] + \
                        AreaRatio
                    spaxel_iflux[cube_index] = spaxel_iflux[cube_index] + 1
# ________________________________________________________________________________


class Error2DPolygon(Exception):
    """ Exception raise when no overlap between the detector pixel and
    output plane is found
    """
    pass
