#
#  Module for defning mask geometry in pupil space
#

import numpy as np
import math

from  .utils import rotate2dccw

m = 1.0
mm = 1.0e-3 * m
um = 1.0e-6 * m


class NRM_mask_definitions():

    def __init__(self, maskname=None, rotdeg=None, holeshape="circ", rescale=False,\
                 chooseholes=None):
        """
        Short Summarry
        --------------
        Set attributes of NRM_mask_definitions class.

        Parameters
        ----------
        maskname: string
            name of mask

        rotdeg: list of floats
            range of rotations to search (degrees)

        holeshape: string
            shape of apertures

        rescale: float
            multiplicative factor to adjust hole sizes and centers in entrance pupil

        chooseholes: list
            None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

        Returns
        -------
        None
        """

        if maskname not in ["gpi_g10s40",  "jwst_g7s6", "jwst_g7s6c", "visir_sam", \
                            "p1640", "keck_nirc2", "pharo", "NIRC2_9NRM"]:
            raise ValueError("mask not supported")
        if holeshape is None:
            holeshape = 'circ'

        if holeshape not in ["circ", "hex",]:
            raise ValueError("Unsupported mask holeshape" + maskname)
        self.maskname = maskname

        if self.maskname == "jwst_g7s6c":
            # activeD and D taken from webbpsf-data/NIRISS/coronagraph/MASK_NRM.fits

            self.hdia, self.ctrs = jwst_g7s6c(chooseholes=chooseholes)
            self.activeD = 6.559*m # webbpsf kwd DIAM  - not a 'circle including all holes'
            self.OD = 6.610645669291339*m # Full pupil file size, incl padding, webbpsf kwd PUPLDIAM
            if rotdeg is not None:
                self.rotdeg = rotdeg

        elif self.maskname == "jwst_g7s6":
            pass # not finished

    def showmask(self):
        """
        Short Summarry
        --------------
        Calculate the diameter of the smallest centered circle (D)
        enclosing the live mask area

        Parameters
        ----------

        Returns
        -------
        Diameter of the smallest centered circle

        """
        radii = []
        for ctr in self.ctrs:
            radii.append(math.sqrt(ctr[0]*ctr[0] + ctr[1]*ctr[1]))

        return 2.0*(max(radii) + 0.5*self.hdia)


def jwst_g7s6_centers_asbuilt(chooseholes=None): # was jwst_g7s6_centers_asdesigned
    """
    Short Summarry
    --------------
    Calculate hole centers with appropriate rotation

    Parameters
    ----------
    chooseholes: list
        None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

    Returns
    -------
    Actual hole centers
    """

    holedict = {} # as_built names, C2 open, C5 closed, but as designed coordinates
    # Assemble holes by actual open segment names (as_built).  Either the full mask or the
    # subset-of-holes mask will be V2-reversed after the as_designed centers  are defined
    # Debug orientations with b4,c6,[c2]
    allholes = ('b4','c2','b5','b2','c1','b6','c6')

    #                                              design  built
    holedict['b4'] = [ 0.00000000,  -2.640000]       #B4 -> B4
    holedict['c2'] = [-2.2863100 ,  0.0000000]       #C5 -> C2
    holedict['b5'] = [ 2.2863100 , -1.3200001]       #B3 -> B5
    holedict['b2'] = [-2.2863100 ,  1.3200001]       #B6 -> B2
    holedict['c1'] = [-1.1431500 ,  1.9800000]       #C6 -> C1
    holedict['b6'] = [ 2.2863100 ,  1.3200001]       #B2 -> B6
    holedict['c6'] = [ 1.1431500 ,  1.9800000]       #C1 -> C6

    # as designed MB coordinates (Mathilde Beaulieu, Peter, Anand).
    # as designed: segments C5 open, C2 closed, meters V2V3 per Paul Lightsey def
    # as built C5 closed, C2 open
    #
    # undistorted pupil coords on PM.  These numbers are considered immutable.
    # as designed seg -> as built seg in comments each ctr entry (no distortion)
    if chooseholes: #holes B4 B5 C6 asbuilt for orientation testing
        holelist = []
        for h in allholes:
            if h in chooseholes:
                holelist.append(holedict[h])
        ctrs_asdesigned = np.array( holelist )
    else:
        # the REAL THING - as_designed 7 hole, m in PM space, no distortion
        # ... shape (7,2)
        ctrs_asdesigned = np.array( [
                [ 0.00000000,  -2.640000],       #B4 -> B4  as-designed -> as-built mapping
                [-2.2863100 ,  0.0000000],       #C5 -> C2
                [ 2.2863100 , -1.3200001],       #B3 -> B5
                [-2.2863100 ,  1.3200001],       #B6 -> B2
                [-1.1431500 ,  1.9800000],       #C6 -> C1
                [ 2.2863100 ,  1.3200001],       #B2 -> B6
                [ 1.1431500 ,  1.9800000]    ] ) #C1 -> C6

    # Preserve ctrs.as-designed (treat as immutable)
    # Reverse V2 axis coordinates to close C5 open C2, and others follow suit...
    # preserve cts.as_built  (treat as immutable)
    ctrs_asbuilt = ctrs_asdesigned.copy()

    # create 'live' hole centers in an ideal, orthogonal undistorted xy pupil space,
    # eg maps open hole C5 in as_designed to C2 as_built, eg C4 unaffacted....
    ctrs_asbuilt[:,0] *= -1

    # LG++ rotate hole centers by 90 deg to match MAST o/p DMS PSF with
    # no affine2d transformations 8/2018 AS
    # LG++ The above aligns the hole patern with the hex analytic FT,
    # flat top & bottom as seen in DMS data. 8/2018 AS
    ctrs_asbuilt = rotate2dccw(ctrs_asbuilt, np.pi/2.0) # overwrites attributes

    # create 'live' hole centers in an ideal, orthogonal undistorted xy pupil space,
    return ctrs_asbuilt * m


def jwst_g7s6c(chooseholes=None):
    """
    Short Summarry
    --------------
    Calculate hole centers with appropriate rotation

    Parameters
    ----------
    chooseholes: list
        None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

    Returns
    -------
    ?, Actual hole centers

    """

    f2f = 0.82 * m # m flat to flat
    return f2f, jwst_g7s6_centers_asbuilt(chooseholes=chooseholes)
