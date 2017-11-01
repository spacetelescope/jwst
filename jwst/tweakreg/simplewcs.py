"""
A module that provides SimpleWCS class - a class that simplifies working with
GWCS objects, in particular as related to WCS transformations sky<->pix
that either include or do not include SIP. Another major purpose
for this class is to allow easy manipulation of WCS parameters related to
standard FITS WCS (CRPIX, CDELT, PC, CRVAL, LONPOLE).

.. warning::
    This class is intended mostly for the internal use by `tweakreg`. This
    class is intended to provide some sort of control in GWCS which, at this
    moment almost completely lacks any kind of standartization. In addition,
    it provides workarounds to some limitations/bugs present in GWCS such as
    the one described here: https://github.com/spacetelescope/gwcs/issues/46

    The API of this class may change in the future as GWCS evolves and bugs
    get fixed and and therefore this class should not be used in external code.

:Authors: Mihai Cara (contact: help@stsci.edu)

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

# STDLIB
import logging
from copy import deepcopy

# THIRD PARTY
import numpy as np
from astropy.modeling.models import (
    Shift, Scale, RotateNative2Celestial, AffineTransformation2D
)
from astropy.modeling.projections import Pix2SkyProjection

import gwcs


__all__ = ['SimpleWCS']

__version__ = '0.8.0'
__vdate__ = '17-April-2016'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class SimpleWCS(object):

    def __init__(self, wcs, copy=False):
        self.check_wcs(wcs)
        if copy:
            self._wcs = deepcopy(wcs)
        else:
            self._wcs = wcs
        self._dirty = False
        self._update_in_progress = False

        # exctract important frames and transformations form the WCS:
        frm1 = self._wcs.available_frames[-2]
        frm2 = self._wcs.available_frames[-1]

        tran = self._wcs.get_transform(frm1, frm2)

        self._shift_x = tran[0].copy()
        self._shift_y = tran[1].copy()
        self._affine = tran[2].copy()
        self._scale_x = tran[3].copy()
        self._scale_y = tran[4].copy()
        self._proj = tran[5].copy()
        self._skyrot = tran[6].copy()

        # estimate pixel scale near CRPIX:
        cd = np.dot(np.diag((self.cdelt1, self.cdelt2)),
                    self._affine.matrix.value)
        self._pscale = self._calc_pixel_scale()

        self.all_pix2world = self._wcs.__call__
        self.all_world2pix = self._wcs.invert

        self.wcs_pix2world = tran
        self.wcs_world2pix = self._wcs.get_transform(frm2, frm1)

    @property
    def wcs(self):
        return self._wcs

    def copy(self):
        return deepcopy(self)

    def check_wcs(self, wcs):
        #TODO: It is not clear (at the time of development) that a WCS will have
        # all the components below. A re-design of how we work with WCS
        # may be needed.

        if len(wcs.available_frames) < 2:
            raise ValueError("Unsupported WCS structure.")

        frm1 = wcs.available_frames[-2]
        frm2 = wcs.available_frames[-1]
        tran = wcs.get_transform(frm1, frm2)
        if len(tran.submodel_names) != 7:
            raise ValueError("Unsupported WCS structure.")

        submod = list(map(tran.__getitem__, tran.submodel_names))

        if (type(submod[0]) != Shift or type(submod[1]) != Shift or
            type(submod[2]) != AffineTransformation2D or
            type(submod[3]) != Scale or type(submod[4]) != Scale or
            not isinstance(submod[5], Pix2SkyProjection) or
            type(submod[6]) != RotateNative2Celestial):
            raise ValueError("Unsupported WCS structure.")

    def begin_update(self):
        self._update_in_progress = True

    def end_update(self):
        if not self._update_in_progress:
            return
        self._update_in_progress = False
        self._update_wcs()

    def _update_wcs(self, force=False):
        if (self._update_in_progress or not self._dirty) and not force:
            return
        #TODO: in the future, pending the outcome of
        # https://github.com/spacetelescope/gwcs/issues/46, we may be able to
        # simplify things. For now we have to re-build the entire transformation
        # for changes to stick.
        frm1 = self._wcs.available_frames[-2]
        frm2 = self._wcs.available_frames[-1]
        new_transform = self._compose_stdwcs_transform()
        self._wcs.set_transform(frm1, frm2, new_transform)
        self._dirty = False

    def _compose_stdwcs_transform(self):
        shift = self._shift_x & self._shift_y
        scale = self._scale_x & self._scale_y
        return (shift | self._affine | scale | self._proj | self._skyrot)

    def strip_nonstd_wcs(self):
        # TODO: once GWCS allows deleting frames "in place" - see
        # https://github.com/spacetelescope/gwcs/issues/48
        # - replace the code below with a code that will remove unneeded
        # frames.
        for frm1, frm2 in zip(self._wcs.available_frames[:-2],
                              self._wcs.available_frames[1:-1]):
            self._wcs.set_transform(frm1, frm2, Identity(2))

        #def remove_frame(wcs, frame_name, tran='before'):
            #npipe_steps = len(wcs.pipeline)
            #if npipe_steps <= 2:
                #raise ValueError("WCS must have more than two frames")
                ## or return unchanged object
            #frmno = [k[0] for k in wcs.pipeline].index(frame_name)
            #if tran == 'before':
                #if frmno == 0:
                    #raise ValueError("Cannot delete first frame when 'tran' is 'before'")
                    ## or return unchanged object
                #prev_frm_name = wcs.pipeline[frmno - 1][0]
                #curr_frm_tran = wcs.pipeline[frmno][1]
                #new_pipe_step = (prev_frm_name, curr_frm_tran)
                #del wcs.pipeline[frmno]
                #del wcs.pipeline[frmno - 1]
                #wcs.pipeline.insert(frmno - 1, new_pipe_step)
            #elif tran == 'after':
                #if frmno == npipe_steps - 1:
                    #raise ValueError("Cannot delete last frame when 'tran' is 'after'")
                    ## or return unchanged object
                #del wcs.pipeline[frmno]
            #else:
                #raise ValueError("Parameter 'tran' must be either 'before' or 'after'")

    def get_std_wcs(self):
        # TODO: this function may become redundant once strip_nonstd_wcs()
        # can remove unneeded frames
        frmname1 = self._wcs.available_frames[0]
        frmname2 = self._wcs.available_frames[-1]
        frm1 = getattr(self._wcs, frmname1)
        frm2 = getattr(self._wcs, frmname2)
        shift = self._shift_x & self._shift_y
        scale = self._scale_x & self._scale_y
        tran = (shift | self._affine | scale | self._proj | self._skyrot)
        return gwcs.wcs.WCS([(frm1, tran), (frm2, None)], name=self._wcs.name)

    @property
    def pscale(self):
        return self._pscale

    @property
    def crpix1(self):
        return -self._shift_x.offset.value

    @crpix1.setter
    def crpix1(self, crpix1):
        self._shift_x.offset = -crpix1
        self._dirty = True
        self._update_wcs()

    @property
    def crpix2(self):
        return -self._shift_y.offset.value

    @crpix2.setter
    def crpix2(self, crpix2):
        self._shift_y.offset = -crpix2
        self._dirty = True
        self._update_wcs()

    @property
    def crpix(self):
        crpix1 = -self._shift_x.offset.value
        crpix2 = -self._shift_y.offset.value
        return np.array([crpix1, crpix2], dtype=np.float64)

    @crpix.setter
    def crpix(self, crpix):
        self._shift_x.offset = -crpix[0]
        self._shift_y.offset = -crpix[1]
        self._dirty = True
        self._update_wcs()

    @property
    def cdelt1(self):
        return self._scale_x.factor.value

    @cdelt1.setter
    def cdelt1(self, cdelt1):
        self._scale_x.factor = cdelt1
        self._dirty = True
        self._update_wcs()

    @property
    def cdelt2(self):
        return self._scale_y.factor.value

    @cdelt2.setter
    def cdelt2(self, cdelt2):
        self._scale_y.factor = cdelt2
        self._dirty = True
        self._update_wcs()

    @property
    def cdelt(self):
        cdelt1 = self._scale_x.factor.value
        cdelt2 = self._scale_y.factor.value
        return np.array([cdelt1, cdelt2], dtype=np.float64)

    @cdelt.setter
    def cdelt(self, cdelt):
        self._scale_x.factor = cdelt[0]
        self._scale_y.factor = cdelt[1]
        self._dirty = True
        self._update_wcs()

    @property
    def lonpole(self):
        return self._skyrot.lon_pole.value

    @lonpole.setter
    def lonpole(self, lonpole):
        self._skyrot.lon_pole = lonpole
        self._dirty = True
        self._update_wcs()

    @property
    def crval1(self):
        return self._skyrot.lon.value

    @crval1.setter
    def crval1(self, crval1):
        self._skyrot.lon = crval1
        self._dirty = True
        self._update_wcs()

    @property
    def crval2(self):
        return self._skyrot.lat.value

    @crval2.setter
    def crval2(self, crval2):
        self._skyrot.lat = crval2
        self._dirty = True
        self._update_wcs()

    @property
    def crval(self):
        crval1 = self._skyrot.lon.value
        crval2 = self._skyrot.lat.value
        return np.array([crval1, crval2], dtype=np.float64)

    @crval.setter
    def crval(self, crval):
        self._skyrot.lon = crval[0]
        self._skyrot.lat = crval[1]
        self._dirty = True
        self._update_wcs()

    @property
    def pc(self):
        return self._affine.matrix.value

    @pc.setter
    def pc(self, pc):
        self._affine.matrix = pc
        self._dirty = True
        self._update_wcs()

    @property
    def cd(self):
        return np.dot(np.diag([self.cdelt1, self.cdelt2]),
                      self._affine.matrix.value)

    def _calc_pixel_scale(self):
        """
        returns pixel scale in arcsec assuming CD is in deg/pix
        """
        # create a wcs that has all transformations except of projection to
        # sky and rotations:
        cd = self.cd
        if len(self._wcs.available_frames) > 2:
            frm1 = self.wcs.available_frames[0]
            frm2 = self.wcs.available_frames[-2]
            tran = self.wcs.get_transform(frm1, frm2)
            d = _deriv(tran, self.crpix1, self.crpix2, dx=2.0, dy=2.0)[0]

            cd = np.dot(cd, d)

        pscale = np.sqrt(np.abs(np.linalg.det(cd))) * 3600.0
        return pscale


def _deriv(tr, x0, y0, dx, dy):
    p = np.asarray([[x0, y0],
                    [x0 - dx, y0],
                    [x0 - dx * 0.5, y0],
                    [x0 + dx * 0.5, y0],
                    [x0 + dx, y0],
                    [x0, y0 - dy],
                    [x0, y0 - dy * 0.5],
                    [x0, y0 + dy * 0.5],
                    [x0, y0 + dy]],
                   dtype=np.float64)

    tp = np.asanyarray(tr(p[:, 0], p[:, 1])).T # transformed points

    # derivative with regard to x:
    u1 = ((tp[1] - tp[4]) + 8 * (tp[3] - tp[2])) / (6 * dx)
    # derivative with regard to y:
    u2 = ((tp[5] - tp[8]) + 8 * (tp[7] - tp[6])) / (6 * dy)

    return (np.asarray([u1, u2]).T, tp[0])
