#
#  Module for defining mask geometry in pupil space
#  Replace old mask_definitions by reading info in nrm_model into a similar object.
#

import numpy as np
import math

from .utils import rotate2dccw

m = 1.0
mm = 1.0e-3 * m
um = 1.0e-6 * m


class NRM_definition():

	def __init__(self, nrm_model, maskname='jwst_ami', rotdeg=None, chooseholes=None):
		"""
		Set attributes of NRM_definition class.

		Parameters
		----------
		nrm_model: NRMModel
			datamodel containing NRM reference file data
		
		rotdeg: list of floats
			range of rotations to search (degrees)

		chooseholes: list
			None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

		Returns
		-------
		None
		"""
		if maskname not in ['jwst_ami','jwst_g7s6c']:
			raise ValueError("Mask name not supported")

		self.maskname = maskname # there's only one mask but this is used in oifits

		self.read_nrm_model(nrm_model, chooseholes=chooseholes)

		if rotdeg is not None:
			self.rotdeg = rotdeg

	def read_nrm_model(self, nrm_model, chooseholes=None):
		"""

		Calculate hole centers with appropriate rotation,
		set these and other mask geometry details from NRMModel as
		attributes.

		Parameters
		----------
		nrm_model: NRMModel
			datamodel containing NRM reference file data 
		chooseholes: list
			None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

		Returns
		-------
		f2f: float
			flat-to-flat distance of mask holes
		ctrs_asbuilt: array
			Actual hole centers [meters]
		"""

		ctrs_asdesigned = np.array([[nrm_model.x_a1, nrm_model.y_a1],        # B4 -> B4
							[nrm_model.x_a2, nrm_model.y_a2],       # C5 -> C2
							[nrm_model.x_a3, nrm_model.y_a3],       # B3 -> B5
							[nrm_model.x_a4, nrm_model.y_a4],       # B6 -> B2
							[nrm_model.x_a5, nrm_model.y_a5],       # C6 -> C1
							[nrm_model.x_a6, nrm_model.y_a6],        # B2 -> B6
							[nrm_model.x_a7, nrm_model.y_a7]])        # C1 -> C6

		self.hdia = nrm_model.flat_to_flat
		self.activeD = nrm_model.diameter
		self.OD = nrm_model.pupil_circumscribed

		holedict = {}  # as_built names, C2 open, C5 closed, but as designed coordinates
		# Assemble holes by actual open segment names (as_built).  Either the full mask or the
		# subset-of-holes mask will be V2-reversed after the as_designed centers are defined
		# Debug orientations with b4,c6,[c2]
		allholes = ('b4', 'c2', 'b5', 'b2', 'c1', 'b6', 'c6')

		for hole, coords in zip(allholes,ctrs_asdesigned):
			holedict[hole] = coords

		if chooseholes:  # holes B4 B5 C6 asbuilt for orientation testing
			holelist = []
			for h in allholes:
				if h in chooseholes:
					holelist.append(holedict[h])
			ctrs_asdesigned = np.array(holelist)

		ctrs_asbuilt = ctrs_asdesigned.copy()

		# create 'live' hole centers in an ideal, orthogonal undistorted xy pupil space,
		# eg maps open hole C5 in as_designed to C2 as_built, eg C4 unaffected....
		ctrs_asbuilt[:, 0] *= -1

		# LG++ rotate hole centers by 90 deg to match MAST o/p DMS PSF with
		# no affine2d transformations 8/2018 AS
		# LG++ The above aligns the hole pattern with the hex analytic FT,
		# flat top & bottom as seen in DMS data. 8/2018 AS
		ctrs_asbuilt = rotate2dccw(ctrs_asbuilt, np.pi / 2.0)  # overwrites attributes

		# create 'live' hole centers in an ideal, orthogonal undistorted xy pupil space,
		self.ctrs = ctrs_asbuilt

	def showmask(self):
		"""
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
			radii.append(math.sqrt(ctr[0] * ctr[0] + ctr[1] * ctr[1]))

		return 2.0 * (max(radii) + 0.5 * self.hdia)

