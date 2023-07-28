#! /usr/bin/env python

"""

RawOifits class: takes fringefitter class, which contains nrm_list and instrument_data attributes,
all info needed to write oifits.
populate structure needed to write out oifits files according to schema.
averaged and multi-integration versions, sigma-clipped stats over ints


calibrate oifits class, that takes two oifitses and makes the final one?
 questions:
 do we want to retain ability to take median or mean over integrations?

TO DO:
remove 'if verbose' statements, use logging when useful

"""
import numpy as np

from scipy.special import comb
from scipy import stats
from astropy.stats import sigma_clipped_stats

from stdatamodels.jwst import datamodels


class RawOifits:
	def __init__(self, fringefitter, nh=7, angunit="radians", method='median'):
		"""
		Short Summary
		-------------
		Class to store AMI data in the format required to write out to OIFITS files
		Based on ObservablesFromText from ImPlaneIA.
		Angular quantities of input are in radians from fringe fitting; converted to degrees for saving.

		Parameters
		----------
		fringefitter: object, containing nrm_list attribute (list of nrm objects)
		nh: default 7, number of holes in the NRM

		kwargs options:

		"""
		self.ff = fringefitter
		self.nh = nh
		self.nslices = len(self.ff.nrm_list)
		self.nbl = int(comb(self.nh, 2))
		self.ncp = int(comb(self.nh, 3))
		self.nca = int(comb(self.nh, 4))

		self.angunit = angunit
		self.method = method

		if verbose:
			print("ImPlaneIA text output angle unit: %s" % angunit)

		if angunit == 'radians':
			print("Will convert all angular quantities to degrees for saving")
			self.degree = 180.0 / np.pi
		else:
			self.degree = 1

		self.ctrs_eqt = self.ff.instrument_data.ctrs_eqt	
		self.ctrs_inst = self.ff.instrument_data.ctrs_inst
		self.pa = self.ff.instrument_data.pav3 # header pav3, not including v3i_yang??


		self.bholes, self.bls = self._makebaselines()
		self.tholes, self.tuv = self._maketriples_all()
		self.qholes, self.quvw = self._makequads_all()

		self.make_obsarrays()
		self.populate_nrm_dict()
		oimodel = self.populate_oimodel()
		return oimodel

	def make_obsarrays(self):
		"""
		Populate arrays of observables
		"""

		# 3d arrays of centered data, models, and residuals (data-model)
		self.ctrd_arr = np.zeros((self.nslices,self.ff.scidata.shape[1],self.ff.scidata.shape[2]))
		self.n_ctrd_arr = np.zeros((self.nslices,self.ff.scidata.shape[1],self.ff.scidata.shape[2]))
		self.model_arr = np.zeros((self.nslices,self.ff.scidata.shape[1],self.ff.scidata.shape[2]))
		self.n_model_arr = np.zeros((self.nslices,self.ff.scidata.shape[1],self.ff.scidata.shape[2]))
		self.resid_arr = np.zeros((self.nslices,self.ff.scidata.shape[1],self.ff.scidata.shape[2]))
		self.n_resid_arr = np.zeros((self.nslices,self.ff.scidata.shape[1],self.ff.scidata.shape[2]))
		# arrays of observables, (nslices,nobservables) shape.
		self.fp = np.zeros((self.nslices, self.nbl))
		self.fa = np.zeros((self.nslices, self.nbl))
		self.cp = np.zeros((self.nslices, self.ncp))
		self.ca = np.zeros((self.nslices, self.nca))
		self.pistons = np.zeros((self.nslices, self.nh))
		# model parameters
		self.solns = np.zeros((self.nslices,44))

		for i,nrmslc in enumerate(self.ff.nrm_list):
			datapeak = nrmslc.reference.max()
			# NOTE: still have to write this out into datamodel
			self.ctrd_arr[i,:,:] = nrmslc.reference
			self.n_ctrd_arr[i,:,:] = nrmslc.reference/datapeak
			self.model_arr[i,:,:] = nrmslc.modelpsf
			self.n_model_arr[i,:,:] = nrmslc.modelpsf/datapeak
			self.resid_arr[i,:,:] = nrmslc.residual
			self.n_resid_arr[i,:,:] = nrmslc.residual/datapeak
			#
			self.fp[i,:] = np.rad2deg(nrmslc.fringephase) # FPs in degrees
			self.fa[i,:] = nrmslc.fringeamp
			self.cp[i,:] = np.rad2deg(nrmslc.redundant_cps) # CPs in degrees
			self.ca[i,:] = nrmslc.redundant_cas
			self.pistons[i,:] = np.rad2deg(nrmslc.fringepistons) # segment pistons in degrees
			self.solns[i,:] = nrmslc.soln

		self.fa2 = self.fa**2 # squared visibilities

	def populate_nrm_dict(self):
		"""
		Remaining manipulations for OIFITS writing
		Produce a dictionary that will be given to whatever actually writes out the files according to datamodels
		"""
		info4oif = self.ff.instrument_data
		t = Time('%s-%s-%s' %
             (info4oif.year, info4oif.month, info4oif.day), format='fits')
	    ins = info4oif.telname
	    filt = info4oif.filt

		# central wavelength, equiv. width from effstims used for fringe fitting
	    wl, e_wl = info4oif.lam_c, info4oif.lam_c*info4oif.lam_w 

	    # Index 0 and 1 reversed to get the good u-v coverage (same fft)
	    ucoord = self.bls[:, 1]
	    vcoord = self.bls[:, 0]

	    D = 6.5  # Primary mirror display

	    theta = np.linspace(0, 2*np.pi, 100)

	    x = D/2. * np.cos(theta)  # Primary mirror display
	    y = D/2. * np.sin(theta)

	    bl_vis = ((ucoord**2 + vcoord**2)**0.5)

	    tuv = self.tuv
	    v1coord = tuv[:, 0, 0]
	    u1coord = tuv[:, 0, 1]
	    v2coord = tuv[:, 1, 0]
	    u2coord = tuv[:, 1, 1]
	    u3coord = -(u1coord+u2coord)
	    v3coord = -(v1coord+v2coord)

	    bl_cp = []
	    n_bispect = len(v1coord)
	    for k in range(n_bispect):
	        B1 = np.sqrt(u1coord[k] ** 2 + v1coord[k] ** 2)
	        B2 = np.sqrt(u2coord[k] ** 2 + v2coord[k] ** 2)
	        B3 = np.sqrt(u3coord[k] ** 2 + v3coord[k] ** 2)
	        bl_cp.append(np.max([B1, B2, B3]))  # rad-1
	    bl_cp = np.array(bl_cp)

	    flagVis = [False] * self.nbl
	    flagT3 = [False] * self.ncp

	    # do the things done by populate_nrm here
	    # average or don't, and get uncertainties
	    # Unwrap phases
	    shift2pi = np.zeros(self.cp.shape)
	    shift2pi[self.cp >= 6] = 2 * np.pi
	    shift2pi[self.cp <= -6] = -2 * np.pi
	    self.cp -= shift2pi

		if self.method not in ['mean','median','multi']:
			self.method = 'median'
		# set these as attributes (some may exist and be overwritten)
	    if self.method == 'multi':
	        self.vis2 = self.fa2.T
	        self.e_vis2 = np.zeros(vis2.shape)
	        self.visamp = self.fa.T
	        self.e_visamp = np.zeros(visamp.shape)
	        self.visphi = self.fp.T
	        self.e_visphi = np.zeros(visphi.shape)
	        self.cp = self.cp.T
	        self.e_cp = np.zeros(cp.shape)
	        self.camp = self.ca.T
	        self.e_camp = np.zeros(camp.shape)
	        self.pist = self.pistons.T
	        self.e_pist = np.zeros(pist.shape)

	    # apply sigma-clipping to uncertainties
	    # sigma_clipped_stats returns mean, median, stddev. nsigma=3, niters=5
	    elif self.method == 'median':
	    	_, self.vis2, self.e_vis2 = sigma_clipped_stats(self.fa2, axis=0)
	    	_, self.visamp, self.e_visamp = sigma_clipped_stats(self.fa, axis=0)
	    	_, self.visphi, self.e_visphi = sigma_clipped_stats(self.fp, axis=0)
	    	_, self.cp, self.e_cp = sigma_clipped_stats(self.cp, axis=0)
	    	_, self.camp, self.e_camp = sigma_clipped_stats(self.ca, axis=0)
	    	_, self.pist, self.e_pist = sigma_clipped_stats(self.pistons, axis=0)

	    else: # take the mean
	    	self.vis2, _, self.e_vis2 = sigma_clipped_stats(self.fa2, axis=0)
	    	self.visamp, _, self.e_visamp = sigma_clipped_stats(self.fa, axis=0)
	    	self.visphi, _, self.e_visphi = sigma_clipped_stats(self.fp, axis=0)
	    	self.cp, _, self.e_cp = sigma_clipped_stats(self.cp, axis=0)
	    	self.camp, _, self.e_camp = sigma_clipped_stats(self.ca, axis=0)
	    	self.pist, _, self.e_pist = sigma_clipped_stats(self.pistons, axis=0)

		# prepare arrays for OI_ARRAY ext
		self.staxy = info4oif.ctrs_inst
		N_ap = len(staxy)
	    tel_name = ['A%i' % x for x in np.arange(N_ap)+1]
	    sta_name = tel_name
	    diameter = [0] * N_ap

	    staxyz = []
	    for x in staxy:
	        a = list(x)
	        line = [a[0], a[1], 0]
	        staxyz.append(line)

    	sta_index = np.arange(N_ap) + 1

    	pscale = info4oif.pscale/1000.  # arcsec
	    isz = info4oif.isz  # Size of the image to extract NRM data
	    fov = [pscale * isz] * N_ap
	    fovtype = info4oif.radius * N_ap

	    # prepare info for OI_TARGET ext
	    name_star = info4oif.objname

	    customSimbad = Simbad()
	    customSimbad.add_votable_fields('propermotions', 'sptype', 'parallax')

	    # Add informations from Simbad:
	    if name_star == 'UNKNOWN':
	        ra, dec, spectyp = [0], [0], ['unknown']
	        pmra, pmdec, plx = [0], [0], [0]
	    else:
	        try:
	            query = customSimbad.query_object(name_star)
	            coord = SkyCoord(query['RA'][0]+' '+query['DEC']
	                             [0], unit=(u.hourangle, u.deg))
	            ra, dec = [coord.ra.deg], [coord.dec.deg]
	            spectyp, plx = query['SP_TYPE'], query['PLX_VALUE']
	            pmra, pmdec = query['PMRA'], query['PMDEC']
	        except TypeError:
	            ra, dec, spectyp = [0], [0], ['unknown']
	            pmra, pmdec, plx = [0], [0], [0]


		self.oifits_dct = {
			'OI_VIS2': {'VIS2DATA': self.vis2,
                       'VIS2ERR': self.e_vis2,
                       'UCOORD': ucoord,
                       'VCOORD': vcoord,
                       'STA_INDEX': self.bholes,
                       'MJD': t.mjd,
                       'INT_TIME': info4oif.itime,
                       'TIME': 0,
                       'TARGET_ID': 1,
                       'FLAG': flagVis,
                       'BL': bl_vis
                       },

           'OI_VIS': {'TARGET_ID': 1,
                      'TIME': 0,
                      'MJD': t.mjd,
                      'INT_TIME': info4oif.itime,
                      'VISAMP': self.visamp,
                      'VISAMPERR': self.e_visamp,
                      'VISPHI': self.visphi,
                      'VISPHIERR': self.e_visphi,
                      'UCOORD': ucoord,
                      'VCOORD': vcoord,
                      'STA_INDEX': self.bholes,
                      'FLAG': flagVis,
                      'BL': bl_vis
                      },

           'OI_T3': {'TARGET_ID': 1,
                     'TIME': 0,
                     'MJD': t.mjd,
                     'INT_TIME': info4oif.itime,
                     'T3PHI': self.cp,
                     'T3PHIERR': self.e_cp,
                     'T3AMP': self.cpamp,
                     'T3AMPERR': self.e_cp,
                     'U1COORD': u1coord,
                     'V1COORD': v1coord,
                     'U2COORD': u2coord,
                     'V2COORD': v2coord,
                     'STA_INDEX': self.tholes,
                     'FLAG': flagT3,
                     'BL': bl_cp
                     },

           'OI_WAVELENGTH': {'EFF_WAVE': wl,
                             'EFF_BAND': e_wl
                             },
            'OI_ARRAY': {'TEL_NAME': tel_name,
            			'STA_NAME': sta_name,
            			'STA_INDEX': sta_index,
            			'DIAMETER': diameter,
            			'STAXYZ': staxyz,
            			'FOV': fov,
            			'FOVTYPE': fovtype,
            			'CTRS_EQT': info4oif.ctrs_eqt, # mask hole coords rotated to equatotial
            			'PISTONS': self.pist, # RAC 2021
                    	'PIST_ERR': self.e_pist
            			},
            'OI_TARGET': {'TARGET_ID':[1],
            			'TARGET': info4oif.objname,
            			'RAEP0': ra,
            			'DECEP0': dec,
            			'EQUINOX': [2000],
            			'RA_ERR': [0],
            			'DEC_ERR': [0],
            			'SYSVEL': [0],
            			'VELTYP': ['UNKNOWN'],
            			'VELDEF': ['OPTICAL'],
            			'PMRA': pmra,
            			'PMDEC': pmdec,
            			'PMRA_ERR': [0],
            			'PMDEC_ERR': [0],
            			'PARALLAX': plx,
            			'PARA_ERR': [0],
            			'SPECTYP': spectyp
            			},

           'info': {'TARGET': info4oif.objname,
                    'CALIB': info4oif.objname,
                    'OBJECT': info4oif.objname,
                    'PROPNAME': info4oif.proposer_name,
                    'FILT': info4oif.filt,
                    'INSTRUME': info4oif.instrument,
                    'ARRNAME': info4oif.arrname,
                    'MASK': info4oif.arrname, # oifits.py looks for dct.info['MASK']
                    'MJD': t.mjd,
                    'DATE-OBS': t.fits,
                    'TELESCOP': info4oif.telname,
                    'OBSERVER': info4oif.pi_name,
                    'INSMODE': info4oif.pupil,
                    'PSCALE': info4oif.pscale_mas,
                    'STAXY': info4oif.ctrs_inst, # as-built mask hole coords
                    'ISZ': self.ff.scidata.shape[1],  # size of the image needed (or fov)
                    'NFILE': 0,
                    'ARRAYX': float(0),
                    'ARRAYY': float(0),
                    'ARRAYZ': float(0),
                    'PA': info4oif.pa,
                    'FRAME': 'SKY',
                    'OI_REVN': 2
                    }
           }

	

	def populate_oimodel(self):
		m = datamodels.AmiOIModel()

		# primary header keywords
	    m.meta.telescope = self.oifits_dct['info']['TELESCOP']
	    m.meta.origin = 'STScI'
	    m.meta.instrument.name = self.oifits_dct['info']['INSTRUME']
	    m.meta.program.pi_name = self.oifits_dct['info']['OBSERVER']
	    m.meta.target.proposer_name = self.oifits_dct['info']['PROPNAME']
	    m.meta.observation.date = self.oifits_dct['info']['DATE-OBS']
	    m.meta.oifits.array_name = self.oifits_dct['info']['MASK']
	    m.meta.oifits.instrument_mode = self.oifits_dct['info']['INSMODE']

	    # oi_array extension data
	    m.array = np.asarray([
	    		self.oifits_dct['OI_ARRAY']['TEL_NAME'],
	    		self.oifits_dct['OI_ARRAY']['STA_NAME']
	    		self.oifits_dct['OI_ARRAY']['STA_INDEX'],
	    		self.oifits_dct['OI_ARRAY']['DIAMTER'],
	    		self.oifits_dct['OI_ARRAY']['STAXYZ'],
	    		self.oifits_dct['OI_ARRAY']['FOV'],
	    		self.oifits_dct['OI_ARRAY']['FOVTYPE'],
	    		self.oifits_dct['OI_ARRAY']['CTRS_EQT'],
	    		self.oifits_dct['OI_ARRAY']['PISTONS'],
	    		self.oifits_dct['OI_ARRAY']['PIST_ERR'],
	    		]).T

	    # oi_target extension data
	    m.target = np.asarray([
	    		self.oifits_dct['OI_TARGET']['TARGET_ID'],
	    		self.oifits_dct['OI_TARGET']['TARGET'],
	    		self.oifits_dct['OI_TARGET']['RAEP0'],
	    		self.oifits_dct['OI_TARGET']['DECEP0'],
	    		self.oifits_dct['OI_TARGET']['EQUINOX'],
	    		self.oifits_dct['OI_TARGET']['RA_ERR'],
	    		self.oifits_dct['OI_TARGET']['DEC_ERR'],
	    		self.oifits_dct['OI_TARGET']['SYSVEL'],
	    		self.oifits_dct['OI_TARGET']['VELTYP'],
	    		self.oifits_dct['OI_TARGET']['VELDEF'],
	    		self.oifits_dct['OI_TARGET']['PMRA'],
	    		self.oifits_dct['OI_TARGET']['PMDEC'],
	    		self.oifits_dct['OI_TARGET']['PMRA_ERR'],
	    		self.oifits_dct['OI_TARGET']['PMDEC_ERR'],
	    		self.oifits_dct['OI_TARGET']['PARALLAX'],
	    		self.oifits_dct['OI_TARGET']['PARA_ERR'],
	    		self.oifits_dct['OI_TARGET']['SPECTYP'],
	    		]).T
	    # oi_vis extension data
	    m.vis = np.asarray([
	    		self.oifits_dct['OI_VIS']['TARGET_ID'],
	    		self.oifits_dct['OI_VIS']['TIME'],
	    		self.oifits_dct['OI_VIS']['MJD'],
	    		self.oifits_dct['OI_VIS']['INT_TIME'],
	    		self.oifits_dct['OI_VIS']['VISAMP'],
	    		self.oifits_dct['OI_VIS']['VISAMPERR'],
	    		self.oifits_dct['OI_VIS']['VISPHI'],
	    		self.oifits_dct['OI_VIS']['VISPHIERR'],
	    		self.oifits_dct['OI_VIS']['UCOORD'],
	    		self.oifits_dct['OI_VIS']['VCOORD'],
	    		self.oifits_dct['OI_VIS']['STA_INDEX'],
	    		self.oifits_dct['OI_VIS']['FLAG'],
	    		]).T
	    # oi_vis2 extension data
	    m.vis2 = np.asarray([
	    		self.oifits_dct['OI_VIS2']['TARGET_ID'],
	    		self.oifits_dct['OI_VIS2']['TIME'],
	    		self.oifits_dct['OI_VIS2']['MJD'],
	    		self.oifits_dct['OI_VIS2']['INT_TIME'],
	    		self.oifits_dct['OI_VIS2']['VIS2DATA'],
	    		self.oifits_dct['OI_VIS2']['VIS2ERR'],
	    		self.oifits_dct['OI_VIS2']['UCOORD'],
	    		self.oifits_dct['OI_VIS2']['VCOORD'],
	    		self.oifits_dct['OI_VIS2']['STA_INDEX'],
	    		self.oifits_dct['OI_VIS2']['FLAG'],
	    		]).T
	    # oi_t3 extension data
	    m.t3 = np.asarray([
	    		self.oifits_dct['OI_T3']['TARGET_ID'],
	    		self.oifits_dct['OI_T3']['TIME'],
	    		self.oifits_dct['OI_T3']['MJD'],
	    		self.oifits_dct['OI_T3']['T3AMP'],
	    		self.oifits_dct['OI_T3']['T3AMPERR'],
	    		self.oifits_dct['OI_T3']['T3PHI'],
	    		self.oifits_dct['OI_T3']['T3PHIERR'],
	    		self.oifits_dct['OI_T3']['U1COORD'],
	    		self.oifits_dct['OI_T3']['V1COORD'],
	    		self.oifits_dct['OI_T3']['U2COORD'],
	    		self.oifits_dct['OI_T3']['V2COORD'],
	    		self.oifits_dct['OI_T3']['STA_INDEX'],
	    		self.oifits_dct['OI_T3']['FLAG'],
	    		]).T
	    # oi_wavelength extension data
	    m.wavelength = np.asarray([
	    		self.oifits_dct['OI_WAVELENGTH']['EFF_WAVE'],
	    		self.oifits_dct['OI_WAVELENGTH']['EFF_BAND'],
	    		]).T

	    return m

	def _maketriples_all(self):
		""" returns int array of triple hole indices (0-based), 
			and float array of two uv vectors in all triangles
		"""
		nholes = self.ctrs_eqt.shape[0]
		tlist = []
		for i in range(nholes):
			for j in range(nholes):
				for k in range(nholes):
					if i < j and j < k:
						tlist.append((i, j, k))
		tarray = np.array(tlist).astype(np.int)
		if self.verbose:
			print("tarray", tarray.shape, "\n", tarray)

		tname = []
		uvlist = []
		# foreach row of 3 elts...
		for triple in tarray:
			tname.append("{0:d}_{1:d}_{2:d}".format(
				triple[0], triple[1], triple[2]))
			if self.verbose:
				print('triple:', triple, tname[-1])
			uvlist.append((self.ctrs_eqt[triple[0]] - self.ctrs_eqt[triple[1]],
						   self.ctrs_eqt[triple[1]] - self.ctrs_eqt[triple[2]]))
		# print(len(uvlist), "uvlist", uvlist)
		if self.verbose:
			print(tarray.shape, np.array(uvlist).shape)
		return tarray, np.array(uvlist)

	def _makebaselines(self):
		"""
		ctrs_eqt (nh,2) in m
		returns np arrays of eg 21 baselinenames ('0_1',...), eg (21,2) baselinevectors (2-floats)
		in the same numbering as implaneia
		"""
		nholes = self.ctrs_eqt.shape[0]
		blist = []
		for i in range(nholes):
			for j in range(nholes):
				if i < j:
					blist.append((i, j))
		barray = np.array(blist).astype(np.int)
		# blname = []
		bllist = []
		for basepair in blist:
			# blname.append("{0:d}_{1:d}".format(basepair[0],basepair[1]))
			baseline = self.ctrs_eqt[basepair[0]] - self.ctrs_eqt[basepair[1]]
			bllist.append(baseline)
		return barray, np.array(bllist)





class CalibOifits:
	def __init__(self,targoifits,caloifits):
		"""
		Short Summary
		-------------
		

		Parameters
		----------
		targoifits: oifits dataodel, target
		caloifits: oifits datamodel, reference star

		kwargs options:

		"""
	def calibrate():
		cp_out = dct_t['OI_T3']['T3PHI'] - dct_c['OI_T3']['T3PHI']
	    sqv_out = dct_t['OI_VIS2']['VIS2DATA'] / dct_c['OI_VIS2']['VIS2DATA']
	    va_out = dct_t['OI_VIS']['VISAMP'] / dct_c['OI_VIS']['VISAMP']
	    # now using correct propagation of error for multiplication/division
	    # which assumes uncorrelated Gaussian errors (not true...?)    
	    cperr_t = dct_t['OI_T3']['T3PHIERR']
	    cperr_c = dct_c['OI_T3']['T3PHIERR']
	    sqverr_c = dct_t['OI_VIS2']['VIS2ERR']
	    sqverr_t = dct_c['OI_VIS2']['VIS2ERR']
	    vaerr_t = dct_t['OI_VIS']['VISAMPERR']
	    vaerr_c = dct_c['OI_VIS']['VISAMPERR']
	    cperr_out = np.sqrt(cperr_t**2. + cperr_c**2.)
	    sqverr_out = sqv_out * np.sqrt((sqverr_t/dct_t['OI_VIS2']['VIS2DATA'])**2. + (sqverr_c/dct_c['OI_VIS2']['VIS2DATA'])**2.)
	    vaerr_out = va_out * np.sqrt((vaerr_t/dct_t['OI_VIS']['VISAMP'])**2. + (vaerr_c/dct_c['OI_VIS']['VISAMP'])**2.)

	    # copy the target dict and modify with the calibrated observables
	    calib_dict = dct_t.copy()
	    calib_dict['OI_T3']['T3PHI'] = cp_out
	    calib_dict['OI_VIS2']['VIS2DATA'] = sqv_out
	    calib_dict['OI_VIS']['VISAMP'] = va_out
	    calib_dict['OI_T3']['T3PHIERR'] = cperr_out
	    calib_dict['OI_VIS2']['VIS2ERR'] = sqverr_out
	    calib_dict['OI_VIS']['VISAMPERR'] = vaerr_out
	    # preserve the name of the calibrator star
	    calib_dict['info']['CALIB'] = dct_c['info']['OBJECT']
	    # include pistons and piston errors from target and calibrator
	    # if old files, raw oifits won't have any pistons
	    if ('PISTONS' in dct_t['OI_ARRAY']) & ('PISTONS' in dct_c['OI_ARRAY']):
	        pistons_t = dct_t['OI_ARRAY']['PISTONS']
	        pisterr_t = dct_t['OI_ARRAY']['PIST_ERR']
	        pistons_c = dct_c['OI_ARRAY']['PISTONS']
	        pisterr_c = dct_c['OI_ARRAY']['PIST_ERR']
	        # sum in quadrature errors from target and calibrator pistons (only if both oifits contain pistons)
	        pisterr_out = np.sqrt(pisterr_t**2 + pisterr_c**2)
	        # populate calibrated dict with pistons 
	        calib_dict['OI_ARRAY']['PISTON_T'] = pistons_t
	        calib_dict['OI_ARRAY']['PISTON_C'] = pistons_c
	        calib_dict['OI_ARRAY']['PIST_ERR'] = pisterr_out
	    # remove plain "pistons" key from dict
	    if 'PISTONS' in calib_dict['OI_ARRAY']:
	        del calib_dict['OI_ARRAY']['PISTONS']

	    return calib_dict


