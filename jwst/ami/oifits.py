#! /usr/bin/env python

"""

RawOifits class: takes fringefitter class, which contains nrm_list and instrument_data attributes,
all info needed to write oifits.
populate structure needed to write out oifits files according to schema.
averaged and multi-integration versions, sigma-clipped stats over ints

CalibOifits class: takes two AmiOIModel datamodels and produces a final calibrated datamodel.
"""
import numpy as np

from scipy.special import comb
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.time.core import Time
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord

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
        fringefitter: FringeFitter object
            Object containing nrm_list attribute (list of nrm objects) 
            and other info needed for OIFITS files
        nh: integer
            default 7, number of holes in the NRM
        angunit: string
            Unit of incoming angular quantities (radians or degrees). 
            Will be converted to degrees for OIFITS. Default radians
        method: string
            Method to average observables: mean or median. Default median.


        """
        self.ff = fringefitter
        self.nh = nh # 7
        self.nslices = len(self.ff.nrm_list) # n ints
        self.nbl = int(comb(self.nh, 2)) # 21
        self.ncp = int(comb(self.nh, 3)) # 35
        self.nca = int(comb(self.nh, 4))

        self.angunit = angunit
        self.method = method


        if angunit == 'radians':
            self.degree = 180.0 / np.pi
        else:
            self.degree = 1

        self.ctrs_eqt = self.ff.instrument_data.ctrs_eqt    
        self.ctrs_inst = self.ff.instrument_data.ctrs_inst
        self.pa = self.ff.instrument_data.pav3 # header pav3, not including v3i_yang??


        self.bholes, self.bls = self._makebaselines()
        self.tholes, self.tuv = self._maketriples_all()
        #self.qholes, self.quvw = self._makequads_all()

        self.make_oifits()

    def make_oifits(self):
        """
        Short Summary
        ------------
        Calls the functions to prepare data for saving and 
        populate AmiOIModel.

        Parameters
        ----------

        Returns
        -------
        oimodel: AmiOIModel object
            Datamodel containing AMI observables
        """
        self.make_obsarrays()
        self.populate_nrm_dict()
        oimodel = self.populate_oimodel()

        return oimodel

    def make_obsarrays(self):
        """
        Short Summary
        ------------
        Make arrays of observables of the correct shape for saving to datamodels

        Parameters
        ----------

        Returns
        -------
        """
        # arrays of observables, (nslices,nobservables) shape.
        self.fp = np.zeros((self.nslices, self.nbl))
        self.fa = np.zeros((self.nslices, self.nbl))
        self.cp = np.zeros((self.nslices, self.ncp))
        self.ca = np.zeros((self.nslices, self.nca))
        self.pistons = np.zeros((self.nslices, self.nh))
        # model parameters
        self.solns = np.zeros((self.nslices,44))

        for i,nrmslc in enumerate(self.ff.nrm_list):
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
        Short Summary
        ------------
        Do all remaining manipulations for OIFITS writing
        Produce a dictionary containing all data to write
        to OIFITS files

        Parameters
        ----------

        Returns
        -------
        
        """
        info4oif = self.ff.instrument_data
        t = Time('%s-%s-%s' %
             (info4oif.year, info4oif.month, info4oif.day), format='fits')

        # central wavelength, equiv. width from effstims used for fringe fitting
        wl, e_wl = info4oif.lam_c, info4oif.lam_c*info4oif.lam_w 

        # Index 0 and 1 reversed to get the good u-v coverage (same fft)
        ucoord = self.bls[:, 1]
        vcoord = self.bls[:, 0]

        # D = 6.5  # Primary mirror display

        # theta = np.linspace(0, 2*np.pi, 100)

        # x = D/2. * np.cos(theta)  # Primary mirror display
        # y = D/2. * np.sin(theta)

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
            self.e_vis2 = np.zeros(self.vis2.shape)
            self.visamp = self.fa.T
            self.e_visamp = np.zeros(self.visamp.shape)
            self.visphi = self.fp.T
            self.e_visphi = np.zeros(self.visphi.shape)
            self.cp = self.cp.T
            self.e_cp = np.zeros(self.cp.shape)
            self.camp = self.ca.T
            self.e_camp = np.zeros(self.camp.shape)
            self.pist = self.pistons.T
            self.e_pist = np.zeros(self.pist.shape)

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
        N_ap = len(self.staxy)
        tel_name = ['A%i' % x for x in np.arange(N_ap)+1]
        sta_name = tel_name
        diameter = [0] * N_ap

        staxyz = []
        for x in self.staxy:
            a = list(x)
            line = [a[0], a[1], 0]
            staxyz.append(line)

        sta_index = np.arange(N_ap) + 1

        pscale = info4oif.pscale_mas/1000.  # arcsec
        isz = self.ff.scidata.shape[1]  # Size of the image to extract NRM data
        fov = [pscale * isz] * N_ap
        fovtype = ['RADIUS'] * N_ap

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
                     'T3AMP': self.camp,
                     'T3AMPERR': self.e_camp,
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
                        'PISTONS': self.pist,
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
                    'ISZ': isz,  # size of the image needed (or fov)
                    'NFILE': 0,
                    'ARRAYX': float(0),
                    'ARRAYY': float(0),
                    'ARRAYZ': float(0),
                    'PA': self.pa,
                    'FRAME': 'SKY',
                    'OI_REVN': 2
                    }
           }

    def init_oimodel_arrays(self, oimodel):
        """
        Short Summary
        ------------
        Set dtypes and initialize shapes for AmiOiModel arrays, 
        depending on if averaged or multi-integration version.

        Parameters
        ----------
        oimodel: AmiOIModel object
            empty model

        Returns
        -------
        
        """
        if self.method == 'multi':
            # update dimensions of arrays for multi-integration oifits
            target_dtype = oimodel.target.dtype 
            wavelength_dtype = np.dtype([('EFF_WAVE', '<f4'), 
                                ('EFF_BAND', '<f4')])
            array_dtype = np.dtype([('TEL_NAME', 'S16'), 
                                ('STA_NAME', 'S16'), 
                                ('STA_INDEX', '<i2'), 
                                ('DIAMETER', '<f4'), 
                                ('STAXYZ', '<f8', (3,)), 
                                ('FOV', '<f8'), 
                                ('FOVTYPE', 'S6'), 
                                ('CTRS_EQT', '<f8', (2,)),
                                ('PISTONS', '<f8', (self.nslices,)), 
                                ('PIST_ERR', '<f8', (self.nslices,))])
            vis_dtype = np.dtype([('TARGET_ID', '<i2'), 
                                ('TIME', '<f8'), 
                                ('MJD', '<f8'), 
                                ('INT_TIME', '<f8'), 
                                ('VISAMP', '<f8', (self.nslices,)), 
                                ('VISAMPERR', '<f8', (self.nslices,)), 
                                ('VISPHI', '<f8', (self.nslices,)), 
                                ('VISPHIERR', '<f8', (self.nslices,)), 
                                ('UCOORD', '<f8'), 
                                ('VCOORD', '<f8'), 
                                ('STA_INDEX', '<i2', (2,)), 
                                ('FLAG', 'i1')])
            vis2_dtype = np.dtype([('TARGET_ID', '<i2'), 
                                ('TIME', '<f8'), 
                                ('MJD', '<f8'), 
                                ('INT_TIME', '<f8'), 
                                ('VIS2DATA', '<f8', (self.nslices,)), 
                                ('VIS2ERR', '<f8', (self.nslices,)), 
                                ('UCOORD', '<f8'), 
                                ('VCOORD', '<f8'), 
                                ('STA_INDEX', '<i2', (2,)), 
                                ('FLAG', 'i1')])
            t3_dtype = np.dtype([('TARGET_ID', '<i2'), 
                                ('TIME', '<f8'), 
                                ('MJD', '<f8'), 
                                ('INT_TIME', '<f8'), 
                                ('T3AMP', '<f8', (self.nslices,)), 
                                ('T3AMPERR', '<f8', (self.nslices,)), 
                                ('T3PHI', '<f8', (self.nslices,)), 
                                ('T3PHIERR', '<f8', (self.nslices,)), 
                                ('U1COORD', '<f8'), 
                                ('V1COORD', '<f8'), 
                                ('U2COORD', '<f8'), 
                                ('V2COORD', '<f8'), 
                                ('STA_INDEX', '<i2', (3,)), 
                                ('FLAG', 'i1')])
        else:
            target_dtype = oimodel.target.dtype 
            wavelength_dtype = oimodel.wavelength.dtype
            array_dtype = np.dtype([('TEL_NAME', 'S16'), 
                                ('STA_NAME', 'S16'), 
                                ('STA_INDEX', '<i2'), 
                                ('DIAMETER', '<f4'), 
                                ('STAXYZ', '<f8', (3,)), 
                                ('FOV', '<f8'), 
                                ('FOVTYPE', 'S6'), 
                                ('CTRS_EQT', '<f8', (2,)),
                                ('PISTONS', '<f8'), 
                                ('PIST_ERR', '<f8'),
                                ])
            vis_dtype = oimodel.vis.dtype
            vis2_dtype = oimodel.vis2.dtype
            t3_dtype = oimodel.t3.dtype
        oimodel.array = np.zeros(self.nh, dtype=array_dtype)
        oimodel.target = np.zeros(1, dtype=target_dtype)
        oimodel.vis = np.zeros(self.nbl, dtype=vis_dtype)
        oimodel.vis2 = np.zeros(self.nbl, dtype=vis2_dtype)
        oimodel.t3 = np.zeros(self.ncp, dtype=t3_dtype)
        oimodel.wavelength = np.zeros(1, dtype=wavelength_dtype)


    def populate_oimodel(self):
        """
        Short Summary
        ------------
        Populate the AmiOIModel with the data.

        Parameters
        ----------

        Returns
        -------
        
        """
        m = datamodels.AmiOIModel()
        self.init_oimodel_arrays(m)

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
        m.array['TEL_NAME'] = self.oifits_dct['OI_ARRAY']['TEL_NAME']
        m.array['STA_NAME'] = self.oifits_dct['OI_ARRAY']['STA_NAME']
        m.array['STA_INDEX'] = self.oifits_dct['OI_ARRAY']['STA_INDEX']
        m.array['DIAMETER'] = self.oifits_dct['OI_ARRAY']['DIAMETER']
        m.array['STAXYZ'] = self.oifits_dct['OI_ARRAY']['STAXYZ']
        m.array['FOV'] = self.oifits_dct['OI_ARRAY']['FOV']
        m.array['FOVTYPE'] = self.oifits_dct['OI_ARRAY']['FOVTYPE']
        m.array['CTRS_EQT'] = self.oifits_dct['OI_ARRAY']['CTRS_EQT']
        m.array['PISTONS'] = self.oifits_dct['OI_ARRAY']['PISTONS']
        m.array['PIST_ERR'] = self.oifits_dct['OI_ARRAY']['PIST_ERR']
        

        # oi_target extension data
        m.target['TARGET_ID'] = self.oifits_dct['OI_TARGET']['TARGET_ID']
        m.target['TARGET'] = self.oifits_dct['OI_TARGET']['TARGET']
        m.target['RAEP0'] = self.oifits_dct['OI_TARGET']['RAEP0']
        m.target['DECEP0'] = self.oifits_dct['OI_TARGET']['DECEP0']
        m.target['EQUINOX'] = self.oifits_dct['OI_TARGET']['EQUINOX']
        m.target['RA_ERR'] = self.oifits_dct['OI_TARGET']['RA_ERR']
        m.target['DEC_ERR'] = self.oifits_dct['OI_TARGET']['DEC_ERR']
        m.target['SYSVEL'] = self.oifits_dct['OI_TARGET']['SYSVEL']
        m.target['VELTYP'] = self.oifits_dct['OI_TARGET']['VELTYP']
        m.target['VELDEF'] = self.oifits_dct['OI_TARGET']['VELDEF']
        m.target['PMRA'] = self.oifits_dct['OI_TARGET']['PMRA']
        m.target['PMDEC'] = self.oifits_dct['OI_TARGET']['PMDEC']
        m.target['PMRA_ERR'] = self.oifits_dct['OI_TARGET']['PMRA_ERR']
        m.target['PMDEC_ERR'] = self.oifits_dct['OI_TARGET']['PMDEC_ERR']
        m.target['PARALLAX'] = self.oifits_dct['OI_TARGET']['PARALLAX']
        m.target['PARA_ERR'] = self.oifits_dct['OI_TARGET']['PARA_ERR']
        m.target['SPECTYP'] = self.oifits_dct['OI_TARGET']['SPECTYP']

        # oi_vis extension data
        m.vis['TARGET_ID'] = self.oifits_dct['OI_VIS']['TARGET_ID']
        m.vis['TIME'] = self.oifits_dct['OI_VIS']['TIME']
        m.vis['MJD'] = self.oifits_dct['OI_VIS']['MJD']
        m.vis['INT_TIME'] = self.oifits_dct['OI_VIS']['INT_TIME']
        m.vis['VISAMP'] = self.oifits_dct['OI_VIS']['VISAMP']
        m.vis['VISAMPERR'] = self.oifits_dct['OI_VIS']['VISAMPERR']
        m.vis['VISPHI'] = self.oifits_dct['OI_VIS']['VISPHI']
        m.vis['VISPHIERR'] = self.oifits_dct['OI_VIS']['VISPHIERR']
        m.vis['UCOORD'] = self.oifits_dct['OI_VIS']['UCOORD']
        m.vis['VCOORD'] = self.oifits_dct['OI_VIS']['VCOORD']
        m.vis['STA_INDEX'] = self.oifits_dct['OI_VIS']['STA_INDEX']
        m.vis['FLAG'] = self.oifits_dct['OI_VIS']['FLAG']

        # oi_vis2 extension data
        m.vis2['TARGET_ID'] = self.oifits_dct['OI_VIS2']['TARGET_ID']
        m.vis2['TIME'] = self.oifits_dct['OI_VIS2']['TIME']
        m.vis2['MJD'] = self.oifits_dct['OI_VIS2']['MJD']
        m.vis2['INT_TIME'] = self.oifits_dct['OI_VIS2']['INT_TIME']
        m.vis2['VIS2DATA'] = self.oifits_dct['OI_VIS2']['VIS2DATA']
        m.vis2['VIS2ERR'] = self.oifits_dct['OI_VIS2']['VIS2ERR']
        m.vis2['UCOORD'] = self.oifits_dct['OI_VIS2']['UCOORD']
        m.vis2['VCOORD'] = self.oifits_dct['OI_VIS2']['VCOORD']
        m.vis2['STA_INDEX'] = self.oifits_dct['OI_VIS2']['STA_INDEX']
        m.vis2['FLAG'] = self.oifits_dct['OI_VIS2']['FLAG']

        # oi_t3 extension data
        m.t3['TARGET_ID'] = self.oifits_dct['OI_T3']['TARGET_ID']
        m.t3['TIME'] = self.oifits_dct['OI_T3']['TIME']
        m.t3['MJD'] = self.oifits_dct['OI_T3']['MJD']
        m.t3['T3AMP'] = self.oifits_dct['OI_T3']['T3AMP']
        m.t3['T3AMPERR'] = self.oifits_dct['OI_T3']['T3AMPERR']
        m.t3['T3PHI'] = self.oifits_dct['OI_T3']['T3PHI']
        m.t3['T3PHIERR'] = self.oifits_dct['OI_T3']['T3PHIERR']
        m.t3['U1COORD'] = self.oifits_dct['OI_T3']['U1COORD']
        m.t3['V1COORD'] = self.oifits_dct['OI_T3']['V1COORD']
        m.t3['U2COORD'] = self.oifits_dct['OI_T3']['U2COORD']
        m.t3['V2COORD'] = self.oifits_dct['OI_T3']['V2COORD']
        m.t3['STA_INDEX'] = self.oifits_dct['OI_T3']['STA_INDEX']
        m.t3['FLAG'] = self.oifits_dct['OI_T3']['FLAG']

        # oi_wavelength extension data
        m.wavelength['EFF_WAVE'] = self.oifits_dct['OI_WAVELENGTH']['EFF_WAVE']
        m.wavelength['EFF_BAND'] = self.oifits_dct['OI_WAVELENGTH']['EFF_BAND']

        return m

    def _maketriples_all(self):
        """
        Short Summary
        ------------
        Calculate all three-hole combinations, baselines

        Parameters
        ----------

        Returns
        -------
        tarray: integer array
            Triple hole indices (0-indexed), 
        float array of two uv vectors in all triangles
        """
        nholes = self.ctrs_eqt.shape[0]
        tlist = []
        for i in range(nholes):
            for j in range(nholes):
                for k in range(nholes):
                    if i < j and j < k:
                        tlist.append((i, j, k))
        tarray = np.array(tlist).astype(int)
        

        tname = []
        uvlist = []
        # foreach row of 3 elts...
        for triple in tarray:
            tname.append("{0:d}_{1:d}_{2:d}".format(
                triple[0], triple[1], triple[2]))

            uvlist.append((self.ctrs_eqt[triple[0]] - self.ctrs_eqt[triple[1]],
                           self.ctrs_eqt[triple[1]] - self.ctrs_eqt[triple[2]]))
        # print(len(uvlist), "uvlist", uvlist)

        return tarray, np.array(uvlist)

    def _makebaselines(self):
        """
        Short Summary
        ------------
        Calculate all hole pairs, baselines

        Parameters
        ----------

        Returns
        -------
        barray: list
            Hole pairs indices, 0-indexed
        float array of baselines
        """
        nholes = self.ctrs_eqt.shape[0]
        blist = []
        for i in range(nholes):
            for j in range(nholes):
                if i < j:
                    blist.append((i, j))
        barray = np.array(blist).astype(int)
        # blname = []
        bllist = []
        for basepair in blist:
            # blname.append("{0:d}_{1:d}".format(basepair[0],basepair[1]))
            baseline = self.ctrs_eqt[basepair[0]] - self.ctrs_eqt[basepair[1]]
            bllist.append(baseline)
        return barray, np.array(bllist)



class CalibOifits:
    def __init__(self,targoimodel,caloimodel):
        """
        Short Summary
        -------------
        Calibrate (normalize) an AMI observation by subtracting closure phases
        of a reference star from those of a target and dividing visibility amplitudes
        of the target by those of the reference star.

        Parameters
        ----------
        targoimodel: AmiOIModlel, target
        caloimodel: AmiOIModlel, reference star (calibrator)

        """
        self.targoimodel = targoimodel
        self.caloimodel = caloimodel
        self.calib_oimodel = targoimodel.copy()

    def update_dtype(self):
        """
        Short Summary
        ------------
        Modify the dtype of OI array to include different pistons columns
        for calibrated OIFITS files
    
        Parameters
        ----------
    
        Returns
        -------
        """
        nrows=7
        modified_dtype = np.dtype([('TEL_NAME', 'S16'), 
                                ('STA_NAME', 'S16'), 
                                ('STA_INDEX', '<i2'), 
                                ('DIAMETER', '<f4'), 
                                ('STAXYZ', '<f8', (3,)), 
                                ('FOV', '<f8'), 
                                ('FOVTYPE', 'S6'), 
                                ('CTRS_EQT', '<f8', (2,)),
                                ('PISTON_T', '<f8'),
                                ('PISTON_C', '<f8'),
                                ('PIST_ERR', '<f8'),
                                ])
        self.calib_oimodel.array = np.zeros(nrows, dtype=modified_dtype)


    def calibrate(self):
        """
        Short Summary
        ------------
        Apply the calibration (normalization) routine to calibrate the 
        target AmiOIModel by the calibrator (reference star) AmiOIModel
    
        Parameters
        ----------
    
        Returns
        -------
        calib_oimodel: AmiOIModel
            Calibrated AMI datamodel
        """
        cp_out = self.targoimodel.t3['T3PHI'] - self.caloimodel.t3['T3PHI']
        sqv_out = self.targoimodel.vis2['VIS2DATA'] / self.caloimodel.vis2['VIS2DATA']
        va_out = self.targoimodel.vis['VISAMP'] / self.caloimodel.vis['VISAMP']
        # using standard propagation of error for multiplication/division
        # which assumes uncorrelated Gaussian errors (questionable)    
        cperr_t = self.targoimodel.t3 ['T3PHIERR']
        cperr_c = self.caloimodel.t3['T3PHIERR']
        sqverr_c = self.targoimodel.vis2['VIS2ERR']
        sqverr_t = self.caloimodel.vis2['VIS2ERR']
        vaerr_t = self.targoimodel.vis['VISAMPERR']
        vaerr_c = self.caloimodel.vis['VISAMPERR']
        cperr_out = np.sqrt(cperr_t**2. + cperr_c**2.)
        sqverr_out = sqv_out * np.sqrt((sqverr_t/self.targoimodel.vis2['VIS2DATA'])**2. + 
                                        (sqverr_c/self.caloimodel.vis2['VIS2DATA'])**2.)
        vaerr_out = va_out * np.sqrt((vaerr_t/self.targoimodel.vis['VISAMP'])**2. + 
                                    (vaerr_c/self.caloimodel.vis['VISAMP'])**2.)

        pistons_t = self.targoimodel.array['PISTONS']
        pisterr_t = self.targoimodel.array['PIST_ERR']
        pistons_c = self.caloimodel.array['PISTONS']
        pisterr_c = self.caloimodel.array['PIST_ERR']
        # sum in quadrature errors from target and calibrator pistons
        pisterr_out = np.sqrt(pisterr_t**2 + pisterr_c**2)

        # update OI array, which is currently all zeros, with input oi array 
        # and updated piston columns
        self.update_dtype()
        self.calib_oimodel.array['TEL_NAME'] = self.targoimodel.array['TEL_NAME']
        self.calib_oimodel.array['STA_NAME'] = self.targoimodel.array['STA_NAME']
        self.calib_oimodel.array['STA_INDEX'] = self.targoimodel.array['STA_INDEX']
        self.calib_oimodel.array['DIAMETER'] = self.targoimodel.array['DIAMETER']
        self.calib_oimodel.array['STAXYZ'] = self.targoimodel.array['STAXYZ']
        self.calib_oimodel.array['FOV'] = self.targoimodel.array['FOV']
        self.calib_oimodel.array['FOVTYPE'] = self.targoimodel.array['FOVTYPE']
        self.calib_oimodel.array['CTRS_EQT'] = self.targoimodel.array['CTRS_EQT']
        self.calib_oimodel.array['PISTON_T'] = pistons_t
        self.calib_oimodel.array['PISTON_C'] = pistons_c
        self.calib_oimodel.array['PIST_ERR'] = pisterr_out

        # update calibrated oimodel arrays with calibrated observables
        self.calib_oimodel.t3['T3PHI'] = cp_out
        self.calib_oimodel.t3['T3PHIERR'] = cperr_out
        self.calib_oimodel.vis2['VIS2DATA'] = sqv_out
        self.calib_oimodel.vis2['VIS2ERR'] = sqverr_out
        self.calib_oimodel.vis['VISAMP'] = va_out
        self.calib_oimodel.vis['VISAMPERR'] = vaerr_out

        self.calib_oimodel.array['PISTON_T'] = pistons_t
        self.calib_oimodel.array['PISTON_C'] = pistons_c
        self.calib_oimodel.array['PIST_ERR'] = pisterr_out

        # add calibrated header keywords
        calname = self.caloimodel.meta.target.proposer_name # name of calibrator star
        self.calib_oimodel.meta.oifits.calib = calname

        return self.calib_oimodel




