from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from astropy.io import fits
from astropy.modeling import models as astmodels
from astropy.utils.data import get_pkg_data_filename

from ...datamodels import ImageModel
from ..import nirspec
from .. import assign_wcs_step

from asdf import AsdfFile
import numpy as np
from numpy.testing.utils import assert_allclose


def test_correct_tilt():
    """
    Example provided by Catarina.
    """
    xtilt = 0.35896975
    ytilt = 0.1343827
    ztilt = None
    corrected_theta_x = 0.02942671219861111
    corrected_theta_y = 0.00018649006677464447
    corrected_theta_z = -0.2523269848788889
    disp = {'gwa_tiltx': {'temperatures': [39.58],
                          'tilt_model': astmodels.Polynomial1D(1, c0=3307.85402614, c1=-9182.87552123),
                          'unit': 'arcsec',
                          'zeroreadings': [0.35972327]},
            'gwa_tilty': {'temperatures': [39.58],
                          'tilt_model': astmodels.Polynomial1D(1, c0=0.0, c1=0.0),
                          'unit': 'arcsec',
                          'zeroreadings': [0.0]},
            'instrument': 'NIRSPEC',
            'reftype': 'DISPERSER',
            'theta_x': 0.02942671219861111,
            'theta_y': -0.0007745488724972222,
            'theta_z': -0.2523269848788889,
            'tilt_x': 0.0,
            'tilt_y': -8.8
            }
    disp_corrected = nirspec.correct_tilt(disp, xtilt, ytilt, ztilt)
    assert(np.isclose(disp_corrected['theta_x'], corrected_theta_x))
    assert(np.isclose(disp_corrected['theta_z'], corrected_theta_z))
    assert(np.isclose(disp_corrected['theta_y'], corrected_theta_y))

"""
Test against the ESA regression data.
Configuration is:
detector: NRS1, NRS2
filter: CLEAR
grating: MIRROR

All tilt angles are 0.
"""

class TestMirror(object):
    def setup_class(self):
        self.nrs = ImageModel()
        self.nrs.meta.instrument.name = "NIRSPEC"
        self.nrs.meta.instrument.filter = "CLEAR"
        self.nrs.meta.instrument.grating = "MIRROR"
        self.nrs.meta.exposure.type = 'NRS_FIXEDSLIT'
        self.nrs.meta.instrument.detector = 'NRS1'
        self.slits = {'S200A1': 'SLIT_A_200_1',
                      'S200A2': 'SLIT_A_200_2',
                      'S400A1': 'SLIT_A_400',
                      'S1600A1': 'SLIT_A_1600',
                      #'S200B1': 'SLIT_B_200',
                      #'IFU': 'IFU_window'
                      }
        self.reference = fits.open(get_pkg_data_filename(
            "data/onSkyAndOnDetectorProjectionSLIT.fits.gz"))
        self.field_names = list(self.reference[1].data.field('name'))
        print(self.field_names)
        step = assign_wcs_step.AssignWcsStep()
        self.model = step.process(self.nrs)
        self.lam = 2e-6
        self.wcs_S200A1 = self.model.meta.wcs_S200A1
        self.wcs_S200A2 = self.model.meta.wcs_S200A2
        self.wcs_S400A1 = self.model.meta.wcs_S400A1
        self.wcs_S1600A1 = self.model.meta.wcs_S1600A1
        #self.wcs_S200B1 = self.model.meta.wcs_S200B1

    def _get_reference_data(self, slit):
        #names = list(self.reference[1].data.field('name'))
        slit_index = self.field_names.index(self.slits[slit])
        return self.reference[1].data[slit_index]

    def test_msa(self):
        for slit in self.slits:
            msa = AsdfFile.open(self.model.meta.ref_file.msa.name).tree[slit]['model']
            name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
                xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data(slit)
            assert_allclose(msa(0, 0), (xmsa, ymsa))

    def test_v23_s200a1(self):
        name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
            xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data('S200A1')
        msa2v23 = self.wcs_S200A1.get_transform('msa', 'v2v3')
        v2, v3, lam = msa2v23(0, 0, self.lam)
        assert_allclose((v2, v3), (xsky, ysky))

    def test_v23_s200a2(self):
        name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
            xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data('S200A2')
        msa2v23 = self.wcs_S200A2.get_transform('msa', 'v2v3')
        v2, v3, lam = msa2v23(0, 0, self.lam)
        assert_allclose((v2, v3), (xsky, ysky))

    def test_v23_s400a1(self):
        name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
            xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data('S400A1')
        msa2v23 = self.wcs_S400A1.get_transform('msa', 'v2v3')
        v2, v3, lam = msa2v23(0, 0, self.lam)
        assert_allclose((v2, v3), (xsky, ysky))

    def test_v23_s1600a1(self):
        name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
            xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data('S1600A1')
        msa2v23 = self.wcs_S1600A1.get_transform('msa', 'v2v3')
        v2, v3, lam = msa2v23(0, 0, self.lam)
        assert_allclose((v2, v3), (xsky, ysky))

    # The FPA coordinates start have origin (1,1)
    def test_sca1_s200a1(self):
        name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
            xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data('S200A1')
        msa2det = self.wcs_S200A1.get_transform('msa', 'detector')
        x, y = msa2det(0, 0, self.lam) + np.array([1, 1])
        assert_allclose((x, y), (xsca1, ysca1))

    def test_sca1_s200a2(self):
        name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
            xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data('S200A2')
        msa2det = self.wcs_S200A2.get_transform('msa', 'detector')
        x, y = msa2det(0, 0, self.lam) + np.array([1, 1])
        assert_allclose((x, y), (xsca1, ysca1))

    def test_sca1_s400a1(self):
        name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
            xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data('S400A1')
        msa2det = self.wcs_S400A1.get_transform('msa', 'detector')
        x, y = msa2det(0, 0, self.lam) + np.array([1, 1])
        assert_allclose((x, y), (xsca1, ysca1))

    def test_sca1_s1600a1(self):
        name, xmsa, ymsa, xsizemsa, ysizemsa, xsky, ysky, xsizesky, yskysize, \
            xsca1, ysca1, xsca2, ysca2, xsizesca, ysizesca = self._get_reference_data('S1600A1')
        msa2det = self.wcs_S1600A1.get_transform('msa', 'detector')
        x, y = msa2det(0, 0, self.lam) + np.array([1, 1])
        assert_allclose((x, y), (xsca1, ysca1))

