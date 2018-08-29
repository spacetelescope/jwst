"""
Test functions for NIRSPEC WCS - all modes.
"""
import os.path
import numpy as np
from numpy.testing.utils import assert_allclose
from astropy.io import fits
from astropy.modeling import models as astmodels
from astropy import wcs as astwcs
from gwcs import wcs
from ... import datamodels
from ...transforms.models import Slit
from .. import nirspec
from .. import assign_wcs_step
from . import data
from ..util import MissingMSAFileError

import pytest

data_path = os.path.split(os.path.abspath(data.__file__))[0]


wcs_kw = {'wcsaxes': 2, 'ra_ref': 165, 'dec_ref': 54,
          'v2_ref': -8.3942412, 'v3_ref': -5.3123744, 'roll_ref': 37,
          'crpix1': 1024, 'crpix2': 1024,
          'cdelt1': .08, 'cdelt2': .08,
          'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
          'pc1_1': 1, 'pc1_2': 0, 'pc2_1': 0, 'pc2_2': 1
          }


slit_fields_num = ["shutter_id", "xcen", "ycen",
                   "ymin", "ymax", "quadrant", "source_id",
                   "stellarity", "source_xpos", "source_ypos"]


slit_fields_str = ["name", "shutter_state", "source_name", "source_alias"]


def _compare_slits(s1, s2):
    for f in slit_fields_num:
        assert_allclose(getattr(s1, f), getattr(s2, f))
    for f in slit_fields_str:
        assert getattr(s1, f) == getattr(s2, f)


def get_file_path(filename):
    """
    Construct an absolute path.
    """
    return os.path.join(data_path, filename)


def create_hdul(detector='NRS1'):
    """
    Create a fits HDUList instance.
    """
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['instrume'] = 'NIRSPEC'
    phdu.header['detector'] = detector
    phdu.header['time-obs'] = '8:59:37'
    phdu.header['date-obs'] = '2016-09-05'

    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    for item in wcs_kw.items():
        scihdu.header[item[0]] = item[1]
    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def create_reference_files(datamodel):
    """
    Create a dict {reftype: reference_file}.
    """
    refs = {}
    step = assign_wcs_step.AssignWcsStep()
    for reftype in assign_wcs_step.AssignWcsStep.reference_file_types:
        refs[reftype] = step.get_reference_file(datamodel, reftype)
    return refs


def create_nirspec_imaging_file():
    image = create_hdul()
    image[0].header['exp_type'] = 'NRS_IMAGE'
    image[0].header['filter'] = 'F290LP'
    image[0].header['grating'] = 'MIRROR'
    return image


def create_nirspec_mos_file():
    image = create_hdul()
    image[0].header['exp_type'] = 'NRS_MSASPEC'
    image[0].header['filter'] = 'F170LP'
    image[0].header['grating'] = 'G235M'
    image[0].header['crval3'] = 0
    image[0].header['wcsaxes'] = 3
    image[0].header['ctype3'] = 'WAVE'
    image[0].header['pc3_1'] = 1
    image[0].header['pc3_2'] = 0

    msa_status_file = get_file_path('SPCB-GD-A.msa.fits.gz')
    image[0].header['MSACONFG'] = msa_status_file
    return image


def create_nirspec_ifu_file(filter, grating, lamp='N/A', detector='NRS1'):
    image = create_hdul(detector)
    image[0].header['exp_type'] = 'NRS_IFU'
    image[0].header['filter'] = filter
    image[0].header['grating'] = grating
    image[1].header['crval3'] = 0
    image[1].header['wcsaxes'] = 3
    image[1].header['ctype3'] = 'WAVE'
    image[0].header['lamp'] = lamp
    image[0].header['GWA_XTIL'] = 0.3318742513656616
    image[0].header['GWA_YTIL'] = 0.1258982867002487
    return image


def create_nirspec_fs_file(grating, filter, lamp="N/A"):
    image = create_hdul()
    image[0].header['exp_type'] = 'NRS_FIXEDSLIT'
    image[0].header['filter'] = filter
    image[0].header['grating'] = grating
    image[0].header['lamp'] = lamp
    image[1].header['crval3'] = 0
    image[1].header['wcsaxes'] = 3
    image[1].header['ctype3'] = 'WAVE'
    image[0].header['GWA_XTIL'] = 0.3316612243652344
    image[0].header['GWA_YTIL'] = 0.1260581910610199
    image[0].header['SUBARRAY'] = "FULL"
    return image


def test_nirspec_imaging():
    """
    Test Nirspec Imaging mode using build 6 reference files.
    """
    #Test creating the WCS
    f = create_nirspec_imaging_file()
    im = datamodels.ImageModel(f)

    refs = create_reference_files(im)

    pipe = nirspec.create_pipeline(im, refs)
    w = wcs.WCS(pipe)
    im.meta.wcs = w
    # Test evaluating the WCS
    im.meta.wcs(1, 2)


def test_nirspec_ifu_against_esa():
    """
    Test Nirspec IFU mode using CV3 reference files.
    """
    ref = fits.open(get_file_path('Trace_IFU_Slice_00_SMOS-MOD-G1M-17-5344175105_30192_JLAB88.fits'))

    # Test NRS1
    pyw = astwcs.WCS(ref['SLITY1'].header)
    hdul = create_nirspec_ifu_file("OPAQUE", "G140M")
    im = datamodels.ImageModel(hdul)
    im.meta.filename = "test_ifu.fits"
    refs = create_reference_files(im)

    pipe = nirspec.create_pipeline(im, refs)
    w = wcs.WCS(pipe)
    im.meta.wcs = w
    # Test evaluating the WCS (slice 0)
    w0 = nirspec.nrs_wcs_set_input(im, 0)

    # get positions within the slit and the coresponding lambda
    slit1 = ref['SLITY1'].data # y offset on the slit
    lam = ref['LAMBDA1'].data
    # filter out locations outside the slit
    cond = np.logical_and(slit1 < .5, slit1 > -.5)
    y, x = cond.nonzero() # 0-based

    x, y = pyw.wcs_pix2world(x, y, 0)
    # The pipeline accepts 0-based cooridnates
    x -= 1
    y -= 1
    sca2world = w0.get_transform('sca', 'msa_frame')
    _, slit_y, lp = sca2world(x, y)

    lp *= 10**-6
    assert_allclose(lp, lam[cond], atol=1e-13)


def test_nirspec_fs_esa():
    """
    Test Nirspec FS mode using build 6 reference files.
    """
    #Test creating the WCS
    filename = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    im = datamodels.ImageModel(filename)
    im.meta.filename = "test_fs.fits"
    refs = create_reference_files(im)

    pipe = nirspec.create_pipeline(im, refs)
    w = wcs.WCS(pipe)
    im.meta.wcs = w
    # Test evaluating the WCS
    w1 = nirspec.nrs_wcs_set_input(im, "S200A1")

    ref = fits.open(get_file_path('Trace_SLIT_A_200_1_V84600010001P0000000002101_39547_JLAB88.fits'))
    pyw = astwcs.WCS(ref[1].header)

    # get positions within the slit and the coresponding lambda
    slit1 = ref[5].data # y offset on the slit
    lam = ref[4].data

    # filter out locations outside the slit
    cond = np.logical_and(slit1 < .5, slit1 > -.5)
    y, x = cond.nonzero() # 0-based

    x, y = pyw.wcs_pix2world(x, y, 0)
    # The pipeline works with 0-based coordinates
    x -= 1
    y -= 1

    sca2world = w1.get_transform('sca', 'v2v3')
    ra, dec, lp = sca2world(x, y)
    # w1 now outputs in microns hence the 1e6 factor
    lp *= 1e-6
    lam = lam[cond]
    nan_cond = ~np.isnan(lp)
    assert_allclose(lp[nan_cond], lam[nan_cond], atol=10**-13)
    ref.close()


def test_correct_tilt():
    """
    Example provided by Catarina.
    """
    disp = datamodels.DisperserModel()
    xtilt = 0.35896975
    ytilt = 0.1343827
    # ztilt = None
    corrected_theta_x = 0.02942671219861111
    corrected_theta_y = 0.00018649006677464447
    # corrected_theta_z = -0.2523269848788889
    disp.gwa_tiltx =  {'temperatures': [39.58],
                       'tilt_model': astmodels.Polynomial1D(1, c0=3307.85402614,
                                                            c1=-9182.87552123),
                       'unit': 'arcsec',
                       'zeroreadings': [0.35972327]}
    disp.gwa_tilty = {'temperatures': [39.58],
                      'tilt_model': astmodels.Polynomial1D(1, c0=0.0, c1=0.0),
                      'unit': 'arcsec',
                      'zeroreadings': [0.0]}
    disp.meta = {'instrument': {'name': 'NIRSPEC', 'detector': 'NRS1'},
                 'reftype': 'DISPERSER'}

    disp.theta_x = 0.02942671219861111
    disp.theta_y = -0.0007745488724972222
    # disp.theta_z = -0.2523269848788889
    disp.tilt_x = 0.0
    disp.tilt_y = -8.8

    disp_corrected = nirspec.correct_tilt(disp, xtilt, ytilt)#, ztilt)
    assert np.isclose(disp_corrected.theta_x, corrected_theta_x)
    # assert(np.isclose(disp_corrected['theta_z'], corrected_theta_z))
    assert np.isclose(disp_corrected.theta_y, corrected_theta_y)


def test_msa_configuration_normal():
    """
    Test the get_open_msa_slits function.
    """

    # Test 1: Reasonably normal as well
    msa_meta_id = 12
    msaconfl = get_file_path('msa_configuration.fits')
    slitlet_info = nirspec.get_open_msa_slits(msaconfl, msa_meta_id)
    ref_slit = Slit(55, 9376, 251, 26, -5.15, 0.55, 4, 1, '1111x', '95065_1', '2122',
                      0.13, -0.31716078999999997, 0.18092266)
    _compare_slits(slitlet_info[0], ref_slit)


def test_msa_configuration_no_background():
    """
    Test the get_open_msa_slits function.
    """
    # Test 2: Two main shutters, not allowed and should fail
    msa_meta_id = 13
    msaconfl = get_file_path('msa_configuration.fits')
    with pytest.raises(ValueError):
        nirspec.get_open_msa_slits(msaconfl, msa_meta_id)


def test_msa_configuration_all_background():
    """
    Test the get_open_msa_slits function.
    """

    # Test 3:  No non-background, not acceptable.
    msa_meta_id = 14
    msaconfl = get_file_path('msa_configuration.fits')
    slitlet_info = nirspec.get_open_msa_slits(msaconfl, msa_meta_id)
    ref_slit = Slit(57, 8646, 251, 24, -2.85, .55, 4, 1, '11x', '95065_1', '2122',
                    0.13, -0.5, 0.5)
    _compare_slits(slitlet_info[0], ref_slit)



def test_msa_configuration_row_skipped():
    """
    Test the get_open_msa_slits function.
    """

    # Test 4: One row is skipped, should be acceptable.
    msa_meta_id = 15
    msaconfl = get_file_path('msa_configuration.fits')
    slitlet_info = nirspec.get_open_msa_slits(msaconfl, msa_meta_id)
    ref_slit = Slit(58, 8646, 251, 24, -2.85, 5.15, 4, 1, '11x1011', '95065_1', '2122',
                      0.130, -0.31716078999999997, 0.18092266)
    _compare_slits(slitlet_info[0], ref_slit)


def test_msa_configuration_multiple_returns():
    """
    Test the get_open_msa_slits function.
    """
    # Test 4: One row is skipped, should be acceptable.
    msa_meta_id = 16
    msaconfl = get_file_path('msa_configuration.fits')
    slitlet_info = nirspec.get_open_msa_slits(msaconfl, msa_meta_id)
    ref_slit1 = Slit(59, 8651, 256, 24, -2.85, 5.15, 4, 1, '11x1011', '95065_1', '2122',
                     0.13000000000000003, -0.31716078999999997, 0.18092266)
    ref_slit2 = Slit(60, 11573, 258, 32, -2.85, 4, 4, 2, '11x111', '95065_2', '172',
                     0.70000000000000007, -0.31716078999999997, 0.18092266)
    _compare_slits(slitlet_info[0], ref_slit1)
    _compare_slits(slitlet_info[1], ref_slit2)


open_shutters = [[24], [23, 24], [22, 23, 25, 27], [22, 23, 25, 27, 28]]
main_shutter = [24, 23, 25, 28]
result = ["x", "x1", "110x01", "110101x"]
test_data = list(zip(open_shutters, main_shutter, result))

@pytest.mark.parametrize(('open_shutters', 'main_shutter', 'result'),
                         test_data)
def test_shutter_state(open_shutters, main_shutter, result):
    shutter_state = nirspec._shutter_id_to_str(open_shutters, main_shutter)
    assert shutter_state == result


def test_slit_projection_on_detector():
    step = assign_wcs_step.AssignWcsStep()

    hdul = create_nirspec_fs_file(grating="G395M", filter="OPAQUE", lamp="ARGON")
    hdul[0].header['DETECTOR'] = 'NRS2'
    im = datamodels.ImageModel(hdul)

    refs = {}
    for reftype in step.reference_file_types:
        refs[reftype] = step.get_reference_file(im, reftype)

    open_slits = nirspec.get_open_slits(im, refs)
    assert len(open_slits) == 1
    assert open_slits[0].name == "S200B1"

    hdul[0].header['DETECTOR'] = 'NRS1'
    im = datamodels.ImageModel(hdul)

    refs = {}
    for reftype in step.reference_file_types:
        refs[reftype] = step.get_reference_file(im, reftype)

    open_slits = nirspec.get_open_slits(im, refs)
    assert len(open_slits) == 4
    names = [s.name for s in open_slits]
    assert "S200A1" in names
    assert "S200A2" in names
    assert "S400A1" in names
    assert "S1600A1" in names


def test_missing_msa_file():
    image = create_nirspec_mos_file()
    model = datamodels.ImageModel(image)

    model.meta.instrument.msa_metadata_file = ""
    with pytest.raises(MissingMSAFileError):
        assign_wcs_step.AssignWcsStep.call(model)

    model.meta.instrument.msa_metadata_file = "missing.fits"
    with pytest.raises(MissingMSAFileError):
        assign_wcs_step.AssignWcsStep.call(model)


def test_open_slits():
    """ Test that get_open_slits works with MSA data.

    Issue #2321
    """
    image = create_nirspec_mos_file()
    model = datamodels.ImageModel(image)
    msaconfl = get_file_path('msa_configuration.fits')

    model.meta.instrument.msa_metadata_file = msaconfl
    model.meta.instrument.msa_metadata_id=12

    slits = nirspec.get_open_slits(model)
    assert len(slits) == 1
