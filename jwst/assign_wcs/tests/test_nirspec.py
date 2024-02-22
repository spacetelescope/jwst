"""
Test functions for NIRSPEC WCS - all modes.
"""
import functools
from math import cos, sin
import os.path

import pytest
import numpy as np
from numpy.testing import assert_allclose
from astropy.io import fits
from astropy.modeling import models as astmodels
from astropy import table
from astropy import wcs as astwcs
import astropy.units as u
import astropy.coordinates as coords
from gwcs import wcs
from gwcs import wcstools

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms import models as trmodels

from jwst.assign_wcs import nirspec
from jwst.assign_wcs import assign_wcs_step
from jwst.assign_wcs.tests import data
from jwst.assign_wcs.util import MSAFileError, in_ifu_slice


data_path = os.path.split(os.path.abspath(data.__file__))[0]


wcs_kw = {'wcsaxes': 2, 'ra_ref': 165, 'dec_ref': 54,
          'v2_ref': -8.3942412, 'v3_ref': -5.3123744, 'roll_ref': 37,
          'crpix1': 1024, 'crpix2': 1024,
          'cdelt1': .08, 'cdelt2': .08,
          'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
          'pc1_1': 1, 'pc1_2': 0, 'pc2_1': 0, 'pc2_2': 1
          }


slit_fields_num = ["shutter_id", "dither_position", "xcen", "ycen",
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


def create_nirspec_mos_file(grating='G235M', filt='F170LP'):
    image = create_hdul()
    image[0].header['exp_type'] = 'NRS_MSASPEC'
    image[0].header['filter'] = filt
    image[0].header['grating'] = grating
    image[0].header['PATT_NUM'] = 1

    msa_status_file = get_file_path('SPCB-GD-A.msa.fits.gz')
    image[0].header['MSAMETFL'] = msa_status_file
    return image


def create_nirspec_ifu_file(filter, grating, lamp='N/A', detector='NRS1',
                            gwa_xtil=0.3318742513656616,
                            gwa_ytil=0.1258982867002487,
                            gwa_tilt=None):
    image = create_hdul(detector)
    image[0].header['exp_type'] = 'NRS_IFU'
    image[0].header['filter'] = filter
    image[0].header['grating'] = grating
    image[1].header['crval3'] = 0
    image[1].header['wcsaxes'] = 3
    image[1].header['ctype3'] = 'WAVE'
    image[0].header['lamp'] = lamp
    image[0].header['GWA_XTIL'] = gwa_xtil
    image[0].header['GWA_YTIL'] = gwa_ytil
    if gwa_tilt is not None:
        image[0].header['GWA_TILT'] = gwa_tilt
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
    image[0].header['FXD_SLIT'] = "S200A1"
    return image


def test_nirspec_imaging():
    """
    Test Nirspec Imaging mode using build 6 reference files.
    """
    # Test creating the WCS
    f = create_nirspec_imaging_file()
    im = datamodels.ImageModel(f)

    refs = create_reference_files(im)

    pipe = nirspec.create_pipeline(im, refs, slit_y_range=[-.5, .5])
    w = wcs.WCS(pipe)
    im.meta.wcs = w
    # Test evaluating the WCS
    im.meta.wcs(1, 2)


def test_nirspec_ifu_against_esa(wcs_ifu_grating):
    """
    Test Nirspec IFU mode using CV3 reference files.
    """
    with fits.open(get_file_path('Trace_IFU_Slice_00_SMOS-MOD-G1M-17-5344175105_30192_JLAB88.fits')) as ref:
        # Test NRS1
        pyw = astwcs.WCS(ref['SLITY1'].header)
        # Test evaluating the WCS (slice 0)
        im, refs = wcs_ifu_grating("G140M", "OPAQUE")
        w0 = nirspec.nrs_wcs_set_input(im, 0)

        # get positions within the slit and the corresponding lambda
        slit1 = ref['SLITY1'].data  # y offset on the slit
        lam = ref['LAMBDA1'].data

    # filter out locations outside the slit
    cond = np.logical_and(slit1 < .5, slit1 > -.5)
    y, x = cond.nonzero()  # 0-based

    x, y = pyw.wcs_pix2world(x, y, 0)
    # The pipeline accepts 0-based coordinates
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
    # Test creating the WCS
    filename = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    im = datamodels.ImageModel(filename)
    im.meta.filename = "test_fs.fits"
    refs = create_reference_files(im)

    pipe = nirspec.create_pipeline(im, refs, slit_y_range=[-.5, .5])
    w = wcs.WCS(pipe)
    im.meta.wcs = w
    # Test evaluating the WCS
    w1 = nirspec.nrs_wcs_set_input(im, "S200A1")

    ref = fits.open(get_file_path('Trace_SLIT_A_200_1_V84600010001P0000000002101_39547_JLAB88.fits'))
    pyw = astwcs.WCS(ref[1].header)

    # get positions within the slit and the corresponding lambda
    slit1 = ref[5].data  # y offset on the slit
    lam = ref[4].data

    # filter out locations outside the slit
    cond = np.logical_and(slit1 < .5, slit1 > -.5)
    y, x = cond.nonzero()  # 0-based

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
    disp.gwa_tiltx = {'temperatures': [39.58],
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

    disp_corrected = nirspec.correct_tilt(disp, xtilt, ytilt)  # , ztilt)
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
    dither_position = 1
    slitlet_info = nirspec.get_open_msa_slits(msaconfl, msa_meta_id, dither_position,
                                              slit_y_range=[-.5, .5])
    ref_slit = trmodels.Slit(55, 9376, 1, 251, 26, -5.6, 1.0, 4, 1, '1111x', '95065_1', '2122',
                             0.13, -0.31716078999999997, -0.18092266)
    _compare_slits(slitlet_info[0], ref_slit)


def test_msa_configuration_no_background():
    """
    Test the get_open_msa_slits function.
    """
    # Test 2: Two main shutters, not allowed and should fail
    msa_meta_id = 13
    msaconfl = get_file_path('msa_configuration.fits')
    dither_position = 1
    with pytest.raises(MSAFileError):
        nirspec.get_open_msa_slits(msaconfl, msa_meta_id, dither_position,
                                   slit_y_range=[-.5, .5])


def test_msa_configuration_all_background():
    """
    Test the get_open_msa_slits function.
    """

    # Test 3:  No non-background, not acceptable.
    msa_meta_id = 14
    msaconfl = get_file_path('msa_configuration.fits')
    dither_position = 1
    slitlet_info = nirspec.get_open_msa_slits(msaconfl, msa_meta_id, dither_position,
                                              slit_y_range=[-.5, .5])
    ref_slit = trmodels.Slit(57, 8281, 1, 251, 23, -2.15, 2.15, 4, 0, '1x1', 'background_57', 'bkg_57',
                             0, -0.5, -0.5)
    _compare_slits(slitlet_info[0], ref_slit)


def test_msa_configuration_row_skipped():
    """
    Test the get_open_msa_slits function.
    """

    # Test 4: One row is skipped, should be acceptable.
    msa_meta_id = 15
    msaconfl = get_file_path('msa_configuration.fits')
    dither_position = 1
    slitlet_info = nirspec.get_open_msa_slits(msaconfl, msa_meta_id, dither_position,
                                              slit_y_range=[-.5, .5])
    ref_slit = trmodels.Slit(58, 8646, 1, 251, 24, -3.3, 5.6, 4, 1, '11x1011', '95065_1', '2122',
                             0.130, -0.31716078999999997, -0.18092266)
    _compare_slits(slitlet_info[0], ref_slit)


def test_msa_configuration_multiple_returns():
    """
    Test the get_open_msa_slits function.
    """
    # Test 4: One row is skipped, should be acceptable.
    msa_meta_id = 16
    msaconfl = get_file_path('msa_configuration.fits')
    dither_position = 1
    slitlet_info = nirspec.get_open_msa_slits(msaconfl, msa_meta_id, dither_position,
                                              slit_y_range=[-.5, .5])
    ref_slit1 = trmodels.Slit(59, 8651, 1, 256, 24, -3.3, 5.6, 4, 1, '11x1011', '95065_1', '2122',
                              0.13000000000000003, -0.31716078999999997, -0.18092266)
    ref_slit2 = trmodels.Slit(60, 11573, 1, 258, 32, -3.3, 4.45, 4, 2, '11x111', '95065_2', '172',
                              0.70000000000000007, -0.31716078999999997, -0.18092266)
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

    hdul = create_nirspec_fs_file(grating="G395M", filter="OPAQUE", lamp="LINE1")
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
    with pytest.raises(MSAFileError):
        assign_wcs_step.AssignWcsStep.call(model)

    model.meta.instrument.msa_metadata_file = "missing.fits"
    with pytest.raises(MSAFileError):
        assign_wcs_step.AssignWcsStep.call(model)


def test_open_slits():
    """ Test that get_open_slits works with MSA data.

    Issue #2321
    """
    image = create_nirspec_mos_file()
    model = datamodels.ImageModel(image)
    msaconfl = get_file_path('msa_configuration.fits')

    model.meta.instrument.msa_metadata_file = msaconfl
    model.meta.instrument.msa_metadata_id = 12

    slits = nirspec.get_open_slits(model)
    assert len(slits) == 1


def test_shutter_size_on_sky():
    """
    Test the size of a MOS shutter on sky is ~ .2 x .4 arcsec.
    """
    image = create_nirspec_mos_file()
    model = datamodels.ImageModel(image)
    msaconfl = get_file_path('msa_configuration.fits')

    model.meta.instrument.msa_metadata_file = msaconfl
    model.meta.instrument.msa_metadata_id = 12

    refs = create_reference_files(model)

    pipe = nirspec.create_pipeline(model, refs, slit_y_range=(-.5, .5))
    w = wcs.WCS(pipe)
    model.meta.wcs = w
    slit = w.get_transform('slit_frame', 'msa_frame').slits[0]
    wslit = nirspec.nrs_wcs_set_input(model, slit.name)
    virtual_corners_x = [-.5, -.5, .5, .5, -.5]
    virtual_corners_y = [-.5, .5, .5, -.5, -.5]
    input_lam = [2e-6] * 5

    slit2world = wslit.get_transform('slit_frame', 'world')
    ra, dec, lam = slit2world(virtual_corners_x,
                              virtual_corners_y,
                              input_lam)
    sky = coords.SkyCoord(ra * u.deg, dec * u.deg)
    sep_x = sky[0].separation(sky[3]).to(u.arcsec)
    sep_y = sky[0].separation(sky[1]).to(u.arcsec)

    assert sep_x.value > 0.193
    assert sep_x.value < 0.194
    assert sep_y.value > 0.45
    assert sep_y.value < 0.46


@pytest.mark.parametrize(('mode'), ['fs', 'msa'])
def test_functional_fs_msa(mode):
    #     """
    #     Compare Nirspec instrument model with IDT model for FS and MSA.
    #     """
    if mode == 'fs':
        model_file = 'fixed_slits_functional_ESA_v4_20180618.txt'
        hdul = create_nirspec_fs_file(grating='G395H', filter='F290LP')
        im = datamodels.ImageModel(hdul)
        refs = create_reference_files(im)
        pipeline = nirspec.create_pipeline(im, refs, slit_y_range=[-0.55, 0.55])
        w = wcs.WCS(pipeline)
        im.meta.wcs = w
        # Use slit S200A1
        slit_wcs = nirspec.nrs_wcs_set_input(im, 'S200A1')

    if mode == 'msa':
        model_file = 'msa_functional_ESA_v2_20180620.txt'
        hdul = create_nirspec_mos_file(grating='G395H', filt='F290LP')
        im = datamodels.ImageModel(hdul)
        refs = create_reference_files(im)
        slit = trmodels.Slit(name=1, shutter_id=4699, xcen=319, ycen=13,
                             ymin=-0.55000000000000004, ymax=0.55000000000000004,
                             quadrant=3, source_id=1, shutter_state='x',
                             source_name='lamp', source_alias='foo',
                             stellarity=100.0, source_xpos=-0.5, source_ypos=0.5)
        open_slits = [slit]
        pipeline = nirspec.slitlets_wcs(im, refs, open_slits)
        w = wcs.WCS(pipeline)
        im.meta.wcs = w
        slit_wcs = nirspec.nrs_wcs_set_input(im, 1)

    ins_file = get_file_path(model_file)
    ins_tab = table.Table.read(ins_file, format='ascii')

    # Setup the test
    slitx = [0] * 5
    slity = [-.5, -.25, 0, .25, .5]
    lam = np.array([2.9, 3.39, 3.88, 4.37, 5]) * 10**-6

    # Slit to MSA absolute
    slit2msa = slit_wcs.get_transform('slit_frame', 'msa_frame')
    msax, msay, _ = slit2msa(slitx, slity, lam)

    assert_allclose(slitx, ins_tab['xslitpos'])
    assert_allclose(slity, ins_tab['yslitpos'])
    assert_allclose(msax, ins_tab['xmsapos'])
    assert_allclose(msay, ins_tab['ymaspos'])

    # Coordinates at Collimator exit
    # Applies the Collimator forward transform to MSa absolute coordinates
    with datamodels.open(refs['collimator']) as col:
        colx, coly = col.model.inverse(msax, msay)
    assert_allclose(colx, ins_tab['xcoll'])
    assert_allclose(coly, ins_tab['ycoll'])

    # After applying directional cosines
    dircos = trmodels.Unitless2DirCos()
    xcolDircosi, ycolDircosi, z = dircos(colx, coly)
    assert_allclose(xcolDircosi, ins_tab['xcolDirCosi'])
    assert_allclose(ycolDircosi, ins_tab['ycolDirCosi'])

    # MSA to GWA entrance
    # This runs the Collimator forward, Unitless to Directional cosine, and
    # 3D Rotation. It uses the corrected GWA tilt value
    with datamodels.DisperserModel(refs['disperser']) as disp:
        disperser = nirspec.correct_tilt(disp, im.meta.instrument.gwa_xtilt,
                                         im.meta.instrument.gwa_ytilt)
    collimator2gwa = nirspec.collimator_to_gwa(refs, disperser)
    x_gwa_in, y_gwa_in, z_gwa_in = collimator2gwa(msax, msay)
    assert_allclose(x_gwa_in, ins_tab['xdispIn'])
    assert_allclose(y_gwa_in, ins_tab['ydispIn'])

    # Slit to GWA out
    slit2gwa = slit_wcs.get_transform('slit_frame', 'gwa')
    x_gwa_out, y_gwa_out, z_gwa_out = slit2gwa(slitx, slity, lam)
    assert_allclose(x_gwa_out, ins_tab['xdispLaw'])
    assert_allclose(y_gwa_out, ins_tab['ydispLaw'])

    # CAMERA entrance (assuming direction is from sky to detector)
    angles = [disperser['theta_x'], disperser['theta_y'],
              disperser['theta_z'], disperser['tilt_y']]
    rotation = trmodels.Rotation3DToGWA(angles, axes_order="xyzy",
                                        name='rotation')
    dircos2unitless = trmodels.DirCos2Unitless()
    gwa2cam = rotation.inverse | dircos2unitless
    x_camera_entrance, y_camera_entrance = gwa2cam(x_gwa_out, y_gwa_out,
                                                   z_gwa_out)
    assert_allclose(x_camera_entrance, ins_tab['xcamCosi'])
    assert_allclose(y_camera_entrance, ins_tab['ycamCosi'])

    # at FPA
    with datamodels.CameraModel(refs['camera']) as camera:
        x_fpa, y_fpa = camera.model.inverse(x_camera_entrance, y_camera_entrance)
    assert_allclose(x_fpa, ins_tab['xfpapos'])
    assert_allclose(y_fpa, ins_tab['yfpapos'])

    # at SCA These are 0-based , the IDT results are 1-based
    slit2sca = slit_wcs.get_transform('slit_frame', 'sca')
    x_sca_nrs1, y_sca_nrs1 = slit2sca(slitx, slity, lam)
    # At NRS2
    with datamodels.FPAModel(refs['fpa']) as fpa:
        x_sca_nrs2, y_sca_nrs2 = fpa.nrs2_model.inverse(x_fpa, y_fpa)
    # expect 1 pix difference
    wvlns_on_nrs1 = slice(2)
    wvlns_on_nrs2 = slice(2, 4)
    assert_allclose(x_sca_nrs1[wvlns_on_nrs1] + 1, ins_tab['i'][wvlns_on_nrs1])
    assert_allclose(y_sca_nrs1[wvlns_on_nrs1] + 1, ins_tab['j'][wvlns_on_nrs1])
    assert_allclose(x_sca_nrs2[wvlns_on_nrs2] + 1, ins_tab['i'][wvlns_on_nrs2])
    assert_allclose(y_sca_nrs2[wvlns_on_nrs2] + 1, ins_tab['j'][wvlns_on_nrs2])

    # at oteip
    slit2oteip = slit_wcs.get_transform('slit_frame', 'oteip')
    x_oteip, y_oteip, _ = slit2oteip(slitx, slity, lam)
    assert_allclose(x_oteip, ins_tab['xOTEIP'])
    assert_allclose(y_oteip, ins_tab['yOTEIP'])

    # at v2, v3 [in arcsec]
    slit2v23 = slit_wcs.get_transform('slit_frame', 'v2v3')
    v2, v3, _ = slit2v23(slitx, slity, lam)
    v2 /= 3600
    v3 /= 3600
    assert_allclose(v2, ins_tab['xV2V3'])
    assert_allclose(v3, ins_tab['yV2V3'])


@pytest.fixture
def wcs_ifu_grating():
    def _create_image_model(grating='G395H', filter='F290LP', **kwargs):
        hdul = create_nirspec_ifu_file(grating=grating, filter=filter, **kwargs)
        im = datamodels.ImageModel(hdul)
        refs = create_reference_files(im)
        pipeline = nirspec.create_pipeline(im, refs, slit_y_range=[-0.5, 0.5])
        w = wcs.WCS(pipeline)
        im.meta.wcs = w
        return im, refs
    return _create_image_model


def test_functional_ifu_grating(wcs_ifu_grating):
    """Compare Nirspec instrument model with IDT model for IFU grating."""

    # setup test
    model_file = 'ifu_grating_functional_ESA_v1_20180619.txt'
    im, refs = wcs_ifu_grating('G395H', 'F290LP', gwa_xtil=0.35986012, gwa_ytil=0.13448857)

    slit_wcs = nirspec.nrs_wcs_set_input(im, 0)  # use slice 0
    ins_file = get_file_path(model_file)
    ins_tab = table.Table.read(ins_file, format='ascii')
    slitx = [0] * 5
    slity = [-.5, -.25, 0, .25, .5]
    lam = np.array([2.9, 3.39, 3.88, 4.37, 5]) * 10**-6
    order, wrange = nirspec.get_spectral_order_wrange(im, refs['wavelengthrange'])
    im.meta.wcsinfo.sporder = order
    im.meta.wcsinfo.waverange_start = wrange[0]
    im.meta.wcsinfo.waverange_end = wrange[1]

    # Slit to MSA entrance
    # This includes the Slicer transform and the IFUFORE transform
    slit2msa = slit_wcs.get_transform('slit_frame', 'msa_frame')
    msax, msay, _ = slit2msa(slitx, slity, lam)
    assert_allclose(slitx, ins_tab['xslitpos'])
    assert_allclose(slity, ins_tab['yslitpos'])
    assert_allclose(msax + 0.0073, ins_tab['xmsapos'], rtol=1e-2)  # expected offset
    assert_allclose(msay + 0.0085, ins_tab['ymaspos'], rtol=1e-2)  # expected offset

    # Slicer
    slit2slicer = slit_wcs.get_transform('slit_frame', 'slicer')
    x_slicer, y_slicer, _ = slit2slicer(slitx, slity, lam)

    # MSA exit
    # Applies the IFUPOST transform to coordinates at the Slicer
    with datamodels.IFUPostModel(refs['ifupost']) as ifupost:
        ifupost_transform = nirspec._create_ifupost_transform(ifupost.slice_0)
    x_msa_exit, y_msa_exit = ifupost_transform(x_slicer, y_slicer, lam)
    assert_allclose(x_msa_exit, ins_tab['xmsapos'])
    assert_allclose(y_msa_exit, ins_tab['ymaspos'])

    # Computations are done using the exact form of the equations in the reports
    # Part I of the Forward IFU-POST transform - the linear transform
    xc_out = 0.0487158154447
    yc_out = 0.00856211956976
    xc_in = 0.000355277216
    yc_in = -3.0089012e-05
    theta = np.deg2rad(-0.129043957046)
    factor_x = 0.100989874454
    factor_y = 0.100405184145

    # Slicer coordinates
    xS = 0.000399999989895
    yS = -0.00600000005215

    x = xc_out + factor_x * (+cos(theta) * (xS - xc_in) + sin(theta) * (yS - yc_in))
    y = yc_out + factor_y * (-sin(theta) * (xS - xc_in) + cos(theta) * (yS - yc_in))

    # Forward IFU-POST II part - non-linear transform
    lam = 2.9e-6
    coef_names = [f'c{x}_{y}' for x in range(6) for y in range(6) if x + y <= 5]
    y_forw = [-82.3492267824, 29234.6982762, -540260.780853, 771881.305018,
              -2563462.26848, 29914272.1164, 4513.04082605, -2212869.44311,
              32875633.0303, -29923698.5288, 27293902.5636, -39820.4434726,
              62431493.9962, -667197265.033, 297253538.182, -1838860.86305,
              -777169857.2, 4514693865.7, 42790637.764, 3596423850.94,
              -260274017.448]
    y_forw_dist = [188531839.97, -43453434864.0, 70807756765.8,
                   -308272809909.0, 159768473071.0, 9712633344590.0,
                   -11762923852.9, 3545938873190.0, -4198643655420.0,
                   12545642983100.0, -11707051591600.0, 173091230285.0,
                   -108534069056000.0, 82893348097600.0, -124708740989000.0,
                   2774389757990.0, 1476779720300000.0, -545358301961000.0,
                   -93101557994100.0, -7536890639430000.0, 646310545048000.0]
    y_coeff = {}
    for i, coef in enumerate(coef_names):
        y_coeff[coef] = y_forw[i] + lam * y_forw_dist[i]
    poly2d = astmodels.Polynomial2D(5, **y_coeff)
    ifupost_y = poly2d(x, y)
    assert_allclose(ifupost_y, ins_tab['ymaspos'][0])
    assert_allclose(ifupost_y, y_msa_exit[0])

    # reset 'lam'
    lam = np.array([2.9, 3.39, 3.88, 4.37, 5]) * 10**-6

    # Coordinates at Collimator exit
    # Applies the Collimator forward transform to coordinates at the MSA exit
    with datamodels.open(refs['collimator']) as col:
        colx, coly = col.model.inverse(x_msa_exit, y_msa_exit)
    assert_allclose(colx, ins_tab['xcoll'])
    assert_allclose(coly, ins_tab['ycoll'])

    # After applying directional cosines
    dircos = trmodels.Unitless2DirCos()
    xcolDircosi, ycolDircosi, z = dircos(colx, coly)
    assert_allclose(xcolDircosi, ins_tab['xcolDirCosi'])
    assert_allclose(ycolDircosi, ins_tab['ycolDirCosi'])

    # Slit to GWA entrance
    # applies the Collimator forward, Unitless to Directional and 3D Rotation to MSA exit coordinates
    with datamodels.DisperserModel(refs['disperser']) as disp:
        disperser = nirspec.correct_tilt(disp, im.meta.instrument.gwa_xtilt,
                                         im.meta.instrument.gwa_ytilt)
    collimator2gwa = nirspec.collimator_to_gwa(refs, disperser)
    x_gwa_in, y_gwa_in, z_gwa_in = collimator2gwa(x_msa_exit, y_msa_exit)
    assert_allclose(x_gwa_in, ins_tab['xdispIn'])

    # Slit to GWA out
    # Runs slit--> slicer --> msa_exit --> collimator --> dircos --> rotation --> angle_from_grating equation
    slit2gwa = slit_wcs.get_transform('slit_frame', 'gwa')
    x_gwa_out, y_gwa_out, z_gwa_out = slit2gwa(slitx, slity, lam)
    assert_allclose(x_gwa_out, ins_tab['xdispLaw'])
    assert_allclose(y_gwa_out, ins_tab['ydispLaw'])

    # CAMERA entrance (assuming direction is from sky to detector)
    angles = [disperser['theta_x'], disperser['theta_y'],
              disperser['theta_z'], disperser['tilt_y']]
    rotation = trmodels.Rotation3DToGWA(angles, axes_order="xyzy", name='rotation')
    dircos2unitless = trmodels.DirCos2Unitless()
    gwa2cam = rotation.inverse | dircos2unitless
    x_camera_entrance, y_camera_entrance = gwa2cam(x_gwa_out, y_gwa_out, z_gwa_out)
    assert_allclose(x_camera_entrance, ins_tab['xcamCosi'])
    assert_allclose(y_camera_entrance, ins_tab['ycamCosi'])

    # at FPA
    with datamodels.CameraModel(refs['camera']) as camera:
        x_fpa, y_fpa = camera.model.inverse(x_camera_entrance, y_camera_entrance)
    assert_allclose(x_fpa, ins_tab['xfpapos'])
    assert_allclose(y_fpa, ins_tab['yfpapos'])

    # at SCA
    slit2sca = slit_wcs.get_transform('slit_frame', 'sca')
    x_sca_nrs1, y_sca_nrs1 = slit2sca(slitx, slity, lam)

    # At NRS2
    with datamodels.FPAModel(refs['fpa']) as fpa:
        x_sca_nrs2, y_sca_nrs2 = fpa.nrs2_model.inverse(x_fpa, y_fpa)
    assert_allclose(x_sca_nrs1[:3] + 1, ins_tab['i'][:3])
    assert_allclose(y_sca_nrs1[:3] + 1, ins_tab['j'][:3])
    assert_allclose(x_sca_nrs2[3:] + 1, ins_tab['i'][3:])
    assert_allclose(y_sca_nrs2[3:] + 1, ins_tab['j'][3:])

    # at oteip
    # Goes through slicer, ifufore, and fore transforms
    slit2oteip = slit_wcs.get_transform('slit_frame', 'oteip')
    x_oteip, y_oteip, _ = slit2oteip(slitx, slity, lam)
    assert_allclose(x_oteip, ins_tab['xOTEIP'])
    assert_allclose(y_oteip, ins_tab['yOTEIP'])

    # at v2, v3 [in arcsec]
    slit2v23 = slit_wcs.get_transform('slit_frame', 'v2v3')
    v2, v3, _ = slit2v23(slitx, slity, lam)
    v2 /= 3600
    v3 /= 3600
    assert_allclose(v2, ins_tab['xV2V3'])
    assert_allclose(v3, ins_tab['yV2V3'])


def test_functional_ifu_prism():
    """Compare Nirspec instrument model with IDT model for IFU prism."""
    # setup test
    model_file = 'ifu_prism_functional_ESA_v1_20180619.txt'
    hdu1 = create_nirspec_ifu_file(grating='PRISM', filter='CLEAR',
                                   gwa_xtil=0.35986012, gwa_ytil=0.13448857,
                                   gwa_tilt=37.1)
    im = datamodels.ImageModel(hdu1)
    refs = create_reference_files(im)
    pipeline = nirspec.create_pipeline(im, refs, slit_y_range=[-0.55, 0.55])
    w = wcs.WCS(pipeline)
    im.meta.wcs = w
    slit_wcs = nirspec.nrs_wcs_set_input(im, 0)  # use slice 0
    ins_file = get_file_path(model_file)
    ins_tab = table.Table.read(ins_file, format='ascii')
    slitx = [0] * 5
    slity = [-.5, -.25, 0, .25, .5]
    lam = np.array([.7e-7, 1e-6, 2e-6, 3e-6, 5e-6])
    order, wrange = nirspec.get_spectral_order_wrange(im, refs['wavelengthrange'])
    im.meta.wcsinfo.sporder = order
    im.meta.wcsinfo.waverange_start = wrange[0]
    im.meta.wcsinfo.waverange_end = wrange[1]

    # Slit to MSA entrance
    # This includes the Slicer transform and the IFUFORE transform
    slit2msa = slit_wcs.get_transform('slit_frame', 'msa_frame')
    msax, msay, _ = slit2msa(slitx, slity, lam)
    assert_allclose(slitx, ins_tab['xslitpos'])
    assert_allclose(slity, ins_tab['yslitpos'])
    assert_allclose(msax + 0.0073, ins_tab['xmsapos'], rtol=1e-2)  # expected offset
    assert_allclose(msay + 0.0085, ins_tab['ymaspos'], rtol=1e-2)  # expected offset

    # Slicer
    slit2slicer = slit_wcs.get_transform('slit_frame', 'slicer')
    x_slicer, y_slicer, _ = slit2slicer(slitx, slity, lam)

    # MSA exit
    # Applies the IFUPOST transform to coordinates at the Slicer
    with datamodels.IFUPostModel(refs['ifupost']) as ifupost:
        ifupost_transform = nirspec._create_ifupost_transform(ifupost.slice_0)
    x_msa_exit, y_msa_exit = ifupost_transform(x_slicer, y_slicer, lam)
    assert_allclose(x_msa_exit, ins_tab['xmsapos'])
    assert_allclose(y_msa_exit, ins_tab['ymaspos'])

    # Coordinates at Collimator exit
    # Applies the Collimator forward transform to coordinates at the MSA exit
    with datamodels.open(refs['collimator']) as col:
        colx, coly = col.model.inverse(x_msa_exit, y_msa_exit)
    assert_allclose(colx, ins_tab['xcoll'])
    assert_allclose(coly, ins_tab['ycoll'])

    # After applying directional cosines
    dircos = trmodels.Unitless2DirCos()
    xcolDircosi, ycolDircosi, z = dircos(colx, coly)
    assert_allclose(xcolDircosi, ins_tab['xcolDirCosi'])
    assert_allclose(ycolDircosi, ins_tab['ycolDirCosi'])

    # Slit to GWA entrance
    # applies the Collimator forward, Unitless to Directional and 3D Rotation to MSA exit coordinates
    with datamodels.DisperserModel(refs['disperser']) as disp:
        disperser = nirspec.correct_tilt(disp, im.meta.instrument.gwa_xtilt,
                                         im.meta.instrument.gwa_ytilt)
    collimator2gwa = nirspec.collimator_to_gwa(refs, disperser)
    x_gwa_in, y_gwa_in, z_gwa_in = collimator2gwa(x_msa_exit, y_msa_exit)
    assert_allclose(x_gwa_in, ins_tab['xdispIn'])

    # Slit to GWA out
    # Runs slit--> slicer --> msa_exit --> collimator --> dircos --> rotation --> angle_from_grating equation
    slit2gwa = slit_wcs.get_transform('slit_frame', 'gwa')
    x_gwa_out, y_gwa_out, z_gwa_out = slit2gwa(slitx, slity, lam)
    assert_allclose(x_gwa_out, ins_tab['xdispLaw'])
    assert_allclose(y_gwa_out, ins_tab['ydispLaw'])

    # CAMERA entrance (assuming direction is from sky to detector)
    angles = [disperser['theta_x'], disperser['theta_y'],
              disperser['theta_z'], disperser['tilt_y']]
    rotation = trmodels.Rotation3DToGWA(angles, axes_order="xyzy", name='rotation')
    dircos2unitless = trmodels.DirCos2Unitless()
    gwa2cam = rotation.inverse | dircos2unitless
    x_camera_entrance, y_camera_entrance = gwa2cam(x_gwa_out, y_gwa_out, z_gwa_out)
    assert_allclose(x_camera_entrance, ins_tab['xcamCosi'])
    assert_allclose(y_camera_entrance, ins_tab['ycamCosi'])

    # at FPA
    with datamodels.CameraModel(refs['camera']) as camera:
        x_fpa, y_fpa = camera.model.inverse(x_camera_entrance, y_camera_entrance)
    assert_allclose(x_fpa, ins_tab['xfpapos'])
    assert_allclose(y_fpa, ins_tab['yfpapos'])

    # at SCA
    slit2sca = slit_wcs.get_transform('slit_frame', 'sca')
    x_sca_nrs1, y_sca_nrs1 = slit2sca(slitx, slity, lam)

    # At NRS2
    with datamodels.FPAModel(refs['fpa']) as fpa:
        x_sca_nrs2, y_sca_nrs2 = fpa.nrs2_model.inverse(x_fpa, y_fpa)
    assert_allclose(x_sca_nrs1 + 1, ins_tab['i'])
    assert_allclose(y_sca_nrs1 + 1, ins_tab['j'])

    # at oteip
    # Goes through slicer, ifufore, and fore transforms
    slit2oteip = slit_wcs.get_transform('slit_frame', 'oteip')
    x_oteip, y_oteip, _ = slit2oteip(slitx, slity, lam)
    assert_allclose(x_oteip, ins_tab['xOTEIP'])
    assert_allclose(y_oteip, ins_tab['yOTEIP'])

    # at v2, v3 [in arcsec]
    slit2v23 = slit_wcs.get_transform('slit_frame', 'v2v3')
    v2, v3, _ = slit2v23(slitx, slity, lam)
    v2 /= 3600
    v3 /= 3600
    assert_allclose(v2, ins_tab['xV2V3'])
    assert_allclose(v3, ins_tab['yV2V3'])


def test_ifu_bbox():
    bbox = {0: ((122.0908542999878, 1586.2584665188083),
                (773.5411133037417, 825.1150258966278)),
            1: ((140.3793485788431, 1606.8904629423566),
                (1190.353197027459, 1243.0853605832503)),
            2: ((120.0139534379125, 1583.9271768905855),
                (724.3249534782219, 775.8104288584977)),
            3: ((142.50252648927454, 1609.3106221382388),
                (1239.4122720740888, 1292.288713688988)),
            4: ((117.88884113088403, 1581.5517394150106),
                (674.9787657901347, 726.3752061973377)),
            5: ((144.57465414462143, 1611.688447569682),
                (1288.4808318659427, 1341.5035313084197)),
            6: ((115.8602297714846, 1579.27471654949),
                (625.7982466386104, 677.1147840452901)),
            7: ((146.7944728147906, 1614.2161842198498),
                (1337.531525654835, 1390.7050687363856)),
            8: ((113.86384530944383, 1577.0293086386203),
                (576.5344359685643, 627.777022204828)),
            9: ((149.0259581360621, 1616.7687282225652),
                (1386.5118806905086, 1439.843598490326)),
            10: ((111.91564190274217, 1574.8351095461135),
                 (527.229828693075, 578.402894851317)),
            11: ((151.3053466801954, 1619.3720722471498),
                 (1435.423685040875, 1488.917203728964)),
            12: ((109.8957204607345, 1572.570246400894),
                 (477.9699083444277, 529.0782087498488)),
            13: ((153.5023503173659, 1621.9005029476564),
                 (1484.38405923062, 1538.0443479389924)),
            14: ((107.98320121613297, 1570.411787034636),
                 (428.6704834494425, 479.7217241891257)),
            15: ((155.77991404913857, 1624.5184927460925),
                 (1533.169633314481, 1586.9984359105376)),
            16: ((106.10212081215678, 1568.286103827344),
                 (379.3860245240618, 430.3780648366697)),
            17: ((158.23149941845386, 1627.305849064835),
                 (1582.0496119714928, 1636.0513450787032)),
            18: ((104.09366374413436, 1566.030231370944),
                 (330.0822744105267, 381.01974582564395)),
            19: ((160.4511021152353, 1629.888830991371),
                 (1630.7797743277185, 1684.9592727079018)),
            20: ((102.25220592881234, 1563.9475099032868),
                 (280.7233309522168, 331.6093009077988)),
            21: ((162.72784286205734, 1632.5257403739463),
                 (1679.6815760587567, 1734.03692957156)),
            22: ((100.40115742738622, 1561.8476640376036),
                 (231.35443588323855, 282.19575854747006)),
            23: ((165.05939163941662, 1635.2270773628682),
                 (1728.511467615387, 1783.0485841263735)),
            24: ((98.45723949658425, 1559.6499479349648),
                 (182.0417295679079, 232.83530870639865)),
            25: ((167.44628840053574, 1637.9923229870349),
                 (1777.2512197664128, 1831.971115503598)),
            26: ((96.56508092457855, 1557.5079027818058),
                 (132.5285162704088, 183.27350269292484)),
            27: ((169.8529496136358, 1640.778485168005),
                 (1826.028691168028, 1880.9336718824313)),
            28: ((94.71390837793813, 1555.4048050512263),
                 (82.94691422559131, 133.63901517357235)),
            29: ((172.3681094850081, 1643.685604697228),
                 (1874.8184744639657, 1929.9072657798927))}

    hdul = create_nirspec_ifu_file("F290LP", "G140M")
    im = datamodels.IFUImageModel(hdul)
    im.meta.filename = "test_ifu.fits"
    refs = create_reference_files(im)

    pipe = nirspec.create_pipeline(im, refs, slit_y_range=[-.5, .5])
    w = wcs.WCS(pipe)
    im.meta.wcs = w

    _, wrange = nirspec.spectral_order_wrange_from_model(im)
    pipe = im.meta.wcs.pipeline

    g2s = pipe[2].transform
    transforms = [pipe[0].transform]
    transforms.append(pipe[1].transform[1:])
    transforms.append(astmodels.Identity(1))
    transforms.append(astmodels.Identity(1))
    transforms.extend([step.transform for step in pipe[4:-1]])

    for sl in range(30):
        transforms[2] = g2s.get_model(sl)
        m = functools.reduce(lambda x, y: x | y, [tr.inverse for tr in transforms[:3][::-1]])
        bbox_sl = nirspec.compute_bounding_box(m, wrange)
        assert_allclose(bbox[sl], bbox_sl)


@pytest.fixture
def ifu_world_coord(wcs_ifu_grating):
    """ Return RA, DEC, LAM for all slices in the NRS IFU."""
    ra_all = []
    dec_all = []
    lam_all = []
    im, refs = wcs_ifu_grating(grating="G140H", filter="F100LP")
    for sl in range(30):
        slice_wcs = nirspec.nrs_wcs_set_input(im, sl)
        x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
        r, d, lam = slice_wcs(x, y)
        ra_all.append(r)
        dec_all.append(d)
        lam_all.append(lam)
    ra_all = np.concatenate([r.flatten() for r in ra_all])
    dec_all = np.concatenate([r.flatten() for r in dec_all])
    lam_all = np.concatenate([r.flatten() for r in lam_all])
    return ra_all, dec_all, lam_all


@pytest.mark.parametrize("slice", [1, 17])
def test_in_slice(slice, wcs_ifu_grating, ifu_world_coord):
    """ Test that the number of valid outputs from a slice forward transform
    equals the valid pixels within the slice from the slice backward transform.
    """
    ra_all, dec_all, lam_all = ifu_world_coord
    im, refs = wcs_ifu_grating("G140H", "F100LP")
    slice_wcs = nirspec.nrs_wcs_set_input(im, slice)
    slicer2world = slice_wcs.get_transform('slicer','world')
    detector2slicer = slice_wcs.get_transform('detector','slicer')
    x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
    onslice_ind = in_ifu_slice(slice_wcs, ra_all, dec_all, lam_all)
    slx, sly, sllam = slicer2world.inverse(ra_all, dec_all, lam_all)
    xinv,yinv = detector2slicer.inverse(slx[onslice_ind], sly[onslice_ind],
                                        sllam[onslice_ind])

    r, d, _ = slice_wcs(x, y)
    assert r[~np.isnan(r)].size == xinv.size
