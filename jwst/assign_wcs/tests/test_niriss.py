"""
Test NIRISS grism WCS transformations.

Notes:
These test the stability of the WFSS transformations based on the
the reference file that is returned from CRDS. The absolute validity
of the results is verified by the team based on the specwcs reference
file.

"""

from numpy.testing import assert_allclose
from astropy.io import fits
from gwcs import wcs

from jwst.datamodels.image import ImageModel

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.assign_wcs import niriss
from jwst.assign_wcs import util


# Allowed settings for niriss
niriss_wfss_frames = ['grism_detector', 'detector', 'v2v3', 'v2v3vacorr', 'world']
niriss_imaging_frames = ['detector', 'v2v3', 'v2v3vacorr', 'world']
niriss_grisms = ['GR150R', 'GR150C']

# default wcs information
wcs_kw = {'wcsaxes': 2, 'ra_ref': 53.16255, 'dec_ref': -27.791461111111,
          'v2_ref': -289.966976, 'v3_ref': -697.723326, 'roll_ref': 0,
          'crval1': 53.16255, 'crval2': -27.791461111111,
          'crpix1': 1024.5, 'crpix2': 1024.5,
          'cdelt1': 0.065398, 'cdelt2': 0.065893,
          'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
          'pc1_1': 1, 'pc1_2': 0,
          'pc2_1': 0, 'pc2_2': 1,
          'cunit1': 'deg', 'cunit2': 'deg',
          }


def create_hdul(detector='NIS', filtername='CLEAR',
                exptype='NIS_IMAGE', pupil='F200W'):
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['telescop'] = "JWST"
    phdu.header['filename'] = "test+" + filtername
    phdu.header['instrume'] = 'NIRISS'
    phdu.header['detector'] = detector
    phdu.header['FILTER'] = filtername
    phdu.header['PUPIL'] = pupil
    phdu.header['time-obs'] = '8:59:37'
    phdu.header['date-obs'] = '2022-09-05'
    phdu.header['exp_type'] = exptype
    phdu.header['FWCPOS'] = 354.2111
    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    scihdu.header.update(wcs_kw)
    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def create_wcsobj(hdul):
    im = ImageModel(hdul)
    ref = get_reference_files(im)
    pipeline = niriss.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def create_wfss_wcs(filtername, pupil='F200W'):
    hdul = create_hdul(filtername=filtername, pupil=pupil, exptype='NIS_WFSS')
    im = ImageModel(hdul)
    ref = get_reference_files(im)
    pipeline = niriss.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def create_imaging_wcs(filtername):
    hdul = create_hdul(filtername=filtername, pupil='CLEAR')
    im = ImageModel(hdul)
    ref = get_reference_files(im)
    pipeline = niriss.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def get_reference_files(datamodel):
    refs = {}
    step = AssignWcsStep()
    for reftype in AssignWcsStep.reference_file_types:
        val = step.get_reference_file(datamodel, reftype)
        print(reftype, val)
        if val.strip() == 'N/A':
            refs[reftype] = None
        else:
            refs[reftype] = val
    print(refs)
    return refs


def test_niriss_wfss_available_frames():
    for f in ['GR150R', 'GR150C']:
        wcsobj = create_wfss_wcs(f)
        available_frames = wcsobj.available_frames
        assert all([a == b for a, b in zip(niriss_wfss_frames, available_frames)])


def traverse_wfss_trace(filtername):
    wcsobj = create_wfss_wcs(filtername)
    detector_to_grism = wcsobj.get_transform('detector', 'grism_detector')
    grism_to_detector = wcsobj.get_transform('grism_detector', 'detector')

    # check the round trip, grism pixel 100,100, source at 110,110,order 1
    xgrism, ygrism, xsource, ysource, orderin = (100, 100, 110, 110, 1)
    x0, y0, lam, order = grism_to_detector(xgrism, ygrism, xsource, ysource, orderin)
    x, y, xdet, ydet, orderdet = detector_to_grism(x0, y0, lam, order)

    assert x0 == xsource
    assert y0 == ysource
    assert order == orderin
    assert xdet == xsource
    assert ydet == ysource
    assert orderdet == orderin


def test_traverse_wfss_grisms():
    """Make sure the trace polynomials roundtrip for both grisms."""
    for f in niriss_grisms:
        traverse_wfss_trace(f)


def test_filter_rotation(theta=[-0.1, 0, 0.5, 20]):
    """Make sure that the filter rotation is reversable."""
    for f in niriss_grisms:
        wcsobj = create_wfss_wcs(f)
        g2d = wcsobj.get_transform('grism_detector', 'detector')
        d2g = wcsobj.get_transform('detector', 'grism_detector')
        for angle in theta:
            d2g.theta = angle
            g2d.theta = -angle
            xsource, ysource, wave, order = (110, 110, 2.3, 1)
            xgrism, ygrism, xs, ys, orderout = d2g(xsource, ysource, wave, order)
            xsdet, ysdet, wavedet, orderdet = g2d(xgrism, ygrism, xs, ys, orderout)
            assert_allclose(xsdet, xsource)
            assert_allclose(ysdet, ysource)
            assert_allclose(wavedet, wave, atol=2e-3)
            assert orderdet == order


def test_imaging_frames():
    """Verify the available imaging mode reference frames."""
    wcsobj = create_imaging_wcs('F200W')
    available_frames = wcsobj.available_frames
    assert all([a == b for a, b in zip(niriss_imaging_frames, available_frames)])


def test_imaging_distortion():
    """Verify that the distortion correction roundtrips."""
    wcsobj = create_imaging_wcs('F200W')
    sky_to_detector = wcsobj.get_transform('world', 'detector')
    detector_to_sky = wcsobj.get_transform('detector', 'world')

    # we'll use the crpix as the simplest reference point
    ra = wcs_kw['crval1']
    dec = wcs_kw['crval2']

    x, y = sky_to_detector(ra, dec)
    raout, decout = detector_to_sky(x, y)

    assert_allclose(raout, ra)
    assert_allclose(decout, dec)


def test_wfss_sip():
    hdul = create_hdul(filtername='GR150R', pupil='F200W', exptype='NIS_WFSS')
    wfss_model = ImageModel(hdul)
    ref = get_reference_files(wfss_model)
    pipeline = niriss.create_pipeline(wfss_model, ref)
    wcsobj = wcs.WCS(pipeline)
    wfss_model.meta.wcs = wcsobj
    util.wfss_imaging_wcs(wfss_model, niriss.imaging, max_pix_error=0.05, bbox=((1, 1024), (1, 1024)))
    for key in ['a_order', 'b_order', 'crpix1', 'crpix2', 'crval1', 'crval2', 'cd1_1']:
        assert key in wfss_model.meta.wcsinfo.instance
