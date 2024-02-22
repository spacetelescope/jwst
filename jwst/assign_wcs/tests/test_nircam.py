"""
Test NIRCAM grism WCS transformations.

Notes:
These test the stability of the WFSS and TSO transformations based on the
the reference file that is returned from CRDS. The absolute validity
of the results is verified by the team based on the specwcs reference
file.

"""
from numpy.testing import assert_allclose
import pytest


from astropy.io import fits
from gwcs import wcs

from stdatamodels.jwst.datamodels import CubeModel, ImageModel

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.assign_wcs import nircam
from jwst.assign_wcs import util


# Allowed settings for nircam
tsgrism_filters = ['F277W', 'F444W', 'F322W2', 'F356W']

nircam_wfss_frames = ['grism_detector', 'detector', 'v2v3', 'v2v3vacorr', 'world']

nircam_tsgrism_frames = ['grism_detector', 'direct_image', 'v2v3', 'v2v3vacorr', 'world']

nircam_imaging_frames = ['detector', 'v2v3', 'v2v3vacorr', 'world']


# Default wcs information
wcs_wfss_kw = {'wcsaxes': 2, 'ra_ref': 53.1423683802, 'dec_ref': -27.8171119969,
               'v2_ref': 86.103458, 'v3_ref': -493.227512, 'roll_ref': 45.04234459270135,
               'crpix1': 1024.5, 'crpix2': 1024.5,
               'crval1': 53.1423683802, 'crval2': -27.8171119969,
               'cdelt1': 1.74460027777777e-05, 'cdelt2': 1.75306861111111e-05,
               'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
               'pc1_1': -1, 'pc1_2': 0,
               'pc2_1': 0, 'pc2_2': 1,
               'cunit1': 'deg', 'cunit2': 'deg',
               }

wcs_tso_kw = {'wcsaxes': 2, 'ra_ref': 86.9875, 'dec_ref': 23.423,
              'v2_ref': 95.043034, 'v3_ref': -556.150466, 'roll_ref': 359.9521,
              'xref_sci': 887.0, 'yref_sci': 35.0,
              'cdelt1': 1.76686111111111e-05, 'cdelt2': 1.78527777777777e-05,
              'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
              'pc1_1': -1, 'pc1_2': 0,
              'pc2_1': 0, 'pc2_2': 1,
              'cunit1': 'deg', 'cunit2': 'deg',
              }


def create_hdul(detector='NRCALONG', channel='LONG', module='A',
                filtername='F444W', exptype='NRC_IMAGE', pupil='GRISMR',
                subarray='FULL', wcskeys=wcs_wfss_kw):
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['telescop'] = "JWST"
    phdu.header['filename'] = "test+" + filtername
    phdu.header['instrume'] = 'NIRCAM'
    phdu.header['channel'] = channel
    phdu.header['detector'] = detector
    phdu.header['FILTER'] = filtername
    phdu.header['PUPIL'] = pupil
    phdu.header['MODULE'] = module
    phdu.header['time-obs'] = '8:59:37'
    phdu.header['date-obs'] = '2023-01-01'
    phdu.header['exp_type'] = exptype
    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    scihdu.header['SUBARRAY'] = subarray
    scihdu.header.update(wcskeys)
    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def create_wfss_wcs(pupil, filtername='F444W'):
    """Help create WFSS GWCS object."""
    hdul = create_hdul(exptype='NRC_WFSS', filtername=filtername, pupil=pupil)
    im = ImageModel(hdul)
    ref = get_reference_files(im)
    pipeline = nircam.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def create_imaging_wcs():
    hdul = create_hdul()
    image = ImageModel(hdul)
    ref = get_reference_files(image)
    pipeline = nircam.create_pipeline(image, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def create_tso_wcs(filtername=tsgrism_filters[0], subarray="SUBGRISM256"):
    """Help create tsgrism GWCS object."""
    hdul = create_hdul(exptype='NRC_TSGRISM', pupil='GRISMR',
                       filtername=filtername, detector='NRCALONG',
                       subarray=subarray, wcskeys=wcs_tso_kw)
    im = CubeModel(hdul)
    ref = get_reference_files(im)
    pipeline = nircam.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def get_reference_files(datamodel):
    """Get the reference files associated with the step."""
    refs = {}
    step = AssignWcsStep()
    for reftype in AssignWcsStep.reference_file_types:
        val = step.get_reference_file(datamodel, reftype)
        if val == 'N/A':
            refs[reftype] = None
        else:
            refs[reftype] = val

    return refs


@pytest.fixture(params=tsgrism_filters)
def tsgrism_inputs(request):
    def _add_missing_key(missing_key=None):
        tso_kw = wcs_tso_kw.copy()

        if missing_key is not None:
            tso_kw[missing_key] = None

        hdu = create_hdul(
            exptype='NRC_TSGRISM',
            pupil='GRISMR',
            filtername=request.param,
            detector='NRCALONG',
            subarray='SUBGRISM256',
            wcskeys=tso_kw,
            channel='LONG',
            module='A',
        )

        image_model = CubeModel(hdu)

        return image_model, get_reference_files(image_model)

    return _add_missing_key


def test_nircam_wfss_available_frames():
    """Make sure that the expected GWCS reference frames are created."""
    for p in ['GRISMR', 'GRISMC']:
        wcsobj = create_wfss_wcs(p)
        available_frames = wcsobj.available_frames
        assert all([a == b for a, b in zip(nircam_wfss_frames, available_frames)])


def test_nircam_tso_available_frames():
    """Make sure that the expected GWCS reference frames for TSO are created."""
    wcsobj = create_tso_wcs()
    available_frames = wcsobj.available_frames
    assert all([a == b for a, b in zip(nircam_tsgrism_frames, available_frames)])


@pytest.mark.parametrize('key', ['xref_sci', 'yref_sci'])
def test_extract_tso_object_fails_without_xref_yref(tsgrism_inputs, key):
    with pytest.raises(ValueError):
        nircam.tsgrism(*tsgrism_inputs(missing_key=key))


def traverse_wfss_trace(pupil):
    """Make sure that the WFSS dispersion polynomials are reversable."""
    wcsobj = create_wfss_wcs(pupil)
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
    """Check the dispersion polynomials for each grism."""
    for pupil in ['GRISMR', 'GRISMC']:
        traverse_wfss_trace(pupil)


@pytest.mark.xfail(reason="Fails due to V2 NIRCam specwcs ref files delivered to CRDS")
def test_traverse_tso_grism():
    """Make sure that the TSO dispersion polynomials are reversable."""
    wcsobj = create_tso_wcs()
    detector_to_grism = wcsobj.get_transform('direct_image', 'grism_detector')
    grism_to_detector = wcsobj.get_transform('grism_detector', 'direct_image')

    # TSGRISM always has same source locations
    # takes x,y,order -> ra, dec, wave, order
    xin, yin, order = (100, 100, 1)

    x0, y0, lam, orderdet = grism_to_detector(xin, yin, order)
    x, y, orderdet = detector_to_grism(x0, y0, lam, order)

    assert x0 == wcs_tso_kw['xref_sci']
    assert y0 == wcs_tso_kw['yref_sci']
    assert order == orderdet
    assert_allclose(x, xin)
    assert y == wcs_tso_kw['yref_sci']


def test_imaging_frames():
    """Verify the available imaging mode reference frames."""
    wcsobj = create_imaging_wcs()
    available_frames = wcsobj.available_frames
    assert all([a == b for a, b in zip(nircam_imaging_frames, available_frames)])


def test_wfss_sip():
    hdul = create_hdul(filtername='F444W', pupil='GRISMR', exptype='NRC_WFSS')
    wfss_model = ImageModel(hdul)
    ref = get_reference_files(wfss_model)
    pipeline = nircam.create_pipeline(wfss_model, ref)
    wcsobj = wcs.WCS(pipeline)
    wfss_model.meta.wcs = wcsobj
    util.wfss_imaging_wcs(wfss_model, nircam.imaging, bbox=((1, 1024), (1, 1024)))
    for key in ['a_order', 'b_order', 'crpix1', 'crpix2', 'crval1', 'crval2', 'cd1_1']:
        assert key in wfss_model.meta.wcsinfo.instance
