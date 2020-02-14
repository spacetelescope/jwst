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

from ..nircam import tsgrism

# Allowed settings for nircam
TSGRISM_FILTERS = ['F277W', 'F444W', 'F322W2', 'F356W']

NIRCAM_WFSS_FRAMES = ['grism_detector', 'detector', 'v2v3', 'world']

NIRCAM_TSGRISM_FRAMES = ['grism_detector', 'full_detector', 'v2v3', 'world']

NIRCAM_IMAGING_FRAMES = ['detector', 'v2v3', 'world']

# Default primary HDU information
DEFAULT_PHDU = {
    'telescop': 'JWST',
    'instrume': 'NIRCAM',
    'channel': 'LONG',
    'detector': 'NRCALONG',
    'FILTER': 'F444W',
    'PUPIL': 'GRISMR',
    'MODULE': 'A',
    'time-obs': '8:59:37',
    'date-obs': '2017-09-05',
    'exp_type': 'NRC_IMAGE'
}

# Default wcs information
DEFAULT_WCS_WFSS = {
    'wcsaxes': 2,
    'ra_ref': 53.1423683802,
    'dec_ref': -27.8171119969,
    'v2_ref': 86.103458,
    'v3_ref': -493.227512,
    'roll_ref': 45.04234459270135,
    'crpix1': 1024.5,
    'crpix2': 1024.5,
    'crval1': 53.1423683802,
    'crval2': -27.8171119969,
    'cdelt1': 1.74460027777777e-05,
    'cdelt2': 1.75306861111111e-05,
    'ctype1': 'RA---TAN',
    'ctype2': 'DEC--TAN',
    'pc1_1': -1,
    'pc1_2': 0,
    'pc2_1': 0,
    'pc2_2': 1,
    'cunit1': 'deg',
    'cunit2': 'deg',
}

DEFAULT_WCS_TSO = {
    'wcsaxes': 2,
    'ra_ref': 86.9875,
    'dec_ref': 23.423,
    'v2_ref': 95.043034,
    'v3_ref': -556.150466,
    'roll_ref': 359.9521,
    'xref_sci': 887.0,
    'yref_sci': 35.0,
    'cdelt1': 1.76686111111111e-05,
    'cdelt2': 1.78527777777777e-05,
    'ctype1': 'RA---TAN',
    'ctype2': 'DEC--TAN',
    'pc1_1': -1,
    'pc1_2': 0,
    'pc2_1': 0,
    'pc2_2': 1,
    'cunit1': 'deg',
    'cunit2': 'deg',
}


@pytest.mark.parametrize('pupil', ('GRISMR', 'GRISMC'))
def test_nircam_wfss_available_frames(pupil, assign_wcs_objects):
    """Make sure that the expected GWCS reference frames are created."""
    sci_hdu_keys = {**{'SUBARRAY': 'FULL'}, **DEFAULT_WCS_WFSS}

    primary_keys = DEFAULT_PHDU.copy()
    primary_keys.update({'PUPIL': pupil, 'exp_type': 'NRC_WFSS'})

    bundle = assign_wcs_objects(
        primary_hdu_keys=primary_keys,
        sci_hdu_keys=sci_hdu_keys,
        model_class_name='ImageModel'
    )

    assert all([a == b for a, b in zip(NIRCAM_WFSS_FRAMES, bundle.wcs_object.available_frames)])


def test_nircam_tso_available_frames(assign_wcs_objects):
    """Make sure that the expected GWCS reference frames for TSO are created."""
    sci_hdu_keys = {**{'SUBARRAY': 'SUBGRISM256'}, **DEFAULT_WCS_TSO}

    primary_keys = DEFAULT_PHDU.copy()
    primary_keys.update({'PUPIL': 'GRISMR', 'detector': 'NRCALONG', 'exp_type': 'NRC_TSGRISM'})

    bundle = assign_wcs_objects(
        primary_hdu_keys=primary_keys,
        sci_hdu_keys=sci_hdu_keys,
        model_class_name='CubeModel'
    )

    assert all([a == b for a, b in zip(NIRCAM_TSGRISM_FRAMES, bundle.wcs_object.available_frames)])


@pytest.mark.parametrize('key', ['xref_sci', 'yref_sci'])
def test_extract_tso_object_fails_without_xref_yref(assign_wcs_objects, key):
    sci_hdu_keys = {**{'SUBARRAY': 'FULL'}, **DEFAULT_WCS_WFSS}
    sci_hdu_keys.update({key: None})

    bundle = assign_wcs_objects(
        primary_hdu_keys=DEFAULT_PHDU,
        sci_hdu_keys=sci_hdu_keys,
        model_class_name='ImageModel'
    )

    with pytest.raises(ValueError):
        tsgrism(bundle.datamodel, bundle.references)


@pytest.mark.parametrize('pupil', ('GRISMR', 'GRISMC'))
def test_traverse_wfss_grisms(pupil, assign_wcs_objects):
    """Make sure that the WFSS dispersion polynomials are reversable."""
    sci_hdu_keys = {**{'SUBARRAY': 'FULL'}, **DEFAULT_WCS_WFSS}

    primary_keys = DEFAULT_PHDU.copy()
    primary_keys.update({'PUPIL': pupil, 'exp_type': 'NRC_WFSS'})

    bundle = assign_wcs_objects(
        primary_hdu_keys=primary_keys,
        sci_hdu_keys=sci_hdu_keys,
        model_class_name='ImageModel'
    )

    detector_to_grism = bundle.wcs_object.get_transform('detector', 'grism_detector')
    grism_to_detector = bundle.wcs_object.get_transform('grism_detector', 'detector')

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


def test_traverse_tso_grism(assign_wcs_objects):
    """Make sure that the TSO dispersion polynomials are reversable."""
    sci_hdu_keys = {**{'SUBARRAY': 'SUBGRISM256'}, **DEFAULT_WCS_TSO}

    primary_keys = DEFAULT_PHDU.copy()
    primary_keys.update({'PUPIL': 'GRISMR', 'detector': 'NRCALONG', 'exp_type': 'NRC_TSGRISM'})

    bundle = assign_wcs_objects(
        primary_hdu_keys=primary_keys,
        sci_hdu_keys=sci_hdu_keys,
        model_class_name='CubeModel'
    )

    detector_to_grism = bundle.wcs_object.get_transform('full_detector', 'grism_detector')
    grism_to_detector = bundle.wcs_object.get_transform('grism_detector', 'full_detector')

    # TSGRISM always has same source locations
    # takes x,y,order -> ra, dec, wave, order
    xin, yin, order = (100, 100, 1)

    x0, y0, lam, orderdet = grism_to_detector(xin, yin, order)
    x, y, orderdet = detector_to_grism(x0, y0, lam, order)

    assert x0 == DEFAULT_WCS_TSO['xref_sci']
    assert y0 == DEFAULT_WCS_TSO['yref_sci']
    assert order == orderdet
    assert_allclose(x, xin)
    assert y == DEFAULT_WCS_TSO['yref_sci']


def test_imaging_frames(assign_wcs_objects):
    """Verify the available imaging mode reference frames."""
    sci_hdu_keys = {**{'SUBARRAY': 'FULL'}, **DEFAULT_WCS_WFSS}

    bundle = assign_wcs_objects(
        primary_hdu_keys=DEFAULT_PHDU,
        sci_hdu_keys=sci_hdu_keys,
        model_class_name='ImageModel'
    )

    assert all([a == b for a, b in zip(NIRCAM_IMAGING_FRAMES, bundle.wcs_object.available_frames)])


@pytest.mark.xfail
def test_imaging_distortion(assign_wcs_objects):
    """Verify that the distortion correction round trips."""
    sci_hdu_keys = {**{'SUBARRAY': 'FULL'}, **DEFAULT_WCS_WFSS}

    bundle = assign_wcs_objects(
        primary_hdu_keys=DEFAULT_PHDU,
        sci_hdu_keys=sci_hdu_keys,
        model_class_name='ImageModel'
    )

    sky_to_detector = bundle.wcs_object.get_transform('world', 'detector')
    detector_to_sky = bundle.wcs_object.get_transform('detector', 'sky')

    # we'll use the crpix as the simplest reference point
    ra = DEFAULT_WCS_WFSS['crval1']
    dec = DEFAULT_WCS_WFSS['crval2']

    x, y = sky_to_detector(ra, dec)
    raout, decout = detector_to_sky(x, y)

    assert_allclose(x, DEFAULT_WCS_WFSS['crpix1'])
    assert_allclose(y, DEFAULT_WCS_WFSS['crpix2'])
    assert_allclose(raout, ra)
    assert_allclose(decout, dec)
