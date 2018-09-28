"""
Test NIRCAM grism extract_2d functionality.

Notes:
No images are needed here to check the location
and size of bounding boxes.

"""
import os
from numpy.testing.utils import assert_allclose
import pytest


from astropy.io import fits
from gwcs import wcs

from ...datamodels.image import ImageModel
from ...datamodels import CubeModel
from ...assign_wcs.util import create_grism_bbox
from ...assign_wcs import AssignWcsStep

from ..extract_2d_step import Extract2dStep
from .. import grisms
from . import data


# Allowed settings for nircam
tsgrism_filters = ['F277W', 'F444W', 'F322W2', 'F356W']

data_path = os.path.split(os.path.abspath(data.__file__))[0]


# Default wcs information
# This is set for a standard nircam image just as an example
# It does not test the validity of the absolute results 
wcs_image_kw = {'wcsaxes': 2, 'ra_ref': 53.1490299775, 'dec_ref': -27.8168745624,
                'v2_ref': -290.1, 'v3_ref': -697.5, 'roll_ref': 0,
                'crpix1': 1024.5, 'crpix2': 1024.5,
                'crval1': 53.1490299775, 'crval2': -27.8168745624,
                'cdelt1': 1.81661111111111e-05, 'cdelt2': 1.8303611111111e-05,
                'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
                'pc1_1': -0.707688557183348, 'pc1_2': 0.7065245261360363,
                'pc2_1': 0.7065245261360363, 'pc2_2': 1.75306861111111e-05,
                'cunit1': 'deg', 'cunit2': 'deg',
                }

wcs_wfss_kw = {'wcsaxes': 2, 'ra_ref': 53.1423683802, 'dec_ref': -27.8171119969,
               'v2_ref': 86.103458, 'v3_ref': -493.227512, 'roll_ref': 45.04234459270135,
               'crpix1': 1024.5, 'crpix2': 1024.5,
               'crval1': 53.1423683802, 'crval2': -27.8171119969,
               'cdelt1': 1.74460027777777e-05, 'cdelt2': 1.75306861111111e-05,
               'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
               'pc1_1': -0.7076885519484576, 'pc1_2': 0.7065245313795517,
               'pc2_1': 0.7065245313795517, 'pc2_2': 0.7076885519484576,
               'cunit1': 'deg', 'cunit2': 'deg',
               }

wcs_tso_kw = {'wcsaxes': 2, 'ra_ref': 86.9875, 'dec_ref': 23.423,
              'v2_ref': 95.043034, 'v3_ref': -556.150466, 'roll_ref': 359.9521,
              'crpix1': 887.0, 'crpix2': 35.0,
              'cdelt1': 1.76686111111111e-05, 'cdelt2': 1.78527777777777e-05,
              'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
              'pc1_1': -1, 'pc1_2': 0,
              'pc2_1': 0, 'pc2_2': 1,
              'cunit1': 'deg', 'cunit2': 'deg',
              }


def get_file_path(filename):
    """
    Construct an absolute path.
    """
    return os.path.join(data_path, filename)


def create_hdul(detector='NRCALONG', channel='LONG', module='A',
                filtername='F335M', exptype='NRC_IMAGE', pupil='GRISMR',
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
    phdu.header['date-obs'] = '2017-09-05'
    phdu.header['exp_type'] = exptype
    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    scihdu.header['SUBARRAY'] = subarray
    scihdu.header.update(wcskeys)
    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def create_wfss_wcs(pupil, filtername='F335M'):
    """Help create WFSS GWCS object."""
    hdul = create_hdul(exptype='NRC_WFSS', filtername=filtername, pupil=pupil)
    im = ImageModel(hdul)
    ref = get_reference_files(im)
    pipeline = nircam.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def create_tso_wcs(filtername="F335M", subarray="SUBGRISM256"):
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
    step = Extract2dStep()
    for reftype in Extract2dStep.reference_file_types:
        refs[reftype] = step.get_reference_file(datamodel, reftype)
    return refs


def test_create_box_fits():
    """Make sure that a box is created around a source catalog object.
    This version allows use of the FITS WCS to translate the source location"""
    source_catalog = get_file_path('step_SourceCatalogStep_cat.ecsv')
    hdul = create_hdul(exptype='NRC_WFSS')
    im = ImageModel(hdul)
    aswcs = AssignWcsStep()
    refs = get_reference_files(im)
    imwcs = aswcs(im)
    imwcs.meta.source_catalog.filename = source_catalog
    test_boxes = create_grism_bbox(imwcs, refs, use_fits_wcs=True, mmag_extract=99.)
    assert len(test_boxes) == 2  # the catalog has 2 objects
    assert test_boxes[0].partial_order is False
    assert test_boxes[1].partial_order is False
    assert test_boxes[0].xcentroid > 0
    assert test_boxes[1].xcentroid > 0
    assert test_boxes[0].ycentroid > 0
    assert test_boxes[1].ycentroid > 0


@pytest.mark.xfail()
def test_create_box_gwcs():
    """Make sure that a box is created around a source catalog object.
    This version allows use of the GWCS to translate the source location.

    This is currently expected to fail because of the distortion 
    reference file. The settings and catalog used should produce
    first order trace boxes for the objects.
    """
    source_catalog = get_file_path('step_SourceCatalogStep_cat.ecsv')
    hdul = create_hdul(exptype='NRC_WFSS')
    im = ImageModel(hdul)
    aswcs = AssignWcsStep()
    imwcs = aswcs(im)
    imwcs.meta.source_catalog.filename = source_catalog
    refs = get_reference_files(im)
    test_boxes = create_grism_bbox(imwcs, refs, use_fits_wcs=False, mmag_extract=99.)
    assert len(test_boxes) == 2  # the catalog has 2 objects
    assert test_boxes[0].partial_order is False
    assert test_boxes[1].partial_order is False
    assert test_boxes[0].xcentroid > 0
    assert test_boxes[1].xcentroid > 0
    assert test_boxes[0].ycentroid > 0
    assert test_boxes[1].ycentroid > 0


def test_extract_specific_order():
    """Test that only the specified order is extracted..
     instead of the default in the reference file.

     TODO:  set use_fits_wcs to False when 
     test_create_box_gwcs stops failing
    """
    extract_orders = [1]  # just extract the first order
    source_catalog = get_file_path('step_SourceCatalogStep_cat.ecsv')
    hdul = create_hdul(exptype='NRC_WFSS')
    im = ImageModel(hdul)
    im.meta.source_catalog.filename = source_catalog
    aswcs = AssignWcsStep()
    imwcs = aswcs(im)
    refs = get_reference_files(im)
    test_boxes = create_grism_bbox(imwcs, refs, use_fits_wcs=True,
                                   mmag_extract=99.,
                                   extract_orders=extract_orders)
    assert len(test_boxes) == 2  # the catalog has 2 objects
    assert test_boxes[0].partial_order is False
    assert test_boxes[1].partial_order is False
    assert 2 not in test_boxes[0].order_bounding.keys()
    assert 2 not in test_boxes[1].order_bounding.keys()


def test_no_zero_boxes():
    """Make sure no boxes of zero size in either dimension exist."""

    pass




