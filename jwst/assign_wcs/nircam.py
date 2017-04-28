from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import Scale
from astropy.table import Table

import gwcs.coordinate_frames as cf
from . import pointing
from .util import not_implemented_mode, subarray_transform
from ..datamodels import DistortionModel


def create_pipeline(input_model, reference_files):
    """Get reference files from crds."""
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def imaging(input_model, reference_files):
    """ The NIRCAM imaging pipeline

    This includes 3 coordinate frames - detector, focal plane and sky

    reference_files={'distortion': 'test.asdf'}
    """
    detector = cf.Frame2D(name='detector',
                          axes_order=(0, 1),
                          unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.deg, u.deg))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')

    subarray2full = subarray_transform(input_model)
    imdistortion = imaging_distortion(input_model, reference_files)
    distortion = subarray2full | imdistortion
    distortion.bounding_box = imdistortion.bounding_box
    del imdistortion.bounding_box
    tel2sky = pointing.v23tosky(input_model)
    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    distortion = DistortionModel(reference_files['distortion']).model
    # Convert to deg - output of distortion models is in arcsec.
    transform = distortion | Scale(1 / 3600) & Scale(1 / 3600)

    try:
        bb = transform.bounding_box
    except NotImplementedError:
        shape = input_model.data.shape
        # Note: Since bounding_box is attached to the model here it's in reverse order.
        transform.bounding_box = ((-0.5, shape[0] - 0.5),
                                  (-0.5, shape[1] - 0.5))
    return transform


def wfss(input_model, reference_files):
    """Create the WCS pipeline for a NIRCAM grism observation.

    Parameters
    ----------
    input_model: datamodels.ImageModel
        The input datamodel, derived from datamodels
    reference_files: dict
        Reference files needed for this pipeline

    Notes
    -----
    reference_files = {
        "grism_config": 'NIRCAM_modA_R.asdf'
    }

    The tree in the grism reference file has a section for each order/beam as
    well as the link to the filter data file, not sure if there will be a
    separate passband reference file needed for the wavelength scaling or the
    wedge offsets. This file is currently created in
    jwreftools/nircam/nircam_reftools

    The direct image the catalog has been made from has been corrected for
    distortion, but the dispersed images have not. This is OK if the trace and
    dispersion solutions are defined with respect to the distortion-corrected
    image. The catalog from the combined direct image has object locations in
    in detector space. The WCS information from the direct image will be used
    to translate this to RA/DEC locations for each of the grism images. The
    grism images will then use their own WCS information to translate to
    detector space. The translations is assumed to be one-to-one for purposes
    of identifying the center of the object trace.

    The extent of the trace for each object can then be calculated based on
    the grism in use (row or column).

    The extraction box in the cross-dispersion direction will be calculated to
    be just larger that the object extent in the segmentation map associated
    with the direct image and the object catalog.

    For each spectral order, the configuration file contains a pair of
    magnitude-cutoff values. Sources with magnitudes fainter than the
    extraction cutoff (MMAG_EXTRACT_X) are will not be extracted, but are
    accounted for when computing the spectral contamination and background
    estimates.

    Sources with magnitudes fainter than the second cutoff (MMAG_MARK_X) are
    completely ignored.  Here, X equals +1, +2, etc., for each spectral order,
    as specified in the configuration file.

    The configuration file also contains the reference to the polynomial models
    which describe the trace dispersion. The sensativity information in this
    configuration file needs to be modified by the passband of the filter used
    for the direct image.


    Step 1: Convert the source catalog from the reference frame of the
            uberimage to that of the dispersed image.  For the Vanilla
            Pipeline we assume that the pointing information in the file
            headers is sufficient.  This will be strictly true if all images
            were obtained in a single visit (same guide stars).
    Step 2: Record source information for each object in the catalog: position
            (RA and Dec), shape (A_IMAGE, B_IMAGE, THETA_IMAGE), and all
            available magnitudes.
    Step 3: Compute the trace and wavelength solutions for each object in the
            catalog and for each spectral order.  Record this information.
    Step 4: Compute the WIDTH of each spectral subwindow, which may be fixed or
            variable (see discussion of optimal extraction, below).  Record
            this information.
    Step 4: Record the MMAG_EXTRACT and MMAG_MARK keywords (described below)
            for each object and spectral order.


    """

    # Where to get the association file from?

    # The input is the grism image
    if not isinstance(datamodel, ImageModel):
        raise TypeError('The input data model must be an ImageModel.')

    # Get the disperser parameters and polynomial trace model
    configuration = AsdfFile.open(reference_files['disperser']).tree

    # the wcs information for the input_model should already be
    # complete. yes? Not sure what processing has been done on the
    # input grism image. We minimally need the basic WCS information
    # here so we can translate the catalog object positions

    # catlog images are distortion corrected before catalog is made
    # but slitless images are not. The polynomials account for this.
    # In order to convert

    # Get the catalog object locations from the direct image

    # Create the model transforms.
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.deg, u.deg))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')
    distortion = imaging_distortion(input_model, reference_files)
    tel2sky = pointing.v23tosky(input_model)

    # for the grism data
    pipeline = [(objects, detector),
                (v2v3, tel2sky),
                (world, None)]


    return pipeline


def ImageXYToGrismXY(catalog, grism_config, polynomials):
    """We want to determine the x, y, and wavelength of pixel containing the
    trace originating from a pixel at x_in, y_in in the un-dispersed frame.

    We want to look at the dispersed pixels that are at location x+dxs on the
    dispersed image. Where dxs can be an array.




    """
    xy_catalog = get_catalog(catalog)



def get_catalog_info(asn_catalog_name, asn_segment_imagename):
    """
    Return the object locations from the direct image catalog
    and the segmentation map that they were created from


    In [360]: !more jw96090_t001_nircam_f322w2_mosaic_cat.ecsv
    # %ECSV 0.9
    # ---
    # datatype:
    # - {name: id, datatype: int64}
    # - {name: xcentroid, unit: pix, datatype: float64}
    # - {name: ycentroid, unit: pix, datatype: float64}
    # - {name: ra_icrs_centroid, unit: deg, datatype: float64}
    # - {name: dec_icrs_centroid, unit: deg, datatype: float64}
    # - {name: area, unit: pix2, datatype: float64}
    # - {name: source_sum, datatype: float32}
    # - {name: source_sum_err, datatype: float32}
    # - {name: semimajor_axis_sigma, unit: pix, datatype: float64}
    # - {name: semiminor_axis_sigma, unit: pix, datatype: float64}
    # - {name: orientation, unit: deg, datatype: float64}
    # - {name: orientation_sky, unit: deg, datatype: float64}
    # - {name: abmag, datatype: float64}
    # - {name: abmag_error, datatype: float32}
    id xcentroid ycentroid ra_icrs_centroid dec_icrs_centroid area source_sum source_sum_err semimajor_axis_sigma semiminor_axis_sigma orientation orientation_sky abmag abmag_error
    1 1478.79057591 8.20729790816 150.158383692 2.30276751202 77.0 20.0082 0.167517 3.25911081809 2.05169541295 -84.7150877614 185.284912239 0.0 0.00909026
    2 1579.15673474 6.91806291859 150.156709574 2.3027460732 56.0 16.1

    Do any sort of extra conversion needed here?
    Where do I get the segmentation map, can I construct the name from the
    name of the direct image in the association file?

    """
    catalog_data = Table.read(catalog, format='ascii.ecsv')
    return catalog_data






exp_type2transform = {'nrc_image': imaging,
                      'nrc_grism': wfss,
                      'nrc_tacq': imaging,
                      'nrc_taconfirm': imaging,
                      'nrc_coron': imaging,
                      'nrc_focus': imaging,
                      'nrc_tsimage': imaging,
                      'nrc_tsgrism': not_implemented_mode,
                      'nrc_led': not_implemented_mode,
                      'nrc_dark': not_implemented_mode,
                      'nrc_flat': not_implemented_mode,
                      }
