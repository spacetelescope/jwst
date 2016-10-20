from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import logging

from asdf import AsdfFile
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import Const1D, Mapping, Scale
from gwcs import wcs
import gwcs.coordinate_frames as cf
from .util import not_implemented_mode
from ..transforms.models import NirissSOSSModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def create_pipeline(input_model, reference_files):
    '''
    get reference files from crds

    '''

    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def niriss_soss_set_input(model, order_number):
    """
    Get the right model given the order number.

    Parameters
    ----------
    model - Input model
    order_number - the, well, order number desired

    Returns
    -------
    WCS - the WCS corresponding to the order_number

    """

    # Make sure the order number is correct.
    if order_number < 1 or order_number > 3:
        raise ValueError('Order must be between 1 and 3')

    # Return the correct transform based on the order_number
    obj = model.meta.wcs.forward_transform.get_model(order_number)

    # use the size of the input subarray7
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,),
                            axes_names=('wavelength',))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS(),
                            axes_names=('ra', 'dec'),
                            axes_order=(0, 1), unit=(u.deg, u.deg), name='sky')
    world = cf.CompositeFrame([sky, spec], name='world')
    pipeline = [(detector, obj),
                (world, None)
                ]

    return wcs.WCS(pipeline)


def niriss_soss(input_model, reference_files):
    """
    The NIRISS SOSS pipeline includes 3 coordinate frames -
    detector, focal plane and sky

    reference_files={'specwcs': 'soss_wavelengths_configuration.asdf'}
    """

    # Get the target RA and DEC, they will be used for setting the WCS RA and DEC based on a conversation
    # with Kevin Volk.
    try:
        target_ra = float(input_model['meta.target.ra'])
        target_dec = float(input_model['meta.target.dec'])
    except:
        # There was an error getting the target RA and DEC, so we are not going to continue.
        raise ValueError('Problem getting the TARG_RA or TARG_DEC from input model {}'.format(input_model))

    # Define the frames
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,),
                            axes_names=('wavelength',))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS(),
                            axes_names=('ra', 'dec'),
                            axes_order=(0, 1), unit=(u.deg, u.deg), name='sky')
    world = cf.CompositeFrame([sky, spec], name='world')
    try:
        with AsdfFile.open(reference_files['specwcs']) as wl:
            wl1 = wl.tree[1].copy()
            wl2 = wl.tree[2].copy()
            wl3 = wl.tree[3].copy()
    except Exception as e:
        raise IOError('Error reading wavelength correction from {}'.format(reference_files['specwcs']))

    cm_order1 = (Mapping((0, 1, 0, 1)) | (Const1D(target_ra) & Const1D(target_dec) & wl1)).rename('Order1')
    cm_order2 = (Mapping((0, 1, 0, 1)) | (Const1D(target_ra) & Const1D(target_dec) & wl2)).rename('Order2')
    cm_order3 = (Mapping((0, 1, 0, 1)) | (Const1D(target_ra) & Const1D(target_dec) & wl3)).rename('Order3')

    # Define the transforms, they should accept (x,y) and return (ra, dec, lambda)
    soss_model = NirissSOSSModel([1, 2, 3], [cm_order1, cm_order2, cm_order3]).rename('3-order SOSS Model')

    # Define the pipeline based on the frames and models above.
    pipeline = [(detector, soss_model),
                (world, None)
                ]

    return pipeline


def imaging(input_model, reference_files):
    """
    The NIRISS imaging pipeline includes 3 coordinate frames -
    detector, focal plane and sky

    reference_files={'distortion': 'jwst_niriss_distortioon_0001.asdf'}
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.arcsec, u.arcsec))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')
    distortion = imaging_distortion(input_model, reference_files)
    tel2sky = pointing.v23tosky(input_model)
    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    distortion = AsdfFile.open(reference_files['distortion']).tree['model']
    # Convert to arcsec
    transform = distortion | Scale(60) & Scale(60)
    return transform


exp_type2transform = {'nis_image': imaging,
                      'nis_wfss': not_implemented_mode,#        ?? WFSS spec
                      'nis_soss': niriss_soss,
                      'nis_ami': not_implemented_mode#        ?? imaging
                      }
