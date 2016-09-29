from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import logging
import numpy as np
from asdf import AsdfFile
from astropy.modeling import models
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits

import gwcs.coordinate_frames as cf
from gwcs import selector
from . import pointing
from ..transforms import models as jwmodels


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_pipeline(input_model, reference_files):
    '''
    get reference files from crds

    '''
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def imaging(input_model, reference_files):
    """
    The MIRI imaging pipeline includes 3 coordinate frames - detector,
    focal plane and sky

    reference_files={'distortion': 'test.asdf', 'filter_offsets': 'filter_offsets.asdf'}
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    focal = cf.Frame2D(name='focal', axes_order=(0, 1), unit=(u.arcmin, u.arcmin))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS())
    distortion = imaging_distortion(input_model, reference_files)
    #fitswcs_transform = pointing.create_fitswcs_transform(input_model)
    pipeline = [(detector, distortion),
                (focal, None)]
                #(sky, None)
                #]
    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create pixe2sky and sky2pixel transformation for the MIRI imager.

    Parameters
    ----------
    model : jwst.datamodels.ImagingModel
        input model
    reference_files : dict
        reference files from CRDS


    using CDP 3 Reference distortion file
    MIRI_FM_MIRIMAGE_F1000W_PSF_03.01.00.fits

    reference files/corrections needed (pixel to sky):

    1. Filter dependent shift in (x,y) (!with an oposite sign to that delievred by the IT)
    2. Apply MI
    3. Apply Ai and BI matrices
    4. Apply the TI matrix (this gives V2/V3 coordinates)
    5. Apply V2/V3 to sky transformation

    ref_file: filter_offset.asdf - (1)
    ref_file: distortion.asdf -(2,3,4)
    """
    distortion = AsdfFile.open(reference_files['distortion']).tree['model']
    filter_offset = AsdfFile.open(reference_files['filteroffset']).tree[input_model.meta.instrument.filter]
    full_distortion = models.Shift(filter_offset['row_offset']) & models.Shift(
        filter_offset['column_offset']) | distortion
    full_distortion = full_distortion.rename('distortion')
    return full_distortion


def lrs(input_model, reference_files):
    """
    Create the WCS pipeline for a MIRI fixed slit observation.

    reference_files = {"specwcs": 'MIRI_FM_MIRIMAGE_P750L_DISTORTION_04.02.00.fits'}
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    focal_spatial = cf.Frame2D(name='focal', axes_order=(0, 1), unit=(u.arcmin, u.arcmin))
    #sky = cf.CelestialFrame(reference_frame=coord.ICRS())
    spec = cf.SpectralFrame(name='wavelength', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    focal = cf.CompositeFrame([focal_spatial, spec])

    ref = fits.open(reference_files['specwcs'])
    ldata = ref[1].data
    if input_model.meta.exposure.type.lower() == 'mir_lrs-fixedslit':
        zero_point = ref[1].header['imx'], ref[1].header['imy']
    elif input_model.meta.exposure.type.lower() == 'mir_lrs-slitless':
        #zero_point = ref[1].header['imysltl'], ref[1].header['imxsltl']
        #zero point in reference file is wrong
        # This should eb moved eventually to the reference file.
        zero_point = [35, 442]#[35, 763] # account for subarray
    lrsdata = np.array([l for l in ldata])
    x0 = lrsdata[:, 3]
    x1 = lrsdata[:, 5]
    y0 = lrsdata[:, 4]
    domain = [{'lower': x0.min() + zero_point[0], 'upper': x1.max() + zero_point[0]},
              {'lower': (y0.min() + zero_point[1]), 'upper': (y0.max() + zero_point[1])}
              ]
    log.info("Setting domain to {0}".format(domain))
    lrs_wav_model = jwmodels.LRSWavelength(lrsdata, zero_point)
    ref.close()
    angle = np.arctan(0.00421924)
    spatial = models.Rotation2D(angle)
    det2focal = models.Mapping((0, 1, 0, 1)) | spatial & lrs_wav_model
    det2focal.meta['domain'] = domain
    pipeline = [(detector, det2focal),
                (focal, None)]
                #(sky, None)
                #]
    return pipeline


def ifu(input_model, reference_files):
    """
    Create the WCS pipeline for a MIRI IFU observation.
    """

    #reference_files = {'distortion': 'jwst_miri_distortion_00001.asdf', #files must hold 2 channels each
                        #'specwcs': 'jwst_miri_specwcs_00001.asdf',
                        #'regions': 'jwst_miri_regions_00001.asdf',
                        #'v2v3': 'jwst_miri_v2v3_00001.asdf'
                        #'wavelengthrange': 'jwst_miri_wavelengthrange_0001.asdf'}
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    alpha_beta = cf.Frame2D(name='alpha_beta_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=('alpha', 'beta'))
    spec_local = cf.SpectralFrame(name='alpha_beta_spectral', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    miri_focal = cf.CompositeFrame([alpha_beta, spec_local], name='alpha_beta')
    xyan_spatial = cf.Frame2D(name='Xan_Yan_spatial', axes_order=(0, 1), unit=(u.arcmin, u.arcmin), axes_names=('v2', 'v3'))
    spec = cf.SpectralFrame(name='Xan_Yan_spectral', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    xyan = cf.CompositeFrame([xyan_spatial, spec], name='Xan_Yan')
    v23_spatial = cf.Frame2D(name='V2_V3_spatial', axes_order=(0, 1), unit=(u.arcmin, u.arcmin), axes_names=('v2', 'v3'))
    spec = cf.SpectralFrame(name='V2_v3_spectral', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    v2v3 = cf.CompositeFrame([v23_spatial, spec], name='V2_V3')
    icrs = cf.CelestialFrame(name='icrs', reference_frame=coord.ICRS(),
                             axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=('RA', 'DEC'))
    sky = cf.CompositeFrame([icrs, spec], name='sky_and_wavelength')
    det2alpha_beta = (detector_to_alpha_beta(input_model, reference_files)).rename(
        "detector_to_alpha_beta")
    ab2xyan = (alpha_beta2XanYan(input_model, reference_files)).rename("alpha_beta_to_Xan_Yan")
    xyan2v23 = models.Identity(1) & (models.Shift(7.8) | models.Scale(-1)) & models.Identity(1)
    #fitswcs_transform = pointing.create_fitswcs_transform(input_model) & models.Identity(1)
    pipeline = [(detector, det2alpha_beta),
                (miri_focal, ab2xyan),
                (xyan, xyan2v23),
                (v2v3, None)]
                 #fitswcs_transform),
                #(sky, None)]
    return pipeline


def detector_to_alpha_beta(input_model, reference_files):
    """
    Create the transform from detector to alpha, beta frame.

    forward transform:
      RegionsSelector
        label_mapper is the regions array
        selector is {slice_number: alphs_model & beta_model & lambda_model}
    backward transform:
      RegionsSelector
        label_mapper is LabelMapperDict
           {channel_wave_range (): LabelMapperDict}
                                   {beta: slice_number}
        selector is {slice_number: x_transform & y_transform}
    """
    band = input_model.meta.instrument.band
    channel = input_model.meta.instrument.channel
    # used to read the wavelength range
    channels = [c + band for c in channel]

    f = AsdfFile.open(reference_files['distortion'])
    # The name of the model indicates the output coordinate
    alpha_model = f.tree['alpha_model']
    beta_model = f.tree['beta_model']
    x_model = f.tree['x_model']
    y_model = f.tree['y_model']
    slice_model = f.tree['slice_model']
    f.close()
    f = AsdfFile.open(reference_files['specwcs'])
    lambda_model = f.tree['model']
    f.close()
    f = AsdfFile.open(reference_files['regions'])
    regions = f.tree['regions'].copy()
    f.close()
    label_mapper = selector.LabelMapperArray(regions)
    transforms = {}

    for sl in alpha_model:
        #chan = str(sl // 100) + band
        forward = models.Mapping([1, 0, 0, 1, 0]) | \
                alpha_model[sl] & beta_model[sl] & lambda_model[sl]
        inv = models.Mapping([2, 0, 2, 0]) | x_model[sl] & y_model[sl]
        forward.inverse = inv
        transforms[sl] = forward

    f = AsdfFile.open(reference_files['wavelengthrange'])
    # the following should go in the asdf reader
    wave_range = f.tree['wavelengthrange'].copy()
    wave_channels = f.tree['channels']
    wr = {}
    for ch, r in zip(wave_channels, wave_range):
        wr[ch] = r
    f.close()
    ch_dict = {}
    for c in channels:
        ch_dict.update({tuple(wr[c]): selector.LabelMapperDict(('alpha', 'beta', 'lam'), slice_model[c],
                                                   models.Mapping([1, ], n_inputs=3))})
    alpha_beta_mapper = selector.LabelMapperRange(('alpha', 'beta', 'lam'), ch_dict,
                                                  models.Mapping((2,)))
    label_mapper.inverse = alpha_beta_mapper


    det2alpha_beta = selector.RegionsSelector(('x', 'y'), ('alpha', 'beta', 'lam'),
                                              label_mapper=label_mapper, selector=transforms)
    return det2alpha_beta


def alpha_beta2XanYan(input_model, reference_files):
    """
    Create the transform from detector to Xan, Yan frame.

    forward transform:
      RegionsSelector
        label_mapper is LabelMapperDict()
        {channel_wave_range (): channel_number}
        selector is {channel_number: ab2Xan & ab2Yan}
    bacward_transform
      RegionsSelector
        label_mapper is LabelMapperDict()
        {channel_wave_range (): channel_number}
        selector is {channel_number: Xan2ab & Yan2ab}
    """
    band = input_model.meta.instrument.band
    channel = input_model.meta.instrument.channel
    # used to read the wavelength range
    channels = [c + band for c in channel]

    f = AsdfFile.open(reference_files['v2v3'])
    v23 = f.tree['model']
    f.close()
    f = AsdfFile.open(reference_files['wavelengthrange'])
    # the following should go in the asdf reader
    wave_range = f.tree['wavelengthrange'].copy()
    wave_channels = f.tree['channels']
    wr = {}
    for ch, r in zip(wave_channels, wave_range):
        wr[ch] = r
    f.close()

    dict_mapper = {}
    sel = {}
    for c in channels:
        ch = int(c[0])
        dict_mapper[tuple(wr[c])] = models.Mapping((2,), name="mapping_lam") | \
                   models.Const1D(ch, name="channel #")
        map1 = models.Mapping((1, 0, 1, 0), name='map2poly')
        map1._outputs = ('alpha', 'beta', 'alpha', 'beta')
        map1._inputs = ('alpha', 'beta')
        map1.inverse = models.Mapping((0, 1))
        ident1 = models.Identity(1, name='identity_lam')
        ident1._inputs = ('lam',)
        chan_v23 = v23[c]
        v23chan_backward = chan_v23.inverse
        del chan_v23.inverse
        v23_spatial = map1 | chan_v23
        v23_spatial.inverse = map1 | v23chan_backward
        v23c = v23_spatial & ident1
        sel[ch] = v23c

    wave_range_mapper = selector.LabelMapperRange(('alpha', 'beta', 'lam'), dict_mapper,
                                                  inputs_mapping=models.Mapping([2, ]))
    wave_range_mapper.inverse = wave_range_mapper.copy()
    ab2xyan = selector.RegionsSelector(('alpha', 'beta', 'lam'), ('v2', 'v3', 'lam'),
                                      label_mapper=wave_range_mapper,
                                      selector=sel)

    return ab2xyan

exp_type2transform = {'mir_image': imaging,
                      'mir_tacq': imaging,
                      'mir_lyot': imaging,
                      'mir_4qpm': imaging,
                      'mir_coroncal': imaging,
                      'mir_lrs-fixedslit': lrs,
                      'mir_lrs-slitless': lrs,
                      'mir_mrs': ifu
                      }
