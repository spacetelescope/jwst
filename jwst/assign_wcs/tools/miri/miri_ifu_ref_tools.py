# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools to create WCS reference files for MIRI MRS using IDT reference
products delivered with CDP-4:

MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_04.02.00.fits

MIRI MRS uses 5 reference files of type:

DISTORTION
SPECWCS
REGIONS
WAVELENGTHRANGE
V2V3

create_cdp4_references() creates all reference files.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from asdf import AsdfFile

from astropy.io import fits
from astropy.modeling import models

from matplotlib import image as mplimage
from matplotlib import pyplot as plt


def create_cdp5_references(fname, ref):
    """
    Create ASDF WCS reference files for MIRI MRS data from a CDP5 reference file.

    Parameters
    ----------
    fname : str
        name of reference file
    ref : dict
        A dictionary {reftype: refname}, e.g.
        {'distortion': 'jwst_miri_distortion_0001.asdf',
         'regions': 'jwst_miri_regioons_0001.asdf',
         'specwcs': 'jwst_miri_specwcs_0001.asdf',
         'v2v3': 'jwst_miri_v2v3_00001.asdf',
         'wavelengthrange': 'jwst_miri_wavelengthrange_0001.asdf'}
         }

    Examples
    --------
    >>> create_cdp4_references('MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_04.02.00.fits', ref)

    """
    with fits.open(fname) as f:
        channel = f[0].header['CHANNEL']
        band = f[0].header['BAND']
        detector = f[0].header['DETECTOR']
        ch1 = 'CH{0}'.format(channel[0])
        ch2 = 'CH{0}'.format(channel[1])
        slices = f[1].data
        fov1 = f[2].data
        fov2 = f[3].data
        alpha1 = f[('Alpha_' + ch1, 1)].data
        lam1 = f[('Lambda_' + ch1, 1)].data
        alpha2 = f[('Alpha_' + ch2, 1)].data
        lam2 = f[('Lambda_' + ch2, 1)].data
        x1 = f[('X_' + ch1, 1)].data
        y1 = f[('Y_' + ch1, 1)].data
        x2 = f[('X_' + ch2, 1)].data
        y2 = f[('Y_' + ch2, 1)].data
        ab_v23 = f[('albe_to_XANYAN', 1)].data.copy()
        v23_ab = f[('XANYAN_to_albe', 1)].data.copy()
        b0_ch1 = f[0].header['B_ZERO' + ch1[2]]
        bdel_ch1 = f[0].header['B_DEL' + ch1[2]]
        b0_ch2 = f[0].header['B_ZERO' + ch2[2]]
        bdel_ch2 = f[0].header['B_DEL' + ch2[2]]
    # Get channel names, e.g. 1LONG, 2LONG
    channels = [c + band for c in channel]

    coeff_names = build_coeff_names_cdp5(alpha1.names)
    amodel1 = create_poly_models(alpha1, int(channel[0]), coeff_names, name='det2local')
    lmodel1 = create_poly_models(lam1, int(channel[0]), coeff_names, name='det2local')
    amodel2 = create_poly_models(alpha2, int(channel[1]), coeff_names, name='det2local')
    lmodel2 = create_poly_models(lam2, int(channel[1]), coeff_names, name='det2local')
    # reverse models

    # 'x' in the report corresponds to y in python and 'y' to x,
    # The x/y models take (lam, alpha)
    xmodel1 = create_xy_models(x1, int(channel[0]), coeff_names, name='x')
    ymodel1 = create_xy_models(y1, int(channel[0]), coeff_names, name='y')
    xmodel2 = create_xy_models(x2, int(channel[1]), coeff_names, name='x')
    ymodel2 = create_xy_models(y2, int(channel[1]), coeff_names, name='y')
    amodel1.update(amodel2)
    xmodel1.update(xmodel2)
    ymodel1.update(ymodel2)
    lmodel1.update(lmodel2)
    bmodel1, smodel1 = create_beta_models(b0_ch1, bdel_ch1, int(channel[0]), len(alpha1))
    bmodel2, smodel2 = create_beta_models(b0_ch2, bdel_ch2, int(channel[1]), len(alpha2))
    bmodel1.update(bmodel2)

    create_distortion_file('DISTORTION', detector, band, channel,
                           (amodel1, bmodel1, xmodel1, ymodel1, smodel1, smodel2), ref['distortion'])

    create_specwcs_file('SPECWCS', detector, band, channel, lmodel1, ref['specwcs'])

    create_v23('V2V3', detector, band, channels, (ab_v23, v23_ab), ref['v2v3'])

    create_regions_file(slices, detector, band, channel, ref['regions'])

    create_wavelengthrange_file(ref['wavelengthrange'])


def create_cdp4_references(fname, ref):
    """
    Create ASDF WCS reference files for MIRI MRS data from a CDP4 reference file.

    Parameters
    ----------
    fname : str
        name of reference file
    ref : dict
        A dictionary {reftype: refname}, e.g.
        {'distortion': 'jwst_miri_distortion_0001.asdf',
         'regions': 'jwst_miri_regioons_0001.asdf',
         'specwcs': 'jwst_miri_specwcs_0001.asdf',
         'v2v3': 'jwst_miri_v2v3_00001.asdf',
         'wavelengthrange': 'jwst_miri_wavelengthrange_0001.asdf'}

         }

    Examples
    --------
    >>> create_cdp4_references('MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_04.02.00.fits', ref)

    """
    with fits.open(fname) as f:
        channel = f[0].header['CHANNEL']
        band = f[0].header['BAND']
        detector = f[0].header['DETECTOR']
        ch1 = 'CH{0}'.format(channel[0])
        ch2 = 'CH{0}'.format(channel[1])
        slices = f[1].data
        fov1 = f[2].data
        fov2 = f[3].data
        alpha1 = f[('Alpha-' + ch1, 1)].data
        lam1 = f[('Lambda-' + ch1, 1)].data
        alpha2 = f[('Alpha-' + ch2, 1)].data
        lam2 = f[('Lambda-' + ch2, 1)].data
        x1 = f[('X-' + ch1, 1)].data
        y1 = f[('Y-' + ch1, 1)].data
        x2 = f[('X-' + ch2, 1)].data
        y2 = f[('Y-' + ch2, 1)].data
        ab_v23 = f[('al,be->V2/V3', 1)].data.copy()
        v23_ab = f[('V2/V3->al,be', 1)].data.copy()
        b0_ch1 = f[0].header['B_ZERO' + ch1[2]]
        bdel_ch1 = f[0].header['B_DEL' + ch1[2]]
        b0_ch2 = f[0].header['B_ZERO' + ch2[2]]
        bdel_ch2 = f[0].header['B_DEL' + ch2[2]]
    # Get channel names, e.g. 1LONG, 2LONG
    channels = [c + band for c in channel]

    coeff_names = build_coeff_names(alpha1.names)
    amodel1 = create_poly_models(alpha1, int(channel[0]), coeff_names, name='det2local')
    lmodel1 = create_poly_models(lam1, int(channel[0]), coeff_names, name='det2local')
    amodel2 = create_poly_models(alpha2, int(channel[1]), coeff_names, name='det2local')
    lmodel2 = create_poly_models(lam2, int(channel[1]), coeff_names, name='det2local')
    # reverse models

    # 'x' in the report corresponds to y in python and 'y' to x,
    # The x/y models take (lam, alpha)
    xmodel1 = create_xy_models(x1, int(channel[0]), coeff_names, name='y')
    ymodel1 = create_xy_models(y1, int(channel[0]), coeff_names, name='x')
    xmodel2 = create_xy_models(x2, int(channel[1]), coeff_names, name='y')
    ymodel2 = create_xy_models(y2, int(channel[1]), coeff_names, name='x')
    amodel1.update(amodel2)
    xmodel1.update(xmodel2)
    ymodel1.update(ymodel2)
    lmodel1.update(lmodel2)
    bmodel1, smodel1 = create_beta_models(b0_ch1, bdel_ch1, int(channel[0]), len(alpha1))
    bmodel2, smodel2 = create_beta_models(b0_ch2, bdel_ch2, int(channel[1]), len(alpha2))
    bmodel1.update(bmodel2)

    create_distortion_file('DISTORTION', detector, band, channel,
                           (amodel1, bmodel1, xmodel1, ymodel1, smodel1, smodel2), ref['distortion'])

    create_specwcs_file('SPECWCS', detector, band, channel, lmodel1, ref['specwcs'])

    create_v23('V2V3', detector, band, channels, (ab_v23, v23_ab), ref['v2v3'])

    create_regions_file(slices, detector, band, channel, ref['regions'])

    create_wavelengthrange_file(ref['wavelengthrange'])


def create_regions_file(slices, detector, band, channel, name):
    tree = create_reffile_header("REGIONS", detector, band, channel)
    f = AsdfFile()
    tree['regions'] = slices
    f.tree = tree
    f.write_to(name)


def create_reffile_header(reftype, detector, band, channel):
    tree = {"detector": detector,
            "instrument": "MIRI",
            "band": band,
            "channel": channel,
            "pedigree": "GROUND",
            "title": "MIRI IFU model - based on CDP-4",
            "reftype": reftype,
            "author": "N. Dencheva",
            "exp_type": "MIR_MRS"
            }
    return tree


def create_v23(reftype, detector, band, channels, data, name):
    """
    Create the transform from MIRI Local to telescope V2/V3 system for all channels.
    """
    channel = "".join([ch[0] for ch in channels])
    tree = {"detector": detector,
            "instrument": "MIRI",
            "band": band,
            "channel": channel,
            "exp_type": "MIR_MRS",
            "pedigree": "GROUND",
            "title": "MIRI IFU model - based on CDP-4",
            "reftype": reftype,
            "author": "N. Dencheva"
            }
    ab_v23 = data[0]
    v23_ab = data[1]
    m = {}
    c0_0, c0_1, c1_0, c1_1 = ab_v23[0][1:]
    ch1_v2 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    c0_0, c0_1, c1_0, c1_1 = v23_ab[0][1:]
    ch1_a = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")

    c0_0, c0_1, c1_0, c1_1 = ab_v23[1][1:]
    ch1_v3 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    c0_0, c0_1, c1_0, c1_1 = v23_ab[1][1:]
    ch1_b = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")
    c0_0, c0_1, c1_0, c1_1 = ab_v23[2][1:]
    ch2_v2 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    c0_0, c0_1, c1_0, c1_1 = v23_ab[2][1:]
    ch2_a = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")

    c0_0, c0_1, c1_0, c1_1 = ab_v23[3][1:]
    ch2_v3 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    c0_0, c0_1, c1_0, c1_1 = v23_ab[3][1:]
    ch2_b = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")
    ch1_for = ch1_v2 & ch1_v3
    ch2_for = ch2_v2 & ch2_v3
    ch1_for.inverse = ch1_a & ch1_b
    ch2_for.inverse = ch2_a & ch2_b
    m[channels[0]] = ch1_for
    m[channels[1]] = ch2_for
    tree['model'] = m

    f = AsdfFile()
    f.tree = tree
    f.write_to(name)


def create_distortion_file(reftype, detector, band, channel, data, name):

    tree = create_reffile_header(reftype, detector, band, channel)

    adata, bdata, xdata, ydata, sdata1, sdata2 = data
    tree['alpha_model'] = adata
    tree['beta_model'] = bdata
    tree['x_model'] = xdata
    tree['y_model'] = ydata
    tree['slice_model'] = {str(channel[0]) + band: sdata1, str(channel[1]) + band: sdata2}
    f = AsdfFile()
    f.tree = tree
    f.write_to(name)


def create_specwcs_file(reftype, detector, band, channel, lmodel, name):
    tree = create_reffile_header(reftype, detector, band, channel)
    tree['model'] = lmodel
    f = AsdfFile()
    f.tree = tree
    f.write_to(name)


def create_poly_models(data, channel, coeff_names, name):
    """
    Create a 2D polynomial model for the transformation
    detector --> local MIRI frame
    Works for alpha and lambda coordinates.
    """
    nslices = len(data)
    sl = channel * 100 + np.arange(1, nslices + 1)

    transforms = {}
    for i in range(nslices):
        sl = channel * 100 + i + 1
        al = data[i]
        xs = al[0]
        coeffs = {}
        for c, val in zip(coeff_names, al[1:]):
            coeffs[c] = val
        transforms[sl] = models.Identity(1) & models.Shift(-xs) | \
                  models.Polynomial2D(8, name=name, **coeffs)
    return transforms


def create_xy_models(data, channel, coeff_names, name):
    """
    Create a 2D polynomial model for the transformation
    local_MIRI --> detector frame.
    """
    nslices = len(data)
    sl = channel * 100 + np.arange(1, nslices + 1)
    shname = "shift_{0}".format(name)
    pname = "polynomial_{0}".format(name)
    transforms = {}
    for i in range(nslices):
        sl = channel * 100 + i + 1
        al = data[i]
        xs = al[0]
        coeffs = {}
        for c, val in zip(coeff_names, al[1:]):
            coeffs[c] = val

        transforms[sl] =  models.Shift(-xs, name=shname) & models.Identity(1) | \
                  models.Polynomial2D(8, name=pname, **coeffs)
    return transforms


def build_coeff_names_cdp5(names):
    names = names[1:]
    names = [name.replace('VAR2_', "c") for name in names]
    return names


def build_coeff_names(names):
    names = names[1:]
    names = [name.replace('VAR2(', "c") for name in names]
    names = [name.replace(')', "") for name in names]
    names = [name.replace(',', "_") for name in names]
    return names


def create_beta_models(b0, bdel, channel, nslices):
    beta = {}
    slices = {}
    for s in range(nslices):
        sl = channel * 100 + s + 1
        beta_s = b0 + s * bdel
        m = models.Const1D(beta_s, name='det2local') #xy2beta and xy2lam
        beta[sl] = m
        inv = models.Const1D(sl)
        slices[beta_s] = models.Mapping([1, ]) | inv
    return beta, slices


def create_wavelengthrange_file(name):
    f = AsdfFile()
    #wavelengthrange = {'1SHORT': (4.88, 5.77),
                        #'1MEDIUM': (5.64, 6.67),
                        #'1LONG': (6.50, 7.70),
                        #'2SHORT': (7.47, 8.83),
                        #'2MEDIUM': (8.63, 10.19),
                        #'2LONG': (9.96, 11.77),
                        #'3SHORT': (11.49, 13.55),
                        #'3MEDIUM': (13.28, 15.66),
                        #'3LONG': (15.34, 18.09),
                        #'4SHORT': (17.60, 21.00),
                        #'4MEDIUM': (20.51, 24.48),
                        #'4LONG': (23.92, 28.55)
                        #}
    # Relaxing the range to match the distortion. The table above
    # comes from the report and is "as designed".
    wavelengthrange = {'1SHORT': (4.68, 5.97),
                        '1MEDIUM': (5.24, 6.87),
                        '1LONG': (6.2, 7.90),
                        '2SHORT': (7.27, 9.03),
                        '2MEDIUM': (8.43, 10.39),
                        '2LONG': (9.76, 11.97),
                        '3SHORT': (11.29, 13.75),
                        '3MEDIUM': (13.08, 15.86),
                        '3LONG': (15.14, 18.29),
                        '4SHORT': (17.40, 21.20),
                        '4MEDIUM': (20.31, 24.68),
                        '4LONG': (23.72, 28.75)
                        }
    channels = ['1SHORT', '1MEDIUM', '1LONG', '2SHORT', '2MEDIUM', '2LONG',
                '3SHORT', '3MEDIUM', '3LONG', '4SHORT', '4MEDIUM', '4LONG']
    tree = {
            "instrument": "MIRI",
            "exp_type": "MIR_MRS",
            "pedigree": "GROUND",
            "title": "MIRI IFU model - based on CDP-4",
            "reftype": "WAVELENGTHRANGE",
            "author": "N. Dencheva"
            }
    tree['channels'] = channels
    f.tree = tree
    vr = np.empty((12, 2), dtype=np.float)
    for i, ch in enumerate(channels):
        vr[i] = wavelengthrange[ch]
    f.tree['wavelengthrange'] = vr
    f.write_to(name)
