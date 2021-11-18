# TODO this file provides plots for development purposes and should be removed along with the devname argument and associated calls when ready.

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.visualization import MinMaxInterval, LogStretch, ImageNormalize


def plot_weights(scidata, scimask=None):

    if scimask is None:
        norm = ImageNormalize(scidata, interval=MinMaxInterval())
    else:
        norm = ImageNormalize(scidata[~scimask], interval=MinMaxInterval())

    plt.figure(figsize=(13, 5))

    gs = gridspec.GridSpec(2, 1, height_ratios=[40, 1])

    plt.subplot(gs[0, 0], aspect=2)
    im = plt.pcolormesh(scidata, norm=norm)

    cax = plt.subplot(gs[1, 0])
    plt.colorbar(im, cax=cax, orientation='horizontal')

    plt.tight_layout()
    plt.show()
    plt.close()

    return


def plot_data(scidata, scimask):

    norm = ImageNormalize(scidata[~scimask], interval=MinMaxInterval(),
                          stretch=LogStretch())

    plt.figure(figsize=(13, 5))

    gs = gridspec.GridSpec(2, 1, height_ratios=[40, 1])

    plt.subplot(gs[0, 0], aspect=2)
    im = plt.pcolormesh(scidata, norm=norm)

    cax = plt.subplot(gs[1, 0])
    plt.colorbar(im, cax=cax, orientation='horizontal')

    plt.tight_layout()
    plt.show()
    plt.close()

    return


def plot_delta(delta, scimask, vlim=None):

    if vlim is None:
        vlim = np.nanmax(np.abs(delta[~scimask]))

    plt.figure(figsize=(13, 5))

    gs = gridspec.GridSpec(2, 1, height_ratios=[40, 1])

    plt.subplot(gs[0, 0], aspect=2)
    im = plt.pcolormesh(delta, vmin=-vlim, vmax=vlim, cmap=plt.cm.coolwarm)

    cax = plt.subplot(gs[1, 0])
    plt.colorbar(im, cax=cax, orientation='horizontal')

    plt.tight_layout()
    plt.show()
    plt.close()

    return


def diagnostic_plot(scidata, scierr, scimask, model_order1, model_order2, devname):
    """Make some simple diagnostic plots."""

    # Plot the error and mask.
    norm = ImageNormalize(scierr[~scimask], interval=MinMaxInterval(),
                          stretch=LogStretch())

    plt.figure(figsize=(13, 8))

    gs = gridspec.GridSpec(4, 1, height_ratios=[40, 1, 40, 1])

    plt.subplot(gs[0, 0], aspect=2)
    im = plt.pcolormesh(scierr, norm=norm)

    cax = plt.subplot(gs[1, 0])
    plt.colorbar(im, cax=cax, orientation='horizontal')

    plt.subplot(gs[2, 0], aspect=2)
    plt.pcolormesh(scimask)

    plt.tight_layout()
    plt.savefig(devname + '_error_and_mask.png')
    plt.close()

    # Plot the data and model image.
    norm = ImageNormalize(scidata[~scimask], interval=MinMaxInterval(),
                          stretch=LogStretch())

    plt.figure(figsize=(13, 8))

    gs = gridspec.GridSpec(4, 1, height_ratios=[40, 1, 40, 1])

    plt.subplot(gs[0, 0], aspect=2)
    im = plt.pcolormesh(scidata, norm=norm)

    cax = plt.subplot(gs[1, 0])
    plt.colorbar(im, cax=cax, orientation='horizontal')

    plt.subplot(gs[2, 0], aspect=2)
    im = plt.pcolormesh(model_order1 + model_order2, norm=norm)

    cax = plt.subplot(gs[3, 0])
    plt.colorbar(im, cax=cax, orientation='horizontal')

    plt.tight_layout()
    plt.savefig(devname + '_data_and_model.png')
    plt.close()

    # Plot the absolute and relative residuals.
    residual = scidata - model_order1 - model_order2

    vlim1 = np.nanmax(np.abs(residual[~scimask]))
    with np.errstate(divide='ignore'):
        vlim2 = np.nanmax(np.abs((residual/scierr)[~scimask]))

    plt.figure(figsize=(13, 8))

    gs = gridspec.GridSpec(4, 1, height_ratios=[40, 1, 40, 1])

    plt.subplot(gs[0, 0], aspect=2)
    im = plt.pcolormesh(residual, vmin=-vlim1, vmax=vlim1, cmap=plt.cm.coolwarm)

    cax = plt.subplot(gs[1, 0])
    plt.colorbar(im, cax=cax, orientation='horizontal')

    plt.subplot(gs[2, 0], aspect=2)
    with np.errstate(divide='ignore'):
        vlim2 = 5
        im = plt.pcolormesh(residual/scierr, vmin=-vlim2, vmax=vlim2, cmap=plt.cm.coolwarm)

    cax = plt.subplot(gs[3, 0])
    plt.colorbar(im, cax=cax, orientation='horizontal')

    plt.tight_layout()
    plt.savefig(devname + '_residuals.png')
    plt.close()

    return


def plot_1d_spectra(wavelengths, fluxes, fluxerrs, npixels, devname):

    plt.figure(figsize=(13, 8))

    plt.subplot(311)

    for order in wavelengths.keys():
        plt.plot(wavelengths[order], fluxes[order])

    plt.xlabel('Wavelength [microns]')
    plt.ylabel('Flux')

    plt.subplot(312)

    for order in wavelengths.keys():
        plt.plot(wavelengths[order], fluxerrs[order])

    plt.xlabel('Wavelength [microns]')
    plt.ylabel('Flux Error')

    plt.subplot(313)

    for order in wavelengths.keys():
        plt.plot(wavelengths[order], npixels[order])

    plt.xlabel('Wavelength [microns]')
    plt.ylabel('# Pixels')

    plt.tight_layout()
    plt.savefig(devname + '_1d_spectra.png')
    plt.close()

    return


def plot_trace(scidata, scimask, xdat, ydat, xrot, yrot):

    norm = ImageNormalize(scidata[~scimask], interval=MinMaxInterval(),
                          stretch=LogStretch())

    plt.figure(figsize=(13, 5))

    plt.subplot(111, aspect=2)

    plt.pcolormesh(scidata, norm=norm)
    plt.plot(xdat, ydat)
    plt.plot(xrot, yrot)

    plt.xlim(0, 2048)
    plt.ylim(0, 256)

    plt.tight_layout()
    plt.show()
    plt.close()

    return
