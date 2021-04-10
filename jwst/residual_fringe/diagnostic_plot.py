import os
#
import numpy as np
import math
#
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#
from astropy.io import fits
from astropy.table import Table
import astropy.units as u


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def diagnostic_plot(type, mode, plot_name, xdata=None, ydata=None, idata=None, save_data=False):
    """Function to take processing data and produce diagnostic plot

    :Parameters:

    type: str, required
        the plot type. Options are: periodogram - the results of the periodogram
                                    mask1d - shows the trace or column used
                                    fits1d - the fits to the background and residual fringes

    mode: str, required
        the fitting mode of ResidualFringeCorrection, either 'columns' or 'iso-alpha'

    plot_name: str, required
        the full path of the output plot

    xdata: list, optional
        list of x arrays if plotting graphs

    ydata: list, optional
        list of y arrays if plotting graphs

    idata: list, optional
        list of images if plotting images

    save_data: boolean, optional
        save the data to a npy file, default=False

    :Returns:

    None

    """
    log.debug("diagnostic_plot: making {} plot".format(type))
    if type == 'fit_quality':
        fig, axs = plt.subplots(1, 1, figsize=(8, 4), sharex=True)

        axs.plot(xdata[0], ydata[0], c='r', markersize=0, linestyle='-', linewidth=1, label='rfc signal')
        axs.plot(xdata[0], ydata[1], c='b', markersize=0, linestyle='-',
                    linewidth=1, label='sine fit')

        axs.grid()
        axs.legend(prop={'size': 8}, loc=0)
        axs.set_ylim(-0.2, 0.2)

        plt.tight_layout(h_pad=0)
        fig.savefig(plot_name, format='pdf', dpi=150)
        fig.clear()
        plt.close(fig)
        del fig
        log.debug(" fit quality pdf saved to {}".format(plot_name))

        if save_data:
            log.debug(" saving fit_quality data to .npy file")
            out_data_name = os.path.splitext(plot_name)[0] + '_data.npy'
            np.save(out_data_name, np.array([xdata[0], ydata[0], ydata[1]]))
            log.debug(" fit quality .npy file saved to {}".format(out_data_name))


    elif type == 'periodogram':
        fig, axs = plt.subplots(2, 1, figsize=(14, 9))
        axs[0].plot(xdata[0], ydata[0], c='b', markersize=0, linestyle='-', linewidth=1)
        axs[0].set_xlabel(r'Wavenumber (cm$^{-1})$ / constant')
        axs[0].set_ylabel(r'Residual fringes')
        axs[0].set_ylim(-0.5, 0.5)
        axs[0].grid()
        axs[1].plot(xdata[1], ydata[1])
        axs[1].plot(xdata[2], ydata[2], "x")
        axs[1].set_xlabel(r'Fringe frequency (cm$^{-1}$) / constant')
        axs[1].set_ylabel(r'Power')
        axs[1].grid()
        plt.tight_layout(h_pad=0)
        fig.savefig(plot_name, format='pdf', dpi=150)
        fig.clear()
        plt.close(fig)
        del fig
        log.debug(" periodogram pdf saved to {}".format(plot_name))

        if save_data:
            log.debug(" saving peridogram data to .npy file")
            out_data_name = os.path.splitext(plot_name)[0] + '_data.npy'
            np.save(out_data_name, ydata[1])
            log.debug(" fits1d .npy file saved to {}".format(out_data_name))



    elif type == 'mask1d':
        if mode == 'columns':
            fig, axs = plt.subplots(1, 1, figsize=(6, 6))
            axs.imshow(idata[0], origin='lower', interpolation='nearest',
                       cmap='gray')
            axs.plot(xdata[0], ydata[0], color='b', markersize=0, linewidth=1, alpha=0.8)
            fig.savefig(plot_name, format='pdf', dpi=150)
            fig.clear()
            plt.close(fig)
            del fig
            log.debug(" mask1d pdf saved to {}".format(plot_name))

        if mode == 'iso-alpha':
            fig, axs = plt.subplots(1, 1, figsize=(6, 6))
            axs.imshow(idata[0], origin='lower', interpolation='nearest',
                       cmap='gray')
            axs.imshow(idata[1], origin='lower', interpolation='nearest', alpha=0.8, cmap='inferno')
            fig.savefig(plot_name, format='pdf', dpi=150)
            fig.clear()
            plt.close(fig)
            del fig
            log.debug(" mask1d pdf saved to {}".format(plot_name))

    elif type == 'fits1d':

        fig, axs = plt.subplots(4, 1, figsize=(8, 6), sharex=True)

        axs[0].plot(xdata[0], ydata[0], c='r', markersize=0, linestyle='-', linewidth=1, label='signal')
        axs[0].plot(xdata[0], ydata[1], c='g', markersize=0, linestyle='-', linewidth=1, label='baseline fit')
        axs[1].plot(xdata[0], ydata[2] * (ydata[3] > 1e-03).astype(int), c='b', markersize=0, linestyle='-',
                    linewidth=1, label='residual fringes')
        axs[1].plot(xdata[0], ydata[4] * (ydata[3] > 1e-03).astype(int), c='orange', markersize=0, linestyle='-',
                    linewidth=1, label='residual fringes fit')
        axs[2].plot(xdata[0], ydata[6] + np.mean(ydata[6]) * 0.75, c='b', markersize=0,
                    linestyle='-', linewidth=1, label='residual fringe corrected + offset')
        axs[2].plot(xdata[0], ydata[0], c='r', markersize=0, linestyle='-', linewidth=1, label='signal')
        axs[3].plot(xdata[0], ydata[3], c='b', markersize=0, linestyle='-', linewidth=1, label='weights')
        axs[3].plot(xdata[0], ydata[5], c='r', markersize=0, linestyle='--', linewidth=1,
                    label='weights with features')

        axs[0].grid()
        axs[1].grid()
        axs[2].grid()
        axs[3].grid()

        axs[0].legend(prop={'size': 8}, loc=0)
        axs[1].legend(prop={'size': 8}, loc=0)
        axs[2].legend(prop={'size': 8}, loc=0)
        axs[3].legend(prop={'size': 8}, loc=0)

        axs[0].set_ylim(np.amin(ydata[0][50:-50]) * 0.95, np.amax(ydata[0][50:-50]) * 2)
        axs[1].set_ylim(-0.5, 0.5)
        axs[2].set_ylim(np.amin(ydata[0][50:-50]) * 0.05, np.amax(ydata[0][50:-50]) * 2)

        plt.tight_layout(h_pad=0)
        fig.savefig(plot_name, format='pdf', dpi=150)
        fig.clear()
        plt.close(fig)
        del fig
        log.debug(" fits1d pdf saved to {}".format(plot_name))

        if save_data:
            log.debug(" saving fits1d data to .npy file")
            out_data = np.array([ydata[0], xdata[0], ydata[3], ydata[2], ydata[5],
                                 ydata[1], ydata[4]])
            out_data_name = os.path.splitext(plot_name)[0] + '_data.npy'
            np.save(out_data_name, out_data)
            log.debug(" fits1d .npy file saved to {}".format(out_data_name))

    else:
        log.debug("diagnostic_plot: plot type {} not recognised, skipping".format(type))

    return None

