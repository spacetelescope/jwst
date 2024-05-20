from __future__ import annotations
from typing import Protocol
import numpy as np
import jwst.datamodels as dm
from jwst.outlier_detection.outlier_detection_ifu import medfilt


def flagger_injector(mode: str) -> BackgroundFlagger:
    """
    Inject the appropriate flagger based on the mode

    Parameters
    ----------
    mode: str
        Mode of the data

    Returns
    -------
    BackgroundFlagger
        Flagging class
    """
    if mode.lower() == 'mrs':
        return MRSMedfiltBackgroundFlagger
    elif mode.lower() == 'mrs polyfit':
        return MRSPolyfitBackgroundFlagger
    else:
        raise ValueError(f"Mode {mode} not recognized")


def badpix_selfcal(input_model: dm.DataModel, flagfrac: float = 0.001) -> dm.DataModel:

    if flagfrac < 0 or flagfrac >=0.5:
        raise ValueError("flagfrac must be between 0 and 0.5. \
                        Note this fraction will be flagged on both high and low ends.")

    shp = input_model.data.shape
    basex,basey = np.meshgrid(np.arange(shp[0]),np.arange(shp[1]))
    _,_,lam=input_model.meta.wcs.transform('detector','world',basex,basey)

    bg = input_model.data
    # to do: update this so it operates on an association
    # where there are multiple exposures; make a 3-D background array
    # and then take the min over the exposures
    # medbg = np.nanmin(bg, axis=0)

    flagged_indices = flagger_injector('MRS').flag_background(bg, lam, flagfrac)
    input_model.data[flagged_indices] = np.nan
    input_model.dq[flagged_indices] = 1

    return input_model


class Channel:

    def __init__(self, 
                 bg: np.ndarray, 
                 wls: np.ndarray, 
                 wl_min: float = 0.0, 
                 wl_max: float = np.inf):
        """
        Cut full input data and wavelengths to the desired wl range

        Parameters
        ----------
        bg: np.ndarray
            Background data

        wls: np.ndarray
            Wavelengths of the data

        wl_min: float
            Minimum wavelength to keep

        wl_max: float
            Maximum wavelength to keep
        """
        self.wl_min = wl_min
        self.wl_max = wl_max
        self.good = (wls >= wl_min) & (wls < wl_max)
        self.wls = wls[self.good]
        self.bg = bg[self.good]


    def remove_polyfit(self, order: int = 2):
        """
        Remove low-order background using polynomial fit

        Parameters
        ----------
        order: int
            Order of the polynomial fit

        Returns
        -------
        np.ndarray
            Background data with the polynomial fit removed
        """
        xtemp=self.wls.ravel()
        ytemp=self.bg.ravel()
        indx=np.where((np.isfinite(xtemp))&(np.isfinite(ytemp)))

        temp=np.polyfit(xtemp[indx], ytemp[indx], order)
        temp2=np.poly1d(temp)
        yfit=temp2(self.wls)
        
        return self.bg - yfit



class BackgroundFlagger(Protocol):

    def flag_background(self, 
                        bg: np.ndarray, 
                        lam: np.ndarray, 
                        flagfrac: float) -> np.ndarray:
        ...


class MRSPolyfitBackgroundFlagger(BackgroundFlagger):

    def __init__(self):

        self.pedestal_indices = (138,911,495,514)


    def flag_background(self, medbg: np.ndarray, lam: np.ndarray, flagfrac: float) -> np.ndarray:
        """
        Flag outlier pixels in the MRS data background

        Parameters
        ----------
        medbg: np.ndarray
            Background data of shape (x, y), i.e., after median has
            already been taken over the number of exposures
        lam: np.ndarray
            Wavelengths of the data, same shape as medbg
        flagfrac: float
            Fraction of outlier pixels to flag on both high and low end
        """

        # subtract inter-slice pedestal dark
        ip = self.pedestal_indices
        pedestal_dark_value = np.nanmedian(medbg[ip[0]:ip[1], ip[2]:ip[3]])
        medbg -= pedestal_dark_value

        # Remove average background separately in low and high channels
        lam_cutoff = np.nanmean(lam)
        channels = [Channel(medbg, lam, wl_max=lam_cutoff), 
                    Channel(medbg, lam, wl_min=lam_cutoff),]
        for channel in channels:
            bg = channel.remove_polyfit()
            medbg[channel.good] = bg

        # Flag outliers using percentile cutoff
        flag_low, flag_high = np.nanpercentile(medbg, [flagfrac*100, (1-flagfrac)*100])
        bad = (medbg > flag_high) | (medbg < flag_low)
        indx=np.where(bad)
        return indx
    

class MRSMedfiltBackgroundFlagger(BackgroundFlagger):

    def __init__(self, kernel_size: int = 15):

        self.kernel_size = kernel_size
        self.pedestal_indices = (138,911,495,514)


    def flag_background(self, medbg: np.ndarray, lam: np.ndarray, flagfrac: float) -> np.ndarray:
        """
        Flag outlier pixels in the MRS data background

        Parameters
        ----------
        medbg: np.ndarray
            Background data of shape (x, y), i.e., after median has
            already been taken over the number of exposures
        lam: np.ndarray
            Wavelengths of the data, same shape as medbg
        flagfrac: float
            Percentile above which to flag
        """

        # subtract inter-slice pedestal dark
        ip = self.pedestal_indices
        pedestal_dark_value = np.nanmedian(medbg[ip[0]:ip[1], ip[2]:ip[3]])
        medbg -= pedestal_dark_value

        # Remove average background separately in low and high channels
        kernel=np.array([1,self.kernel_size])
        smoothed = medfilt(medbg, kernel)
        medbg = medbg - smoothed

        # Flag outliers using percentile cutoff
        # TO DO: add negative outliers
        # Flag outliers using percentile cutoff
        flag_low, flag_high = np.nanpercentile(medbg, [flagfrac*100, (1-flagfrac)*100])
        bad = (medbg > flag_high) | (medbg < flag_low)
        indx=np.where(bad)
        return indx