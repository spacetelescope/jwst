import numpy as np


def make_lowpass_filter(f_half_power, w_cutoff, n, d=1.0):
    """
    Make a lowpass Fourier filter
    
    Parameters: f_half_power, number
                  Half power frequency
                w_cutoff, number
                  Width of cosine cutoff. The response transitions
                  from 1x to 0x over this range of frequencies
                n, int
                  Number of samples in the timeseries
                d, float
                  Sample spacing (inverse of the sampling rate). Defaults to 1.
    """
    
    # Make frequencies vector
    freq = np.fft.rfftfreq(n, d=d)
    
    # Build a cosine wave that is approriately shifted
    cos = (1 + np.cos(np.pi * (freq - f_half_power) / w_cutoff + np.pi / 2)) / 2
    
    # Construct low-pass filter with cosine rolloff
    filt = np.where(freq <= f_half_power - w_cutoff / 2, 1, cos)
    filt = np.where(freq <= f_half_power + w_cutoff / 2, filt, 0)
    
    # Done
    return(filt)
