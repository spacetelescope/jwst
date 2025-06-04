Step Arguments
==============
The ``emicorr`` step has the following step-specific arguments.

``--algorithm`` (string, default='joint')
    Algorithm for fitting EMI noise in ramps.  If 'sequential', ramps are fit
    and then EMI noise is fit to the residuals.  If 'joint', ramps and noise
    are fit simultaneously.  The sequential algorithm can be used without
    a reference waveform, generating a new reference file on-the-fly, but it
    requires 10 or more groups for a reliable fit.  The joint algorithm
    requires a reference waveform but can successfully fit EMI noise in ramps
    with 3 or more groups.

``--nints_to_phase`` (integer, default=None)
    Number of integrations to phase, when `algorithm` is 'sequential'.

``--nbins`` (integer, default=None)
    Number of bins in one phased wave, when `algorithm` is 'sequential'.

``--scale_reference`` (boolean, default=True)
    If True and `algorithm` is 'sequential', the reference wavelength will be scaled
    to the data's phase amplitude.

``--onthefly_corr_freq``  (list, default=None)
    Frequency values to use to create a correction on-the-fly.  If provided,
    any input EMICORR reference model is ignored and the `algorithm` is set to
    'sequential'.

``--use_n_cycles`` (integer, default=3)
    Number of cycles to use to calculate the phase when `algorithm` is 'sequential'.
    To use all cycles, set to None.

``--fit_ints_separately`` (boolean, default=False)
    If True, fit and remove EMI noise for each integration separately, when
    `algorithm` is 'joint'.

``--save_intermediate_results``  (string, default=False)
    If True, and input frequencies are provided in `onthefly_corr_freq`,
    save a reference file with the fit phase amplitudes for the provided frequencies
    to disk. The file will be in ASDF output with the same format as an
    EMICORR reference file.
