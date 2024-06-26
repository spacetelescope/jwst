Step Arguments
==============
The ``emicorr`` step has the following step-specific arguments.

``--nints_to_phase`` (integer, default=None)
    Number of integrations to phase.

``--nbins`` (integer, default=None)
    Number of bins in one phased wave.

``--scale_reference`` (boolean, default=True)
    If True, the reference wavelength will be scaled to the
    data's phase amplitude.

``--user_supplied_reffile`` (boolean, default=None)
    This is to specify an ASDF-format user-created reference file.

``--save_intermediate_results``  (string, default=False)
    This is a boolean flag to specify whether to write a step output
    file with the EMI correction, and a reference file with all the
    on-the-fly frequencies' phase amplitudes. The file will have an
    ASDF output with the same format as an EMI reference file.

``--onthefly_corr_freq``  (list, default=None)
    This is to tell the code to do correction for the frequencies in
    the list with a reference file created on-the-fly instead of CRDS.

``--use_n_cycles`` (integer, default=3)
    Number of cycles to use to calculate the phase. To use all
    integrations, set to None.
