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
    file with the EMI correction. Additionaly, if the flag
    ``user_supplied_reffile`` is None and no CRDS reference file was
    found, all the on-the-fly frequencies phase amplitudes will be
    saved to an ASDF output with the same format as an EMI reference
    file.

