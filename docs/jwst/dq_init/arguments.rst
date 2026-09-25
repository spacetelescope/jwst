Step Arguments
==============
The ``dq_init`` step has the following step-specific argument.

``--user_supplied_dq`` (string, default=None)
    User-supplied DQ reference file pre-populated with desired DQ flags to
    be inserted into the data being processed. This file must be a simple
    single-extension FITS containing 2D image HDU with the same dimension
    as the science frame. The data type of the DQ array in this file
    must be integer. It is up to the user to ensure the supplied flags are
    valid as no such checks will be done before applying them in the
    pipeline (see :ref:`Data Quality Flags`).
