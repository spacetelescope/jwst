import logging
import numpy as np

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def apply_linearity_func(ramparr, dqarr, coeffarr, dq_flag):
    """
    Short Summary
    -------------
    Apply linearity correction to individual groups in ramparr where dq
    has not been flagged as saturated in the saturation step with dq_flag.

    Reference:
    http://www.stsci.edu/hst/wfc3/documents/handbooks/currentDHB/wfc3_Ch606.html
    WFC3 Data Handbook, Chapter 6, WFC3-IR Error Sources, 6.5 Detector
    Nonlinearity issues, 2.1 May 2011

    Scorr = Smeasured(1+A+B*Smeasured+C*Smeasured^2+D*Smeasured^3+...)

    Parameters
    ----------
    ramparr: 4D array containing ramp data

    dqarr: 4D array containing DQ information

    coeffarr: 3D array containing pixel-by-pixel linearity coefficient values
              for each term in the polynomial fit.

    dq_flag: Integer value representing the saturation flag value that was
             applied to dqarr by the saturation step

    Returns
    -------
    ramparr: Linearity corrected 4D array containing ramp data

    """

    # Retrieve the ramp data cube characteristics
    nints, ngroups, nrows, ncols = ramparr.shape

    # Number of coeffs is equal to the number of planes in coeff cube
    ncoeffs = coeffarr.shape[0]

    # Apply the linearity correction one integration at a time.
    for ints in range(nints):

        # Apply the linearity correction one group at a time
        for plane in range(ngroups):

            # Accumulate the polynomial terms into the corrected counts
            scorr = coeffarr[ncoeffs - 1] * ramparr[ints, plane]
            for j in range(ncoeffs - 2, 0, -1):
                scorr = (scorr + coeffarr[j]) * ramparr[ints, plane]
            scorr = coeffarr[0] + scorr

           # Only use the corrected signal where the original signal value
           # has not been flagged by the saturation step.
           # Otherwise use the original signal.
            ramparr[ints, plane, :, :] = \
                np.where(np.bitwise_and(dqarr[ints, plane, :, :], dq_flag),\
                        ramparr[ints, plane, :, :], scorr)

    del scorr

    return ramparr
