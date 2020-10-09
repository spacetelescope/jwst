#
#  Module for averaging LG results for a set of AMI exposures
#
import logging
from .. import datamodels

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def average_LG(lg_products):
    """
    Short Summary
    -------------
    Averages the LG results for a set of AMI exposures

    Parameters
    ----------
    lg_products: file list
        List of LG product file names to be averaged

    Returns
    -------
    output_model: Fringe model object
        Averaged fringe data
    """

    # Create the ouput model as a copy of the first input model
    log.debug(' Create output as copy of %s', lg_products[0])
    output_model = datamodels.AmiLgModel(lg_products[0]).copy()

    # Loop over remaining list of inputs, adding their values to the output
    for file in lg_products[1:]:

        log.debug(' Accumulate data from %s', file)
        prod = datamodels.AmiLgModel(file)

        output_model.fit_image += prod.fit_image
        output_model.resid_image += prod.resid_image
        output_model.closure_amp_table['coeffs'] += \
                     prod.closure_amp_table['coeffs']
        output_model.closure_phase_table['coeffs'] += \
                     prod.closure_phase_table['coeffs']
        output_model.fringe_amp_table['coeffs'] += \
                     prod.fringe_amp_table['coeffs']
        output_model.fringe_phase_table['coeffs'] += \
                     prod.fringe_phase_table['coeffs']
        output_model.pupil_phase_table['coeffs'] += \
                     prod.pupil_phase_table['coeffs']
        output_model.solns_table['coeffs'] += prod.solns_table['coeffs']

        prod.close()

    # Take the average of the accumulated results
    log.debug(' Divide accumulated results by %d', len(lg_products))
    output_model.fit_image /= len(lg_products)
    output_model.resid_image /= len(lg_products)
    output_model.closure_amp_table['coeffs'] /= len(lg_products)
    output_model.closure_phase_table['coeffs'] /= len(lg_products)
    output_model.fringe_amp_table['coeffs'] /= len(lg_products)
    output_model.fringe_phase_table['coeffs'] /= len(lg_products)
    output_model.pupil_phase_table['coeffs'] /= len(lg_products)
    output_model.solns_table['coeffs'] /= len(lg_products)

    # Return the averaged model
    return output_model
