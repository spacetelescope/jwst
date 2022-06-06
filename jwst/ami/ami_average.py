#
#  Module for averaging LG results for a set of AMI exposures
#
import logging
from .. import datamodels

log = logging.getLogger(__name__)
#log.addHandler(logging.NullHandler())
log.setLevel(logging.DEBUG)


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


    # Create the output model as a copy of the first input model
    log.debug(' Create output as copy of %s', lg_products[0])
    output_model = datamodels.AmiLgModel(lg_products[0]).copy()

    # Find the input product with the smallest fit_image attribute
    min_size = 2048
    for input in lg_products:
        prod = datamodels.AmiLgModel(input)
        if prod.fit_image.shape[0] < min_size:
            min_size = prod.fit_image.shape[0]
        prod.close()
    log.debug(' minimum size of fit_image=%d', min_size)

    # Loop over inputs, adding their values to the output
    for prod_num, input in enumerate(lg_products):

        log.debug(' Accumulate data from %s', input)
        prod = datamodels.AmiLgModel(input)

        prod_size = prod.fit_image.shape[0]
        if prod_size > min_size:
            trim = int((prod_size - min_size) / 2)
            log.debug(' trim fit and resid images by %d pixels', trim)
            fit_image = prod.fit_image[trim:-trim, trim:-trim]
            resid_image = prod.resid_image[trim:-trim, trim:-trim]
            if prod_num == 0:
                output_model.fit_image = fit_image
                output_model.resid_image = resid_image
        else:
            fit_image = prod.fit_image
            resid_image = prod.resid_image
        if prod_num > 0:
            output_model.fit_image += fit_image
            output_model.resid_image += resid_image
            output_model.closure_amp_table['coeffs'] += prod.closure_amp_table['coeffs']
            output_model.closure_phase_table['coeffs'] += prod.closure_phase_table['coeffs']
            output_model.fringe_amp_table['coeffs'] += prod.fringe_amp_table['coeffs']
            output_model.fringe_phase_table['coeffs'] += prod.fringe_phase_table['coeffs']
            output_model.pupil_phase_table['coeffs'] += prod.pupil_phase_table['coeffs']
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
