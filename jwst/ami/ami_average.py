#
#  Module for averaging LG results for a set of AMI exposures
#
import logging

from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)
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
    log.debug('create output as copy of %s', lg_products[0])
    output_model = datamodels.AmiLgModel(lg_products[0]).copy()

    # Find the input product with the smallest fit_image image size
    sizes = []
    for input in lg_products:
        prod = datamodels.AmiLgModel(input)
        sizes.append(prod.fit_image.shape[0])
        prod.close()
    min_size = min(sizes)
    log.debug('minimum size of fit_image=%d', min_size)

    # Loop over inputs, adding their values to the output
    for prod_num, input in enumerate(lg_products):

        log.info('Accumulate data from %s', input)
        prod = datamodels.AmiLgModel(input)

        prod_size = prod.fit_image.shape[0]
        if prod_size > min_size:
            # If the fit_image in this input is bigger than the minimum,
            # trim extra rows/cols from the edges
            trim = int((prod_size - min_size) / 2)
            log.debug('trim fit and resid images by %d pixels', trim)
            fit_image = prod.fit_image[trim:-trim, trim:-trim]
            resid_image = prod.resid_image[trim:-trim, trim:-trim]

            # If this is the first input, which was copied to create
            # the output, replace the original images with the trimmed
            # versions.
            if prod_num == 0:
                output_model.fit_image = fit_image
                output_model.resid_image = resid_image
        else:
            # Otherwise use the input images as is
            fit_image = prod.fit_image
            resid_image = prod.resid_image

        # For inputs past the first, accumulate the data into the output
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
    log.debug('Divide accumulated results by %d', len(lg_products))
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
