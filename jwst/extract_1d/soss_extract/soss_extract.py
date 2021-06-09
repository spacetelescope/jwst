import logging

from stdatamodels import DataModel

from .soss_syscor import make_background_mask, soss_background
from .soss_solver import solve_transform, apply_transform
from .soss_engine import ExtractionEngine

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def extract_image(scidata, scierr, scimask, transform=None, tikfac=None):  # TODO how best to padd reference files and additional parameters (e.g. threshold)?

    # Perform background correction.
    bkg_mask = make_background_mask(scidata)
    scidata_bkg, col_bkg, npix_bkg = soss_background(scidata, scimask, bkg_mask=bkg_mask)

    # TODO add 1/f correction?

    if transform is None:

        # Use the Solver on the image.
        transform = solve_transform(scidata, scimask, xref, yref, subarray)  # TODO xref, yref, subarray do not exist.

    # Use the transformation on the reference files.
    trans_map = apply_transform(transform, ref_map, ovs, pad)  # TODO adjust and duplicate for all 4 maps.

    # Initialize the Engine.
    engine = TrpzOverlap(p_list, lam_list, t_list=t_list, c_list=c_list, n_os=n_os, thresh=thresh)  # TODO update call signature

    if tikfac is None:

        # Find the tikhonov factor.
        # Initial pass 14 orders of magnitude.
        factors = np.logspace(-25, -12, 14)
        tiktests = engine.get_tikho_tests(factors, data=image_bkg, sig=imerr)  # TODO use mask here?
        tikfac = engine.best_tikho_factor(test=tiktests)

        # Refine across 4 orders of magnitude.
        tikfac = np.log10(tikfac)
        factors = np.logspace(tikfac - 2, tikfac + 2, 20)
        tiktests = engine.get_tikho_tests(factors, data=image_bkg, sig=imerr)  # TODO use mask here?
        tikfac = engine.best_tikho_factor(test=tiktests)

    # Run the extract method of the Engine.
    f_k = engine.extract(data=image_bkg, sig=imerr, mask=mask, tikhonov=True, factor=tikfac)  # TODO update call signature.

    # Re-construct orders 1 and 2.
    model_order_1 = engine.rebuild(f_k, i_orders=[1])  # TODO update call signature
    model_order_2 = engine.rebuild(f_k, i_orders=[2])

    # Run a box extraction on the order 1/2 subtracted image.
    # TODO still being written.

    return wavelength, flux, fluxerr, transform, tikfac  # TODO update return proucts when box extraction finished.


def run_extract1d(input_model: DataModel,
                  spectrace_ref_name: str,
                  wavemap_ref_name: str,
                  specprofile_ref_name: str,
                  speckernel_ref_name: str):  # TODO What other parameters are needed: threshold, oversampling etc.

    log.info('Our glorious NIRISS SOSS extraction algorithm will go here.')

    # Read the reference files.

    # Check if input_model is a single image or a cube.

    # Build deepstack out of max N images and use this to find the transform and tikhonov factor.

    # If cube loop over each image.

    # Place the results into the output_model.
    # Include simple transform, tikfac.

    # TODO placeholder.
    output_model = []

    return output_model
