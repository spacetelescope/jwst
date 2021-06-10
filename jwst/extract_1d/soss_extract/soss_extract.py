import logging

from stdatamodels import DataModel

from ... import datamodels
from .soss_syscor import make_background_mask, soss_background
from .soss_solver import solve_transform, apply_transform
# from .soss_engine import ExtractionEngine  # TODO placeholder, pending review of the Engine.

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def extract_image(scidata, scierr, scimask, transform=None, tikfac=None):  # TODO how best to pass reference files and additional parameters (e.g. threshold)?

    # Perform background correction.
    bkg_mask = make_background_mask(scidata)
    scidata_bkg, col_bkg, npix_bkg = soss_background(scidata, scimask, bkg_mask=bkg_mask)

    # TODO add 1/f correction?

    if transform is None:

        # Use the Solver on the image.
        transform = solve_transform(scidata, scimask, xref, yref, subarray)  # TODO xref, yref, subarray need to be passed to extract_image().

    # Use the transformation on the reference files.
    trans_map = apply_transform(transform, ref_map, ovs, pad)  # TODO adjust and duplicate for all 4 maps. Add as method to the relevant datamodels?

    # Initialize the Engine.
    engine = TrpzOverlap(p_list, lam_list, t_list=t_list, c_list=c_list, n_os=n_os, thresh=thresh)  # TODO placeholder, pending review of the Engine.

    if tikfac is None:

        # Find the tikhonov factor.
        # Initial pass 14 orders of magnitude.
        factors = np.logspace(-25, -12, 14)
        tiktests = engine.get_tikho_tests(factors, data=image_bkg, sig=imerr)  # TODO placeholder, pending review of the Engine.
        tikfac = engine.best_tikho_factor(test=tiktests)

        # Refine across 4 orders of magnitude.
        tikfac = np.log10(tikfac)
        factors = np.logspace(tikfac - 2, tikfac + 2, 20)
        tiktests = engine.get_tikho_tests(factors, data=image_bkg, sig=imerr)  # TODO placeholder, pending review of the Engine.
        tikfac = engine.best_tikho_factor(test=tiktests)

    # Run the extract method of the Engine.
    f_k = engine.extract(data=image_bkg, sig=imerr, mask=mask, tikhonov=True, factor=tikfac)  # TODO placeholder, pending review of the Engine.

    # Re-construct orders 1 and 2.
    model_order_1 = engine.rebuild(f_k, i_orders=[1])  # TODO placeholder, pending review of the Engine.
    model_order_2 = engine.rebuild(f_k, i_orders=[2])

    # Run a box extraction on the order 1/2 subtracted image.
    # TODO still being written.

    # TODO Placeholders for return values. Dictionaries with values for each order.
    wavelenghts = dict()
    fluxes = dict()
    fluxerrs = dict()

    return wavelengths, fluxes, fluxerrs, transform, tikfac  # TODO update return products when box extraction finished.


def run_extract1d(input_model: DataModel,
                  spectrace_ref_name: str,
                  wavemap_ref_name: str,
                  specprofile_ref_name: str,
                  speckernel_ref_name: str):  # TODO What other parameters are needed: threshold, oversampling etc.

    # Read the reference files.
    spectrace_ref = datamodels.SpecTraceModel(spectrace_ref_name)
    wavemap_ref = datamodels.WaveMapModel(wavemap_ref_name)
    specprofile_ref = datamodels.SpecProfileModel(specprofile_ref_name)
    speckernel_ref = datamodels.SpecKernelModel(speckernel_ref_name)

    if isinstance(input_model, datamodels.ImageModel):

        # Received a single 2D image.
        scidata = input_model.data
        scierr = input_model.err
        scimask = input_model.dq == 0

        # Perform the extraction.
        wavelengths, fluxes, fluxerrs, transform, tikfac = extract_image(scidata, scierr, scimask)

        # Initialize the output model.
        output_model = datamodels.MultiSpecModel()  # TODO is this correct for ImageModel input?
        output_model.update(input_model)  # Copy meta data from input to output.

        # Copy spectral data for each order into the output model.
        # TODO how to include parameters like transform and tikfac in the output.
        for order in wavelengths.keys():
            wavelength = wavelengths[order]
            flux = fluxes[order]
            fluxerr = fluxerrs[order]
            surf_bright = np.zeros_like(flux)
            sb_error = np.zeros_like(flux)
            dq = np.zeros_like(flux, dtype=np.uint32)
            background = np.zeros_like(flux)  # TODO we do compute a background but not per order.
            berror = np.zeros_like(flux)
            npixels = np.zeros_like(flux)

            out_table = np.array(list(zip(wavelength,
                                          flux, fluxerr,
                                          surf_bright, sb_error,
                                          dq, background, berror, npixels)),
                                 dtype=datamodels.SpecModel().spec_table.dtype)
            spec = datamodels.SpecModel(spec_table=out_table)
            output_model.spec.append(spec)

    elif isinstance(input_model, datamodels.CubeModel):

        nimages = len(input_model.data)  # TODO Do this or use meta.exposure.

        # Build deepstack out of max N images OPTIONAL.
        # TODO making a deepstack could be used to get a more robust transform and tikfac

        # Set transform and tikfac, will be computed first iteration only.
        transform = None
        tikfac = None

        # Initialize the output model.
        output_model = datamodels.MultiSpecModel()  # TODO is this correct for CubeModel input?
        output_model.update(input_model)  # Copy meta data from input to output.

        # Loop over images.
        for i in range(nimages):

            # Unpack the i-th image.
            scidata = input_model.data[i]
            scierr = input_model.err[i]
            scimask = input_model.dq[i] == 0

            # Perform the extraction.
            wavelengths, fluxes, fluxerrs, transform, tikfac = extract_image(scidata, scierr, scimask, transform=transform, tikfac=tikfac)

            # Copy spectral data for each order into the output model.
            # TODO how to include parameters like transform and tikfac in the output.
            for order in wavelengths.keys():
                wavelength = wavelengths[order]
                flux = fluxes[order]
                fluxerr = fluxerrs[order]
                surf_bright = np.zeros_like(flux)
                sb_error = np.zeros_like(flux)
                dq = np.zeros_like(flux, dtype=np.uint32)
                background = np.zeros_like(flux)  # TODO we do compute a background but not per order.
                berror = np.zeros_like(flux)
                npixels = np.zeros_like(flux)

                out_table = np.array(list(zip(wavelength,
                                              flux, fluxerr,
                                              surf_bright, sb_error,
                                              dq, background, berror, npixels)),
                                     dtype=datamodels.SpecModel().spec_table.dtype)
                spec = datamodels.SpecModel(spec_table=out_table)
                output_model.spec.append(spec)

    return output_model
