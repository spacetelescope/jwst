import logging

import numpy as np

from stdatamodels import DataModel

from ... import datamodels
from .soss_syscor import make_background_mask, soss_background
from .soss_solver import solve_transform, apply_transform
from .soss_engine import ExtractionEngine
from .engine_utils import ThroughputSOSS, WebbKernel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def get_ref_file_args(ref_files, transform):
    """Prepare the reference files for the extraction engine.

    :param ref_files: A dictionary of the reference file DataModels. # TODO not final?
    :param transform: A 3-elemnt list or array describing the rotation and
        translation to apply to the reference files in order to match the
        observation.

    :type ref_files: dict
    :type transform: array_like

    :returns: The reference file args used with the extraction engine.
    :rtype: Tuple(wavemaps, specprofiles, throughputs, kernels)
    """

    # The wavelength maps for order 1 and 2.
    wavemap_ref = ref_files['wavemap']

    ovs = wavemap_ref.map[0].oversampling
    pad = wavemap_ref.map[0].padding

    wavemap_o1 = apply_transform(transform, wavemap_ref.map[0].data, ovs, pad)
    wavemap_o2 = apply_transform(transform, wavemap_ref.map[1].data, ovs, pad)

    # The spectral profiles for order 1 and 2.
    specprofile_ref = ref_files['specprofile']
    ovs = specprofile_ref.profile[0].oversampling
    pad = specprofile_ref.profile[0].padding

    specprofile_o1 = apply_transform(transform, specprofile_ref.profile[0].data, ovs, pad, norm=True)
    specprofile_o2 = apply_transform(transform, specprofile_ref.profile[1].data, ovs, pad, norm=True)

    # The throughput curves for order 1 and 2.
    spectrace_ref = ref_files['spectrace']

    throughput_o1 = ThroughputSOSS(spectrace_ref.trace[0].data['WAVELENGTH'], spectrace_ref.trace[0].data['THROUGHPUT'])
    throughput_o2 = ThroughputSOSS(spectrace_ref.trace[1].data['WAVELENGTH'], spectrace_ref.trace[1].data['THROUGHPUT'])

    # The spectral kernels.
    speckernel_ref = ref_files['speckernel']
    ovs = speckernel_ref.meta.spectral_oversampling
    n_pix = 2*speckernel_ref.meta.halfwidth + 1

    kernels_o1 = WebbKernel(speckernel_ref.wavelengths, speckernel_ref.kernels, wavemap_o1, ovs, n_pix)
    kernels_o2 = WebbKernel(speckernel_ref.wavelengths, speckernel_ref.kernels, wavemap_o2, ovs, n_pix)

    return [wavemap_o1, wavemap_o2], [specprofile_o1, specprofile_o2], [throughput_o1, throughput_o2], [kernels_o1, kernels_o2]


# TODO how best to pass reference files and additional parameters (e.g. threshold)?
def extract_image(scidata, scierr, scimask, ref_files, transform=None,
                  tikfac=None, n_os=5, threshold=1e-4, width=40):
    """Perform the spectral extraction on a single image.

    :param scidata: A single NIRISS SOSS detector image.
    :param scierr: The uncertainties corresponding to the detector image.
    :param scimask: Pixel that should be masked from the detectot image.
    :param ref_files: A dictionary of the reference file DataModels. # TODO not final?
    :param transform: A 3-elemnt list or array describing the rotation and
        translation to apply to the reference files in order to match the
        observation. If None the transformation is computed.
    :param tikfac: The Tikhonov regularization factor used when solving for
        the uncontaminated flux.
    :param n_os: The oversampling factor of the wavelength grid used when
        solving for the uncontaminated flux.
    :param threshold: The threshold value for using pixels based on the spectral
        profile.
    :param width: The width of the aperture used to extract the un-contaminated spectrum.

    :type scidata: array[float]
    :type scierr: array[float]
    :type scimask: array[float]
    :type ref_files: dict
    :type transform: array_like
    :type tikfac: float
    :type n_os: int
    :type threshold: float
    :type width: float

    :returns: TODO TBD
    :rtype: TODO TBD
    """

    # Perform background correction.
    bkg_mask = make_background_mask(scidata, width=40)
    scidata_bkg, col_bkg, npix_bkg = soss_background(scidata, scimask, bkg_mask=bkg_mask)

    # Some error values are 0, we need to mask those pixels for the ectraction engine.
    scimask = scimask | ~(scierr > 0)

    # TODO add 1/f correction?

    # TODO placing the tranform (and the call to get_ref_file_args()) in run_extract_1d might be better.
    if transform is None:

        log.info('Solving for the transformation parameters.')

        # Unpack the expected order 1 positions.
        spectrace_ref = ref_files['spectrace']
        xref = spectrace_ref.trace[0].data['X']
        yref = spectrace_ref.trace[0].data['Y']
        subarray = spectrace_ref.meta.subarray.name  # TODO better way of propagating the subarray?

        # Use the Solver on the image.
        transform = solve_transform(scidata, scimask, xref, yref, subarray)

    log.info('Using transformation parameters {}'.format(transform))

    # Prepare the reference file arguments.
    ref_file_args = get_ref_file_args(ref_files, transform)

    # Initialize the Engine.
    # TODO set c_kwargs?
    engine = ExtractionEngine(*ref_file_args, n_os=n_os, threshold=threshold)

    if tikfac is None:

        log.info('Solving for the optimal Tikhonov factor.')

        # Find the tikhonov factor.
        # Initial pass 14 orders of magnitude.
        factors = np.logspace(-25, -12, 14)
        tiktests = engine.get_tikho_tests(factors, data=scidata_bkg, error=scierr, mask=scimask)
        tikfac = engine.best_tikho_factor(tests=tiktests)

        # Refine across 4 orders of magnitude.
        tikfac = np.log10(tikfac)
        factors = np.logspace(tikfac - 2, tikfac + 2, 20)
        tiktests = engine.get_tikho_tests(factors, data=scidata_bkg, error=scierr, mask=scimask)
        tikfac = engine.best_tikho_factor(tests=tiktests)

    log.info('Using a Tikhonov factor of {}'.format(tikfac))

    # Run the extract method of the Engine.
    f_k = engine.extract(data=scidata_bkg, error=scierr, mask=scimask, tikhonov=True, factor=tikfac)

    # Compute the log-likelihood of the best fit.
    logl = engine.compute_likelihood(f_k, same=False)

    log.info('Optimal solution has a log-likelihood of {}'.format(logl))

    # Re-construct orders 1 and 2.
    # TODO get this mask from the DQ array? Ensure it works for all subarrays.
    border_mask = np.zeros_like(scidata)
    border_mask[-4:] = True
    border_mask[:, :4] = True
    border_mask[:, -4:] = True

    model = ExtractionEngine(*ref_file_args, wave_grid=engine.wave_grid, threshold=1e-5, global_mask=border_mask)

    model_order_1 = model.rebuild(f_k, i_orders=[0])
    model_order_2 = model.rebuild(f_k, i_orders=[1])

    # Run a box extraction on the order 1/2 subtracted image for orders 1,2 and 3.
    # TODO still being written, use width here.

    # TODO Placeholders for return values.
    wavelengths = dict()
    fluxes = dict()
    fluxerrs = dict()

    # TODO update return products when box extraction finished.
    return wavelengths, fluxes, fluxerrs, transform, tikfac


def run_extract1d(input_model: DataModel,
                  spectrace_ref_name: str,
                  wavemap_ref_name: str,
                  specprofile_ref_name: str,
                  speckernel_ref_name: str,
                  soss_kwargs: dict):
    """Run the spectral extraction on NIRISS SOSS data.

    :param input_model:
    :param spectrace_ref_name:
    :param wavemap_ref_name:
    :param specprofile_ref_name:
    :param speckernel_ref_name:

    :type input_model:
    :type spectrace_ref_name:
    :type wavemap_ref_name:
    :type specprofile_ref_name:
    :type speckernel_ref_name:

    :returns: An output_model containing the extracted spectra.
    :rtype:
    """

    # Read the reference files.
    spectrace_ref = datamodels.SpecTraceModel(spectrace_ref_name)
    wavemap_ref = datamodels.WaveMapModel(wavemap_ref_name)
    specprofile_ref = datamodels.SpecProfileModel(specprofile_ref_name)
    speckernel_ref = datamodels.SpecKernelModel(speckernel_ref_name)

    ref_files = dict()
    ref_files['spectrace'] = spectrace_ref
    ref_files['wavemap'] = wavemap_ref
    ref_files['specprofile'] = specprofile_ref
    ref_files['speckernel'] = speckernel_ref

    if isinstance(input_model, datamodels.ImageModel):

        log.info('Input is an ImageModel, processing a single integration.')

        # Received a single 2D image.
        scidata = input_model.data.astype('float64')  # TODO eewww.
        scierr = input_model.err.astype('float64')
        scimask = input_model.dq > 0  # Mask bad pixels with True.

        # Perform the extraction.
        result = extract_image(scidata, scierr, scimask, ref_files, transform=None, tikfac=soss_kwargs['tikfac'], n_os=soss_kwargs['n_os'], threshold=soss_kwargs['threshold'])
        wavelengths, fluxes, fluxerrs, transform, soss_kwargs['tikfac'] = result

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

        log.info('Input is a CubeModel containing {} integrations.'.format(nimages))

        # Build deepstack out of max N images OPTIONAL.
        # TODO making a deepstack could be used to get a more robust transform and tikfac, 1/f.

        # Set transform and tikfac, will be computed first iteration only.
        transform = None

        # Initialize the output model.
        output_model = datamodels.MultiSpecModel()  # TODO is this correct for CubeModel input?
        output_model.update(input_model)  # Copy meta data from input to output.

        # Loop over images.
        for i in range(nimages):

            log.info('Processing integration {} of {}.'.format(i + 1, nimages))

            # Unpack the i-th image.
            scidata = input_model.data[i].astype('float64')  # TODO eewww.
            scierr = input_model.err[i].astype('float64')
            scimask = input_model.dq[i] > 0

            # Perform the extraction.
            result = extract_image(scidata, scierr, scimask, ref_files, transform=transform, tikfac=soss_kwargs['tikfac'], n_os=soss_kwargs['n_os'], threshold=soss_kwargs['threshold'])
            wavelengths, fluxes, fluxerrs, transform, soss_kwargs['tikfac'] = result

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

    else:
        msg = "Only ImageModel and CubeModel are implemented for the NIRISS SOSS extraction."
        raise ValueError(msg)

    return output_model
