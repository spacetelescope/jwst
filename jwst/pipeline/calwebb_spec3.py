#!/usr/bin/env python
import os

from ..stpipe import Pipeline
from .. import datamodels
from ..exp_to_source import exp_to_source
from ..resample import blend

# step imports
from ..skymatch import skymatch_step
from ..outlier_detection import outlier_detection_step
from ..resample import resample_spec_step
from ..cube_build import cube_build_step
from ..extract_1d import extract_1d_step

__version__ = "0.7.1"

# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class Spec3Pipeline(Pipeline):
    """
    Spec3Pipeline: Processes JWST spectroscopic exposures from Level 2b to 3.

    Included steps are:
    MIRI MRS background matching (skymatch)
    outlier detection (outlier_detection)
    2-D spectroscopic resampling (resample_spec)
    3-D spectroscopic resampling (cube_build)
    1-D spectral extraction (extract_1d)
    """

    spec = """
    """

    # Define aliases to steps
    step_defs = {
        'skymatch': skymatch_step.SkyMatchStep,
        'outlier_detection': outlier_detection_step.OutlierDetectionStep,
        'resample_spec': resample_spec_step.ResampleSpecStep,
        'cube_build': cube_build_step.CubeBuildStep,
        'extract_1d': extract_1d_step.Extract1dStep
    }

    # Main processing
    def process(self, input):
        """Entrypoint for this pipeline

        Parameters
        ----------
        input: str, Level3 Association, or DataModel
            The exposure or association of exposures to process
        """
        log.info('Starting calwebb_spec3 ...')

        # Retrieve the inputs:
        # could either be done via LoadAsAssociation and then manually
        # load input members into models and ModelContainer, or just
        # do a direct open of all members in ASN file, e.g.
        input_models = datamodels.open(input)

        # For the first round of development we will assume that the input
        # is ALWAYS an ASN. There's no use case for anyone ever running a
        # single exposure through.

        # Once data are loaded, store a few things for future use;
        # some of this is here only for the purpose of creating fake
        # products until the individual tasks work and do it themselves
        exptype = input_models[0].meta.exposure.type
        pool_name = input_models.meta.asn_table.asn_pool
        asn_file = input
        prod_template = input_models.meta.asn_table.products[0].name
        prog = input_models.meta.asn_table.program
        acid = input_models.meta.asn_table.asn_id

        # Convert to source-based products, when necessary
        if exptype in ['NRS_MSASPEC', 'NRC_GRISM', 'NIS_WFSS']:
            filter = input_models[0].meta.instrument.filter
            grating = input_models[0].meta.instrument.grating

            source_file_template = \
                'jw{prog}-{acid}_s{source}_nirspec_{filter}_{grating}_cal.fits'

            source_files = exp_to_source(input_models)

            for source in source_files:
                input_sources = source_files[source]

                # Save to a source-based product
                source_file_name = source_file_template.format(
                    prog=prog, acid=acid, source=source.lower(),
                    filter=filter.lower(), grating=grating.lower())
                input_sources.save(source_file_name)

                # Convert source model to a ModelContainer for input to
                # subsequent steps
                # Extracts data from each of the Exposures in MultiExposureModel
                # and packs into a ModelContainer
                raise RuntimeError('No code to convert MultiSource to ModelContainer')

                # Call outlier detection
                source_models = self.outlier_detection(source_models)

                # Need to somehow now translate the updated DQ flags contained
                # in the returned source_models back into DQ arrays of original
                # MultiExposureModel and save as a pseudo level-2c product
                raise RuntimeError('No code to update DQ arrays in exposure-based data')

                # Save updated inputs to level-2c products

                # Now send updated source_models that came back from outlier_detection,
                # which have CR DQ flags in them, into the resample_spec step. Note
                # that we will only ever call resample_spec on this branch, because
                # this branch only deals with multi-source 2D spectra, not IFU data.
                s2d_result = self.resample_spec(source_models)

                # Save this result as a 's2d' source-based product
                raise RuntimeError('No code to save s2d result')

                # Finally call extract_1d on resampled s2d product
                x1d_result = self.extract_1d_step(s2d_result)

                # Save this result to a source-based 'x1d' product
                raise RuntimeError('No code to save x1d result')

                # Done with this source-based input product. Loop to the next ...


        else:

            # We're back to the non-source-based branch of operations

            # Call the skymatch step for MIRI MRS data
            if exptype in ['MIR_MRS']:
                input_models = self.skymatch(input_models)

            # Call outlier detection
            input_models = self.outlier_detection(input_models)

            # Save updated inputs to level-2c products
            generate_2c_names(input_models)
            input_models.save()

            # Call 3-D resampling for IFU data
            if exptype in ['MIR_MRS', 'NRS_IFU']:
                resampled = self.cube_build(input_models)
                suffix = 's3d'

                # Temporary hack to create a fake resampled model
                log.warning('Creating fake resampled product until step is available')
                resampled = datamodels.IFUCubeModel()
                naxis = input_models[0].data.shape[0]
                array = np.ones((naxis, naxis, naxis))
                resampled.data = array
                resampled.dq = array
                resampled.err = array
                resampled.weightmap = array
                #resampled.update(input_models[0])
                resampled.meta = input_models[0].meta

            # Call 2-D resampling for non-IFU data
            else:
                resampled = self.resample_spec(input_models)
                suffix = 's2d'

                # Temporary hack to create a fake resample model
                resampled = datamodels.DrizProductModel()
                resampled.data = input_models[0].data
                resampled.wht = input_models[0].err
                resampled.con = input_models[0].dq
                resampled.relsens = input_models[0].relsens
                #resampled.update(input_models[0])
                resampled.meta = input_models[0].meta

            # Save resampled product
            resampled_name = mk_prodname(self.output_dir, prod_template, suffix)
            log.debug('Blending metadata for {}'.format(resampled_name))
            cal_files = get_2b_names(input_models)
            blend.blendfitsdata(cal_files, resampled)
            resampled.meta.asn.pool_name = pool_name
            resampled.meta.asn.table_name = asn_file
            resampled.meta.cal_step.outlier_detection = 'COMPLETE'
            resampled.meta.cal_step.resample = 'COMPLETE'
            if isinstance(resampled, datamodels.IFUCubeModel):
                resampled.meta.model_type = 'IFUCubeModel'
            else:
                resampled.meta.model_type = 'DrizProductModel'
            log.info('Saving resampled product to %s', resampled_name)
            resampled.save(resampled_name)

            # Do 1-D spectral extraction
            output_1d = self.extract_1d(resampled)

            # Save extracted product
            output_1d.meta.asn.pool_name = pool_name
            output_1d.meta.asn.table_name = asn_file
            extracted_name = mk_prodname(self.output_dir, prod_template, 'x1d')
            log.info('Saving extracted product to %s', extracted_name)
            output_1d.save(extracted_name)

        # We're done
        log.info('Ending calwebb_spec3')
        return


def generate_2c_names(input_models):
    """Update the names of the input files with Level 2C names
    """

    if hasattr(input_models.meta, 'asn_id'):
        asn_id = input_models.meta.asn_table.asn_id
    else:
        asn_id = "a3001"

    for i in input_models:
        i.meta.filename = i.meta.filename.replace('.fits', '-{}.fits'.format(asn_id))


def get_2b_names(input_models):
    """Construct a list of the input level-2b file names
    """

    file_list = []
    for i in input_models:
        file_list.append(i.meta.filename)

    return file_list


def mk_prodname(output_dir, filename, suffix):
    """
    Build a file name based on an ASN product name template.

    The input ASN product name is used as a template. A user-specified
    output directory path is prepended to the root of the product name.
    The input product type suffix is appended to the root of the input
    product name, preserving any existing file name extension
    (e.g. ".fits").

    Args:
        output_dir (str): The output_dir requested by the user
        filename (str): The input file name, to be reworked
        suffix (str): The desired file type suffix for the new file name

    Returns:
        string: The new file name

    Examples:
        For output_dir='/my/path', filename='jw12345_nrca_cal.fits', and
        suffix='i2d', the returned file name will be
        '/my/path/jw12345_nrca_cal_i2d.fits'
    """

    # If the user specified an output_dir, replace any existing
    # path with output_dir
    if output_dir is not None:
        dirname, filename = os.path.split(filename)
        filename = os.path.join(output_dir, filename)

    # Now append the new suffix to the root name, preserving
    # any existing extension
    base, ext = os.path.splitext(filename)
    if len(ext) == 0:
        ext = ".fits"
    return base + '_' + suffix + ext
