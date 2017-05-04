#!/usr/bin/env python
import os

from ..stpipe import Pipeline
from .. import datamodels

# calwebb Image3 step imports
from ..resample import resample_step
from ..skymatch import skymatch_step
from ..outlier_detection import outlier_detection_step
from ..source_catalog import source_catalog_step
from ..tweakreg_catalog import tweakreg_catalog_step
from ..tweakreg import tweakreg_step

__version__ = "0.7.0"
# Define logging
import logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)

class Image3Pipeline(Pipeline):
    """

    Image3Pipeline: Applies level 3 processing to imaging-mode data from
                    any JWST instrument.

    Included steps are:
        tweakreg_catalog
        tweakreg
        skymatch
        outlier_detection
        resample
        source_catalog

    """

    # Define alias to steps
    step_defs = {'resample': resample_step.ResampleStep,
                 'skymatch': skymatch_step.SkyMatchStep,
                 'outlier_detection': outlier_detection_step.OutlierDetectionStep,
                 'tweakreg': tweakreg_step.TweakRegStep,
                 'source_catalog': source_catalog_step.SourceCatalogStep,
                 'tweakreg_catalog': tweakreg_catalog_step.TweakregCatalogStep
                 }

    def process(self, input):

        log.info('Starting calwebb_image3 ...')

        input_models = datamodels.open(input)

        is_container = (type(input_models) == type(datamodels.ModelContainer()))
        if is_container and len(input_models.group_names) > 1:

            generate_2c_names(input_models)

            # perform full outlier_detection of ASN data
            log.info("Generating source catalogs for alignment...")
            input_models = self.tweakreg_catalog(input_models)
            log.info("Aligning input images...")
            input_models = self.tweakreg(input_models)
            log.info("Matching sky values across all input images...")
            input_models = self.skymatch(input_models)
            log.info("Performing outlier detection on input images...")
            input_models = self.outlier_detection(input_models)
            
            # Now clean up intermediate products which no are no longer needed
            for i in input_models:
                try:
                    catalog_name = i.meta.tweakreg_catalog.filename
                    os.remove(catalog_name)
                except:
                    pass

            log.info("Resampling {} to create combined "
                "product: {}".format(input, input_models.meta.resample.output))

        # Setup output file name for subsequent use
        # TODO: fix single resample to do what outlier detection does
        # in updating meta.resample.*
        output_file = mk_prodname(self.output_dir,
            input_models.meta.resample.output, 'i2d')
        input_models.meta.resample.output = output_file

        # Resample step always returns ModelContainer,
        # yet we only need the DataModel result
        output = self.resample(input_models)

        # create final source catalog from resampled output
        out_catalog = self.source_catalog(output)

        # Save the final image product
        log.info('Saving final image product to %s', output_file)
        output.save(output_file)
        log.info('... ending calwebb_image3')

        return

def generate_2c_names(input_models):
    """ Update the names of the input files with Level 2C names
    """

    if hasattr(input_models.meta,'asn_id'):
        asn_id = input_models.meta.asn_table.asn_id
    else:
        asn_id = "a3001"
        
    for i in input_models:
        i.meta.filename = i.meta.filename.replace('.fits','-{}.fits'.format(asn_id))

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
