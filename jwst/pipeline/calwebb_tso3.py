from __future__ import unicode_literals, absolute_import

import os

from astropy.table import vstack

from ..stpipe import Pipeline
from .. import datamodels

from ..outlier_detection import outlier_detection_tso_step
from ..outlier_detection import outlier_detection_scaled_step
from ..tso_photometry import tso_photometry_step
from ..extract_1d import extract_1d_step
from ..white_light import white_light_step

__version__ = "0.7.8"


class TSO3Pipeline(Pipeline):
    """
    TSO3Pipeline: Applies level 3 processing to TSO-mode data from
                    any JWST instrument.

    Included steps are:
        outlier_detection
        tso_photometry
        extract_1d
        white_light
    """

    spec = """
        suffix = string(default='i2d')
        scale_detection = boolean(default=False)
    """

    # Define alias to steps
    step_defs = {'outlier_detection': 
                        outlier_detection_tso_step.OutlierDetectionTSOStep,
                 'outlier_detection_scaled': 
                        outlier_detection_scaled_step.OutlierDetectionScaledStep,
                 'tso_photometry': tso_photometry_step.TSOPhotometryStep,
                 'extract_1d': extract_1d_step.Extract1dStep,
                 'white_light': white_light_step.WhiteLightStep
                 }

    def process(self, input):
        """
        Run the TSO3Pipeline

        Parameters
        ----------
        input: Level3 Association, or ModelContainer
            The exposures to process
        """

        self.log.info('Starting calwebb_tso3 ...')
        print("scale_detection set to : {}".format(self.scale_detection))
        
        input_asn = datamodels.open(input)
        product_name = input_asn.meta.asn_table.products[0].name

        # Setup output creation
        #self.output_basename = product['name']
        self.output_basename = product_name 
        
        # Input may consist of multiple exposures, so loop over each of them
        for cube in input_asn:
            # Convert CubeModel into ModelContainer of 2-D DataModels
            input_models = datamodels.ModelContainer()
            for i in range(cube.data.shape[0]):
                # convert each plane of data cube into it own array
                # for outlier detection...
                image = datamodels.ImageModel(data=cube.data[i],
                        err=cube.err[i], dq=cube.dq[i])
                image.meta = cube.meta
                input_models.append(image)
            
            if not self.scale_detection:
                self.log.info("Performing outlier detection on input images...")
                results = self.outlier_detection(input_models)
            else:
                self.log.info("Performing scaled outlier detection on input images...")
                results = self.outlier_detection_scaled(input_models)
                
            # Transfer updated DQ values to original input observation
            for i in range(cube.data.shape[0]):
                # preserve output filename 
                orig_filename = cube.meta.filename
                # Update DQ arrays with those from outlier_detection step
                cube.dq[i] = input_models[i].dq
                # reset output filename to original value
                cube.meta.filename = orig_filename

        self.log.info("Writing Level 2c images with updated DQ arrays...")
        suffix_2c = 'crfints'
        for cube in input_asn:
            self.output_basename = cube.meta.filename
            self.save_model(cube, suffix=suffix_2c)

        # Create final photometry results as a single output
        # regardless of how many members their may be...
        phot_result_list = []
        if 'image' in results[0].meta.exposure.type.lower():
            # Create name for extracted photometry (Level 3) products
            phot_tab_name = "{}_phot.ecsv".format(product_name)

            # For each cube, extract photometry...
            for cube in input_asn:
                # Extract Photometry from imaging data
                phot_result_list.append(self.tso_photometry(cube))
        else:
            # Working with spectroscopic TSO data...
            # define output for x1d (level 2c/3) products
            x1d_models = datamodels.ModelContainer()
            x1d_models.update(input_asn)
            #x1d_models.meta = input_asn.meta

            # Create name for extracted white-light (Level 3) products
            phot_tab_name = "{}_whtlt.ecsv".format(product_name)

            # For each exposure in the TSO...
            for cube in input_asn:
                # Process spectroscopic TSO data
                # extract 1D
                self.log.info("Extracting 1-D spectra...")
                self.extract_1d.suffix = 'x1dints'
                result = self.extract_1d(cube)
                x1d_models.append(result)
                #
                # perform white-light photometry on 1d extracted data
                self.log.info("Performing white-light photometry...")
                phot_result_list.append(self.white_light(result))
                
            # Define suffix for stack of extracted 1d (level 2b) products
            #self.save_model(x1d_models,suffix="x1dints") # did not work!!!
            x1d_prod_name = "{}_x1dints.fits".format(product_name)
            self.log.info("Writing Level 3 X1DINTS product {}...".format(x1d_prod_name))
            x1d_models.write(x1d_prod_name)
                        
        phot_results = vstack(phot_result_list)
        phot_results.write(phot_tab_name, format='ascii.ecsv')
            
        return 
