import os.path as op

from astropy.table import vstack

from ..stpipe import Pipeline
from .. import datamodels

from ..outlier_detection import outlier_detection_step
from ..tso_photometry import tso_photometry_step
from ..extract_1d import extract_1d_step
from ..white_light import white_light_step

__all__ = ['Tso3Pipeline']


class Tso3Pipeline(Pipeline):
    """
    TSO3Pipeline: Applies level 3 processing to TSO-mode data from
                    any JWST instrument.

    Included steps are:

        * outlier_detection
        * tso_photometry
        * extract_1d
        * white_light
    """

    spec = """
        scale_detection = boolean(default=False)
    """

    # Define alias to steps
    step_defs = {'outlier_detection':
                 outlier_detection_step.OutlierDetectionStep,
                 'tso_photometry': tso_photometry_step.TSOPhotometryStep,
                 'extract_1d': extract_1d_step.Extract1dStep,
                 'white_light': white_light_step.WhiteLightStep
                 }
    image_exptypes = ['NRC_TSIMAGE']
    reference_file_types = ['gain', 'readnoise']

    def process(self, input):
        """
        Run the TSO3Pipeline

        Parameters
        ----------
        input: Level3 Association, json format
            The exposures to process
        """

        self.log.info('Starting calwebb_tso3...')
        input_models = datamodels.open(input)

        if self.output_file is None:
            self.output_file = input_models.meta.asn_table.products[0].name
        self.asn_id = input_models.meta.asn_table.asn_id

        input_exptype = None
        # Input may consist of multiple exposures, so loop over each of them
        for cube in input_models:
            if input_exptype is None:
                input_exptype = cube.meta.exposure.type
            # Convert CubeModel into ModelContainer of 2-D DataModels
            input_2dmodels = datamodels.ModelContainer()
            for i in range(cube.data.shape[0]):
                # convert each plane of data cube into it own array
                # for outlier detection...
                image = datamodels.ImageModel(data=cube.data[i],
                                              err=cube.err[i], dq=cube.dq[i])
                image.update(cube)
                image.meta.wcs = cube.meta.wcs
                input_2dmodels.append(image)

            if not self.scale_detection:
                msg = "Performing outlier detection on input images..."
                self.log.info(msg)
                input_2dmodels = self.outlier_detection(input_2dmodels)

                # Transfer updated DQ values to original input observation
                for i in range(cube.data.shape[0]):
                    # Update DQ arrays with those from outlier_detection step
                    cube.dq[i] = input_2dmodels[i].dq
                cube.meta.cal_step.outlier_detection = \
                    input_2dmodels[0].meta.cal_step.outlier_detection

            else:
                msg = "Performing scaled outlier detection on input images..."
                self.log.info(msg)
                self.outlier_detection.scale_detection = True
                cube = self.outlier_detection(cube)

        if input_models[0].meta.cal_step.outlier_detection == 'COMPLETE':
            self.log.info("Writing Level 2c cubes with updated DQ arrays...")
            for cube in input_models:
                # preserve output filename
                original_filename = cube.meta.filename
                self.save_model(
                    cube, output_file=original_filename, suffix='crfints',
                    asn_id=input_models.meta.asn_table.asn_id
                )
                cube.meta.filename = original_filename

        # Create final photometry results as a single output
        # regardless of how many members there may be...
        phot_result_list = []
        if input_exptype in self.image_exptypes:
            # Create name for extracted photometry (Level 3) product
            phot_tab_suffix = 'phot'

            for cube in input_models:
                # Extract Photometry from imaging data
                phot_result_list.append(self.tso_photometry(cube))
        else:
            # Create name for extracted white-light (Level 3) product
            phot_tab_suffix = 'whtlt'

            # Working with spectroscopic TSO data...
            # define output for x1d (level 3) products
            x1d_result = datamodels.MultiSpecModel()
            # TODO: check to make sure the following line is working
            x1d_result.update(input_models[0])

            # For each exposure in the TSO...
            for cube in input_models:
                # Process spectroscopic TSO data
                # extract 1D
                self.log.info("Extracting 1-D spectra...")
                result = self.extract_1d(cube)
                x1d_result.spec.extend(result.spec)

                # perform white-light photometry on 1d extracted data
                self.log.info("Performing white-light photometry...")
                phot_result_list.append(self.white_light(result))

            # Update some metadata from the association
            x1d_result.meta.asn.pool_name = \
                input_models.meta.asn_table.asn_pool
            x1d_result.meta.asn.table_name = op.basename(input)

            # Save the final x1d Multispec model
            self.save_model(x1d_result, suffix='x1dints')

        phot_results = vstack(phot_result_list)
        phot_tab_name = self.make_output_path(suffix=phot_tab_suffix, ext='ecsv')
        self.log.info("Writing Level 3 photometry catalog {}...".format(
                      phot_tab_name))
        phot_results.write(phot_tab_name, format='ascii.ecsv')

        return
