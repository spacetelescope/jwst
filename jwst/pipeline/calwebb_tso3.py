from pathlib import Path

import numpy as np
from astropy.table import vstack

from stdatamodels.jwst import datamodels


from jwst.stpipe import Pipeline

from jwst.outlier_detection import outlier_detection_step
from jwst.tso_photometry import tso_photometry_step
from jwst.extract_1d import extract_1d_step
from jwst.white_light import white_light_step
from jwst.photom import photom_step
from jwst.pixel_replace import pixel_replace_step
from jwst.lib.pipe_utils import is_tso
from astropy.io.fits import FITS_rec

__all__ = ["Tso3Pipeline"]


class Tso3Pipeline(Pipeline):
    """
    Apply level 3 processing to TSO-mode data from any JWST instrument.

    Included steps are:

        * outlier_detection
        * tso_photometry
        * pixel_replace
        * extract_1d
        * photom
        * white_light
    """

    class_alias = "calwebb_tso3"

    spec = """
    """  # noqa: E501

    # Define alias to steps
    step_defs = {
        "outlier_detection": outlier_detection_step.OutlierDetectionStep,
        "tso_photometry": tso_photometry_step.TSOPhotometryStep,
        "pixel_replace": pixel_replace_step.PixelReplaceStep,
        "extract_1d": extract_1d_step.Extract1dStep,
        "photom": photom_step.PhotomStep,
        "white_light": white_light_step.WhiteLightStep,
    }
    reference_file_types = ["gain", "readnoise"]

    def process(self, input_data):
        """
        Run the TSO3Pipeline on the input association.

        Parameters
        ----------
        input_data : Level3 Association, json format
            The exposures to process.
        """
        self.log.info("Starting calwebb_tso3...")
        asn_exptypes = ["science"]

        input_models = datamodels.open(input_data, asn_exptypes=asn_exptypes)

        # Sanity check the input data
        input_tsovisit = is_tso(input_models[0])
        if not input_tsovisit:
            self.log.error("INPUT DATA ARE NOT TSO MODE. ABORTING PROCESSING.")
            return

        if self.output_file is None:
            self.output_file = input_models.asn_table["products"][0]["name"]

        # This asn_id assignment is important as it allows outlier detection
        # to know the asn_id since that step receives the cube as input.
        self.asn_id = input_models.asn_table["asn_id"]
        self.outlier_detection.mode = "tso"

        # Input may consist of multiple exposures, so loop over each of them
        input_exptype = None
        for i, cube in enumerate(input_models):
            if input_exptype is None:
                input_exptype = cube.meta.exposure.type

            # Can't do outlier detection if there isn't a stack of images
            if len(cube.data.shape) < 3:
                self.log.warning("Input data are 2D; skipping outlier_detection")
                break

            self.log.info("Performing outlier detection on input images ...")
            cube = self.outlier_detection.run(cube)

            # Save crfints products
            if cube.meta.cal_step.outlier_detection == "COMPLETE":
                self.log.info("Saving crfints products with updated DQ arrays ...")
                # preserve output filename
                original_filename = cube.meta.filename

                # ensure output filename will not have duplicate asn_id
                if "_" + self.asn_id in original_filename:
                    original_filename = original_filename.replace("_" + self.asn_id, "")
                self.save_model(
                    cube, output_file=original_filename, suffix="crfints", asn_id=self.asn_id
                )
                cube.meta.filename = original_filename
            input_models[i] = cube

        # Create final photometry results as a single output
        # regardless of how many input members there may be
        phot_result_list = []

        # Imaging
        if input_exptype == "NRC_TSIMAGE" or (input_exptype == "MIR_IMAGE" and input_tsovisit):
            # Create name for extracted photometry (Level 3) product
            phot_tab_suffix = "phot"

            for cube in input_models:
                # Extract Photometry from imaging data
                phot_result_list.append(self.tso_photometry.run(cube))

        # Spectroscopy
        else:
            # Create name for extracted white-light (Level 3) product
            phot_tab_suffix = "whtlt"

            # Working with spectroscopic TSO data;
            # define output for x1d (level 3) products

            x1d_result = datamodels.TSOMultiSpecModel()
            x1d_result.update(input_models[0], only="PRIMARY")
            x1d_result.int_times = FITS_rec.from_columns(
                input_models[0].int_times.columns, nrows=input_models[0].meta.exposure.nints
            )

            # Remove source_type from the output model, if it exists, to prevent
            # the creation of an empty SCI extension just for that keyword.
            x1d_result.meta.target.source_type = None

            # For each exposure in the TSO...
            for cube in input_models:
                # interpolate pixels that have a NaN value or are flagged
                # as DO_NOT_USE or NON_SCIENCE.
                cube = self.pixel_replace.run(cube)
                state = cube.meta.cal_step.pixel_replace
                # Process spectroscopic TSO data
                # extract 1D
                self.log.info("Extracting 1-D spectra ...")
                result = self.extract_1d.run(cube)
                for row in cube.int_times:
                    # Subtract one to assign 1-indexed int_nums to int_times array locations
                    x1d_result.int_times[row[0] - 1] = row

                # SOSS F277W may return None - don't bother with that
                if (result is None) or (result.meta.cal_step.extract_1d == "SKIPPED"):
                    continue

                if input_exptype == "NIS_SOSS":
                    # SOSS data have yet to be photometrically calibrated
                    # Calibrate 1D spectra here.
                    result = self.photom.run(result)

                x1d_result.spec.extend(result.spec)

            # perform white-light photometry on all 1d extracted data
            self.log.info("Performing white-light photometry ...")
            phot_result_list.append(self.white_light.run(x1d_result))

            # Update some metadata from the association
            x1d_result.meta.asn.pool_name = input_models.asn_table["asn_pool"]
            x1d_result.meta.asn.table_name = Path(input_data).name

            # Save the final x1d Multispec model
            x1d_result.meta.cal_step.pixel_replace = state
            if len(x1d_result.spec) == 0:
                self.log.warning("extract_1d step could not be completed for any integrations")
                self.log.warning("x1dints products will not be created.")
            else:
                self.save_model(x1d_result, suffix="x1dints")

        # Done with all the inputs
        input_models.close()

        # Check for all null photometry results before saving
        all_none = np.all([(x is None) for x in phot_result_list])
        if all_none:
            self.log.warning("Could not create a photometric catalog; all results are null")
        else:
            # Otherwise, save results to a photometry catalog file
            phot_results = vstack(phot_result_list)
            phot_results.meta["number_of_integrations"] = len(phot_results)
            phot_tab_name = self.make_output_path(suffix=phot_tab_suffix, ext="ecsv")
            self.log.info(f"Writing Level 3 photometry catalog {phot_tab_name}")
            phot_results.write(phot_tab_name, format="ascii.ecsv", overwrite=True)

        # All done. Nothing to return, because all products have
        # been created here.
        return
