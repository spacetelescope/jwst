#!/usr/bin/env python
from collections import defaultdict
from pathlib import Path

from stdatamodels.jwst import datamodels

from jwst.associations.load_as_asn import LoadAsLevel2Asn
from jwst.stpipe import Pipeline

# calwebb IMAGE2 step imports
from jwst.background import background_step
from jwst.assign_wcs import assign_wcs_step
from jwst.flatfield import flat_field_step
from jwst.photom import photom_step
from jwst.resample import resample_step


__all__ = ["Image2Pipeline"]


class Image2Pipeline(Pipeline):
    """
    Process JWST imaging-mode slope data from Level-2a to Level-2b.

    Included steps are:
    bkg_subtract, assign_wcs, flat_field, photom and resample.
    """

    class_alias = "calwebb_image2"

    spec = """
        save_bsub = boolean(default=False) # Save background-subtracted science
    """  # noqa: E501

    # Define alias to steps
    step_defs = {
        "bkg_subtract": background_step.BackgroundStep,
        "assign_wcs": assign_wcs_step.AssignWcsStep,
        "flat_field": flat_field_step.FlatFieldStep,
        "photom": photom_step.PhotomStep,
        "resample": resample_step.ResampleStep,
    }

    # List of normal imaging exp_types
    image_exptypes = ["MIR_IMAGE", "NRC_IMAGE", "NIS_IMAGE", "FGS_IMAGE"]

    def process(self, input_data):
        """
        Run the Image2Pipeline on the input data.

        Parameters
        ----------
        input_data : str, Level2 Association, or DataModel
            The exposure or association of exposures to process.

        Returns
        -------
        list[JWSTDataModel]
            The calibrated data models.
        """
        self.log.info("Starting calwebb_image2 ...")

        # Retrieve the input(s)
        asn = LoadAsLevel2Asn.load(input_data, basename=self.output_file)
        if len(asn["products"]) > 1 and self.output_file is not None:
            self.log.warning(
                "Multiple products in input association. Output file name will be ignored."
            )
            self.output_file = None

        # Each exposure is a product in the association.
        # Process each exposure.
        results = []
        for product in asn["products"]:
            self.log.info("Processing product {}".format(product["name"]))
            if (self.save_results) & (self.output_file is None):
                self.output_file = product["name"]
            if not hasattr(asn, "filename"):
                asn.filename = "singleton"

            result = self.process_exposure_product(
                product, asn["asn_pool"], Path(asn.filename).name
            )

            # Save result
            suffix = "cal"
            if isinstance(result, datamodels.CubeModel):
                suffix = "calints"
            result.meta.filename = self.make_output_path(basepath=self.output_file, suffix=suffix)
            results.append(result)
            self.output_file = None

        self.log.info("... ending calwebb_image2")

        self.output_use_model = True
        self.suffix = False
        return results

    # Process each exposure
    def process_exposure_product(
        self,
        exp_product,
        pool_name=" ",
        asn_file=" ",
    ):
        """
        Calibrate an exposure found in the association product.

        Parameters
        ----------
        exp_product : dict
            A Level2b association product.
        pool_name : str
            The pool file name. Used for recording purposes only.
        asn_file : str
            The name of the association file. Used for recording purposes only.

        Returns
        -------
        JWSTDataModel
            The final calibrated product.
        """
        # Find all the member types in the product
        members_by_type = defaultdict(list)
        for member in exp_product["members"]:
            members_by_type[member["exptype"].lower()].append(member["expname"])

        # Get the science member. Technically there should only be
        # one. We'll just get the first one found.
        science = members_by_type["science"]
        if len(science) != 1:
            self.log.warning(
                "Wrong number of science files found in {}".format(exp_product["name"])
            )
            self.log.warning("    Using only first one.")
        science = science[0]

        self.log.info("Working on input %s ...", science)
        if isinstance(science, datamodels.JwstDataModel):
            input_data = science
        else:
            input_data = datamodels.open(science)

        # Record ASN pool and table names in output
        input_data.meta.asn.pool_name = pool_name
        input_data.meta.asn.table_name = asn_file
        input_data.meta.filename = self.make_output_path(basepath=self.output_file)

        # Do background processing, if necessary
        if len(members_by_type["background"]) > 0:
            # Setup for saving
            self.bkg_subtract.suffix = "bsub"
            if isinstance(input_data, datamodels.CubeModel):
                self.bkg_subtract.suffix = "bsubints"

            # Backwards compatibility
            if self.save_bsub:
                self.bkg_subtract.save_results = True

            # Call the background subtraction step
            input_data = self.bkg_subtract.run(input_data, members_by_type["background"])

        # work on slope images
        input_data = self.assign_wcs.run(input_data)
        input_data = self.flat_field.run(input_data)
        input_data = self.photom.run(input_data)

        # Resample individual exposures, but only if it's one of the
        # regular 2D science image types
        if (
            input_data.meta.exposure.type.upper() in self.image_exptypes
            and len(input_data.data.shape) == 2
        ):
            self.resample.save_results = self.save_results
            self.resample.suffix = "i2d"
            self.resample.run(input_data)

        # That's all folks
        self.log.info("Finished processing product {}".format(exp_product["name"]))
        return input_data
