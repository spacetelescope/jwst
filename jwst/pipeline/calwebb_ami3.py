#!/usr/bin/env python
import logging
from pathlib import Path

from jwst.stpipe import Pipeline

# step imports
from jwst.ami import ami_analyze_step
from jwst.ami import ami_normalize_step

__all__ = ["Ami3Pipeline"]

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class Ami3Pipeline(Pipeline):
    """
    Apply all level-3 calibration steps to an association of level-2b AMI exposures.

    Included steps are:
    ami_analyze (fringe detection)
    ami_normalize (normalize results by reference target)
    """

    class_alias = "calwebb_ami3"

    # Define aliases to steps
    step_defs = {
        "ami_analyze": ami_analyze_step.AmiAnalyzeStep,
        # 'ami_average': ami_average_step.AmiAverageStep,
        "ami_normalize": ami_normalize_step.AmiNormalizeStep,
    }

    def process(self, input_data):
        """
        Run the ami3 pipeline on the input association.

        Parameters
        ----------
        input_data : str, Association, or any other type loadable as an association
            The input association to be processed. If a string, a filename
            pointing to an asn.json file is assumed.
        """
        log.info("Starting calwebb_ami3")

        # Load the input association table
        asn = self.load_as_level3_asn(input_data)

        # We assume there's one final product defined by the association
        asn_id = asn["asn_id"]
        prod = asn["products"][0]
        self.output_file = prod.get("name", self.output_file)

        # Construct lists of all the PSF and science target members
        # so that we know what we've been given to work with
        psf_files = [m["expname"] for m in prod["members"] if m["exptype"].upper() == "PSF"]
        targ_files = [m["expname"] for m in prod["members"] if m["exptype"].upper() == "SCIENCE"]

        # Make sure we found some science target members
        if len(targ_files) == 0:
            log.error("No science target members found in association table")
            log.error("Calwebb_ami3 processing will be aborted")
            return

        # If there aren't any PSF images, we can't do the normalize step
        if len(psf_files) == 0:
            log.info("No PSF reference members found in association table")
            log.info("ami_normalize step will be skipped")

        # Run ami_analyze on all the target members
        targ_lg = []
        for input_file in targ_files:
            # Do the LG analysis for this image
            log.debug("Do LG processing for member %s", input_file)
            result1, result2, result3 = self.ami_analyze.run(input_file)

            # Save the averaged LG analysis results to a file
            result1.meta.asn.pool_name = asn["asn_pool"]
            result1.meta.asn.table_name = Path(asn.filename).name
            self.save_model(result1, output_file=input_file, suffix="ami-oi", asn_id=asn_id)

            # Save the result for use as input to ami_average
            targ_lg.append(result1)

        # Run ami_analyze on all the psf members
        psf_lg = []
        for input_file in psf_files:
            # Do the LG analysis for this image
            log.debug("Do LG processing for member %s", input_file)
            result1, result2, result3 = self.ami_analyze.run(input_file)

            # Save the LG analysis results to a file
            result1.meta.asn.pool_name = asn["asn_pool"]
            result1.meta.asn.table_name = Path(asn.filename).name
            self.save_model(result1, output_file=input_file, suffix="psf-ami-oi", asn_id=asn_id)

            # Save the result for use as input to ami_average
            psf_lg.append(result1)

        # This zip operation matches science exposures to reference star exposures
        # in a one-to-one fashion, truncating a list if it is longer than the
        # other. This is not intended behavior and should be fixed by JP-3978.
        if (len(psf_files) > 0) & (len(targ_files) > 0):
            for targ, psf in zip(targ_lg, psf_lg, strict=False):
                result = self.ami_normalize.run(targ, psf)
                # Save the result
                result.meta.asn.pool_name = asn["asn_pool"]
                result.meta.asn.table_name = Path(asn.filename).name

                # Perform blending of metadata for all inputs to this output file
                # self.log.info('Blending metadata for PSF normalized target')
                self.save_model(result, suffix="aminorm-oi")
                result.close()
            del psf_lg
            del targ_lg

        # We're done
        log.info("... ending calwebb_ami3")

        return
