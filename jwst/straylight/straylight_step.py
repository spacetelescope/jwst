import logging

from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from jwst.straylight import straylight

__all__ = ["StraylightStep"]

log = logging.getLogger(__name__)


class StraylightStep(Step):
    """Correct for straylight caused by cross-artifact effect or residual cosmic rays."""

    class_alias = "straylight"

    spec = """
        clean_showers = boolean(default=False)  # Clean up straylight from residual cosmic ray showers
        shower_plane = integer(default=3)  # Throughput plane for identifying inter-slice regions
        shower_x_stddev = float(default=18) # X standard deviation for shower model
        shower_y_stddev = float(default=5) # Y standard deviation for shower model
        shower_low_reject = float(default=0.1) # Low percentile of pixels to reject
        shower_high_reject = float(default=99.9) # High percentile of pixels to reject
        save_shower_model = boolean(default=False) # Save the shower model
    """  # noqa: E501

    reference_file_types = ["mrsxartcorr", "regions"]

    def process(self, input_data):
        """
        Correct MIRI MRS data for the cross-artifact effect or residual cosmic rays.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.IFUImageModel`
            Input file name or datamodel to be corrected.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.IFUImageModel`
            Straylight corrected data.
        """
        output_model = self.prepare_output(input_data)

        # check the data is an IFUImageModel (not TSO)
        if isinstance(output_model, (datamodels.ImageModel, datamodels.IFUImageModel)):
            # Check for a valid mrsxartcorr reference file
            self.straylight_name = self.get_reference_file(output_model, "mrsxartcorr")

            if self.straylight_name == "N/A":
                log.warning("No MRSXARTCORR reference file found")
                log.warning("Straylight step will be skipped")
                output_model.meta.cal_step.straylight = "SKIPPED"
                return output_model

            log.info("Using mrsxartcorr reference file %s", self.straylight_name)

            with datamodels.MirMrsXArtCorrModel(self.straylight_name) as modelpars:
                # Apply the cross artifact correction
                output_model = straylight.correct_xartifact(output_model, modelpars)

            # Apply the cosmic ray droplets correction if desired
            if self.clean_showers:
                self.regions_name = self.get_reference_file(output_model, "regions")
                with datamodels.RegionsModel(self.regions_name) as f:
                    allregions = f.regions
                    output_model, output_shower_model = straylight.clean_showers(
                        output_model,
                        allregions,
                        self.shower_plane,
                        self.shower_x_stddev,
                        self.shower_y_stddev,
                        self.shower_low_reject,
                        self.shower_high_reject,
                        self.save_shower_model,
                    )
                    if self.save_shower_model and output_shower_model:
                        shower_path = self.make_output_path(
                            basepath=output_model.meta.filename, suffix="shower_model"
                        )
                        log.info(f"Saving shower model file {shower_path}")
                        output_shower_model.save(shower_path)
                        output_shower_model.close()

            output_model.meta.cal_step.straylight = "COMPLETE"

        else:
            if not isinstance(output_model, (datamodels.ImageModel, datamodels.IFUImageModel)):
                log.warning("Straylight correction not defined for datatype %s", output_model)
            log.warning("Straylight step will be skipped")
            output_model.meta.cal_step.straylight = "SKIPPED"

        return output_model
