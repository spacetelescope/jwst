import logging
import numpy as np

from stdatamodels.jwst import datamodels

from jwst.persistence import persistence
from jwst.stpipe import Step

__all__ = ["PersistenceStep"]

log = logging.getLogger(__name__)


class PersistenceStep(Step):
    """Correct a science image for persistence."""

    class_alias = "persistence"

    spec = """
        input_trapsfilled = string(default="") # Name of the most recent trapsfilled file for the current detector
        flag_pers_cutoff = float(default=40.) # Pixels with persistence correction >= this value in DN will be flagged in the DQ
        save_persistence = boolean(default=False) # Save subtracted persistence to an output file with suffix '_output_pers'
        save_trapsfilled = boolean(default=True) # Save updated trapsfilled file with suffix '_trapsfilled'
        modify_input = boolean(default=False)
        persistence_time = integer(default=None) # Time, in seconds, to use for persistence window
        persistence_array = list(default=None) # A 2-D array or none.
        persistence_dnu = boolean(default=False) # If True the set the DO_NOT_USE flag with PERSISTENCE
    """  # noqa: E501

    reference_file_types = ["trapdensity", "trappars", "persat"]

    def process(self, step_input):
        """
        Execute the persistence correction step.

        Parameters
        ----------
        step_input : `~stdatamodels.jwst.datamodels.RampModel` or str
            Input datamodel or file to be corrected

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.RampModel`
            The persistence corrected datamodel
        """
        if self.input_trapsfilled is not None:
            if (self.input_trapsfilled == "None") or (len(self.input_trapsfilled) == 0):
                self.input_trapsfilled = None

        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)

        self.process_persistence_options(result)

        trap_density_filename = self.get_reference_file(result, "trapdensity")
        trappars_filename = self.get_reference_file(result, "trappars")
        persat_filename = self.get_reference_file(result, "persat")

        # Is any reference file missing?
        missing = False
        missing_reftypes = []
        if persat_filename == "N/A":
            missing = True
            missing_reftypes.append("PERSAT")
        if trap_density_filename == "N/A":
            missing = True
            missing_reftypes.append("TRAPDENSITY")
        if trappars_filename == "N/A":
            missing = True
            missing_reftypes.append("TRAPPARS")
        if missing:
            if len(missing_reftypes) == 1:
                msg = "Missing reference file type:  " + missing_reftypes[0]
            else:
                msg = "Missing reference file types: "
                for name in missing_reftypes:
                    msg += " " + name
            log.warning("%s", msg)
            result.meta.cal_step.persistence = "SKIPPED"
            return result, None

        if self.input_trapsfilled is None:
            traps_filled_model = None
        else:
            traps_filled_model = datamodels.TrapsFilledModel(self.input_trapsfilled)
        trap_density_model = datamodels.TrapDensityModel(trap_density_filename)
        trappars_model = datamodels.TrapParsModel(trappars_filename)
        persat_model = datamodels.PersistenceSatModel(persat_filename)

        pers_a = persistence.DataSet(
            result,
            traps_filled_model,
            self.flag_pers_cutoff,
            self.save_persistence,
            trap_density_model,
            trappars_model,
            persat_model,
            self.persistence_time,
            self.persistence_array,
            self.persistence_dnu,
        )
        (result, traps_filled, output_pers, skipped) = pers_a.do_all()
        if skipped:
            result.meta.cal_step.persistence = "SKIPPED"
        else:
            result.meta.cal_step.persistence = "COMPLETE"

        if traps_filled_model is not None:  # input traps_filled
            del traps_filled_model
        if traps_filled is not None:  # output traps_filled
            # Save the traps_filled image with suffix 'trapsfilled'.
            self.save_model(traps_filled, "trapsfilled", force=self.save_trapsfilled)
            del traps_filled

        if output_pers is not None:  # output file of persistence
            self.save_model(output_pers, suffix="output_pers", force=self.save_persistence)
            del output_pers

        persistence_list = None
        if pers_a.persistence_array is not None:
            persistence_list = self.persistence_array.tolist()

        # Cleanup
        del trap_density_model
        del trappars_model
        del persat_model

        # XXX Instead of returning persistence_list, could store it in ASDF.
        return result, persistence_list

    def process_persistence_options(self, result):
        """
        Processing  persistence_time, persistence_array, and persistence_dnu as the inputs.

        Parameters
        ----------
        result : RampModel
            The RampModel on which to process the persistence flag.
        """
        # Could make less than or equal to frametime.
        if self.persistence_time is None or self.persistence_time <= 0.0:
            self.persistence_time = None
            # XXX raise error or log info?
            return  # No persistence option chosen

        # XXX think about using ASDF for persistence_array input and output
        _, _, nrows, ncols = result.groupdq.shape
        if self.persistence_array is not None:
            self.persistence_array_create = False 
            self.persistence_array = np.array(self.persistence_array)

            # Make sure array has correct dimensions
            dims = self.persistence_array.shape 
            if len(dims) != 2 or dims[0] != nrows or dims[1] != ncols:
                raise ValueError("'persistence_array' needs to be a 2-D list with dimensions (nrows, ncols)")
        else:
            self.persistence_array_create = True
            self.persistence_array = np.zeros(shape=(nrows, ncols), dtype=np.float64)
