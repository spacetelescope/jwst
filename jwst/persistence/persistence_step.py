import asdf
import datetime
import logging
import numpy as np
import os

from stdatamodels.jwst import datamodels

from jwst.lib.suffix import KNOW_SUFFIXES
from jwst.persistence import persistence
from jwst.stpipe import Step

__all__ = ["PersistenceStep"]

log = logging.getLogger(__name__)


class PersistenceStep(Step):
    """Correct a science image for persistence."""

    class_alias = "persistence"

    spec = """
        save_persistence = boolean(default=False) # Save subtracted persistence to an output file with suffix '_output_pers'
        persistence_time = integer(default=None) # Time, in seconds, to use for persistence window
        persistence_array_file = string(default=None) # A path to an ASDF file containing a 2-D array of persistence times per pixel
        persistence_dnu = boolean(default=False) # If True the set the DO_NOT_USE flag with PERSISTENCE
        skip = boolean(default=False) # Skip the persistence step entirely
    """  # noqa: E501

    def process(self, step_input):
        """
        Execute the persistence correction step.

        Parameters
        ----------
        step_input : `~stdatamodels.jwst.datamodels.RampModel` or str
            Input datamodel or file to be corrected

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.RampModel`
            The persistence corrected datamodel
        """
        result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)
        if self.skip:
            log.info("Skipping persistence step as requested.")
            result.meta.cal_step.persistence = "SKIPPED"
            return result

        self.process_persistence_options(result)


        pers_a = persistence.DataSet(
            result,
            self.save_persistence,
            self.persistence_time,
            self.persistence_array,
            self.persistence_dnu,
        )
        result, skipped = pers_a.do_all()

        if skipped:
            result.meta.cal_step.persistence = "SKIPPED"
        else:
            result.meta.cal_step.persistence = "COMPLETE"

        if pers_a.save_persistence:
            self.write_persistence_array(result)

        return result

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
            self.persistence_array = None
            return  # No persistence option chosen

        _, _, nrows, ncols = result.groupdq.shape
        if self.persistence_array_file is not None:
            self.persistence_array_create = False 

            with asdf.open(self.persistence_array_file) as af:
                self.persistence_array = af.tree["persistence_data"].copy()

            # Make sure array has correct dimensions
            dims = self.persistence_array.shape 
            if len(dims) != 2 or dims[0] != nrows or dims[1] != ncols:
                raise ValueError("'persistence_array' needs to be a 2-D list with dimensions (nrows, ncols)")
        else:
            self.persistence_array_create = True
            self.persistence_array = np.zeros(shape=(nrows, ncols), dtype=np.float64)

    def write_persistence_array(self, result):
        """
        Write the persistence array to an ASDF file.

        Parameters
        ----------
        result : RampModel
            The RampModel on which to process the persistence flag.
        """
        # Setup persistence array filename with time suffix to avoid overwriting existing files
        now = datetime.datetime.now()
        time_fmt = "%Y%m%d%H%M%S%f"
        time_str = now.strftime(time_fmt)
        pers_suffix = f"pers{time_str}"

        # persistence_array_file always gets set if the persistence options are processed.
        if self.persistence_array_file is None:
            filename = result.meta.filename
        else:
            filename = self.persistence_array_file

        filename = self.make_output_path(basepath=filename, suffix=pers_suffix, ext="asdf")

        # Write persistence array to ASDF file
        tree = {"persistence_data": self.persistence_array}
        with asdf.AsdfFile(tree) as af:
            af.write_to(filename)
