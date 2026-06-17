import logging
from pathlib import Path

import asdf
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
        save_persistence = string(default=None) # Name of ASDF output file to save the persistence array
        persistence_time = integer(default=None) # Time, in seconds, to use for persistence window
        persistence_array_file = string(default=None) # A path to an ASDF file containing a 2-D array of persistence times per pixel
        persistence_dnu = boolean(default=False) # If True the set the DO_NOT_USE flag with PERSISTENCE
        skip = boolean(default=True) # By default, skip the step.
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

        if self.process_persistence_options(result) == "Failed":
            log.info("Persistence step failed due to invalid persistence_time.")
            result.meta.cal_step.persistence = "FAILED"
            return result

        pers_a = persistence.DataSet(
            result,
            self.save_persistence,
            self.persistence_time,
            self.persistence_array,
            self.persistence_dnu,
        )
        result, skipped = pers_a.do_all()

        result.meta.cal_step.persistence = "COMPLETE"
        if pers_a.save_persistence is not None:
            self.write_persistence_array(result, pers_a.save_persistence)

        return result

    def process_persistence_options(self, result):
        """
        Process persistence_time, persistence_array, and persistence_dnu as the inputs.

        Parameters
        ----------
        result : RampModel
            The RampModel on which to process the persistence flag.

        Returns
        -------
        ret : str or None
            "Failed" if invalid persistence_time; otherwise NoneType.
        """
        # Could make less than or equal to frametime.
        if self.persistence_time is None or self.persistence_time <= 0.0:
            self.persistence_time = None
            self.persistence_array = None
            ret = "Failed"
            return ret

        _, _, nrows, ncols = result.groupdq.shape
        if self.persistence_array_file is not None:
            self.get_persistence_array_from_file(nrows, ncols)
        else:
            self.persistence_array = np.zeros(shape=(nrows, ncols), dtype=np.float64)

        return None

    def write_persistence_array(self, result, filename):
        """
        Write the persistence array to an ASDF file.

        Parameters
        ----------
        result : RampModel
            The RampModel on which to process the persistence flag.
        """
        ext = str(Path(filename).suffix)
        stem = Path(filename).stem
        parent = Path(filename).parent
        root = str(parent / stem)
        if ext != ".asdf":
            filename = f"{root}.asdf"

        # Write persistence array to ASDF file
        rows, cols = np.nonzero(self.persistence_array)
        vals = self.persistence_array[rows, cols]
        tree = {
            "filename": result.meta.filename,
            "rows": rows,
            "cols": cols,
            "vals": vals,
            "pers_time": self.persistence_time,
        }

        with asdf.AsdfFile(tree) as af:
            af.write_to(filename)

    def get_persistence_array_from_file(self, nrows, ncols):
        """
        Get the persistence array from an ASDF file.

        Parameters
        ----------
        nrows : int
            The number of rows in the RampModel data.

        ncols : int
            The number of columns in the RampModel data.
        """
        with asdf.open(self.persistence_array_file) as pers_file:
            if pers_file["pers_time"] != self.persistence_time:
                raise ValueError("Invalid persistence file. Mismatch of persistence time.")

            rows = pers_file["rows"]
            cols = pers_file["cols"]
            vals = pers_file["vals"]

            self.persistence_array = np.zeros((nrows, ncols), dtype=vals.dtype)
            self.persistence_array[rows, cols] = vals
