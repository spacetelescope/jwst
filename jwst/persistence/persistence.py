#
#  Module for applying persistence files


import numpy as np
import logging
from jwst import datamodels

log = logging.getLogger()
log.setLevel(logging.DEBUG)


class DataSet():
    """
    Input dataset to which persisence will be applied

    Parameters
   ----------
    """
    def __init__(self, input_DM):
        """
        Short Summary
        -------------
        Set file name of persistence file

        Parameters
        ----------
        input_DM: data model object
            input Data Model object

        """
        try:
            self.input_file = input_DM
            model = models.open(input_DM)
            # If model comes back as generic DataModel, reopen as MultiSlit
            if isinstance(model, models.CubeModel) or isinstance(model, models.ImageModel):
                pass
            elif isinstance(model, models.DataModel):
                model.close()
                model = models.MultiSlitModel(input_DM)
            self.input = model
        except Exception as errmess:
            log.error('Error opening %s', input_DM)
            self.input = None


    def do_all(self):
        """
        Short Summary
        -------------
        Execute all tasks for Persistence Correction - which
            is a no-op before Build 1

        Parameters
        ----------

        Returns
        -------
        self.input: input fits file
            persistence-applied input file data

        """
        self.apply_persistence(get_pers_file_name())

        return self.input


    def apply_persistence(self, persistence):
        """
        Short Summary
        -------------
        Persistence Correction: may eventually apply persistence, but until
            Build 1 will be a no-op

        Parameters
        ----------
        persistence: persistence object
            instance of object

        Returns
        -------

        """
        pass


def get_pers_file_name():
    """
    Short Summary
    -------------
    Retrieve the particular persistence reference file name. Will be
    done via call(s) to CRDS eventually

    Parameters
    ----------

    Returns
    -------
    pers_file: string
        name of persistence file (None for pre Build 1)

    """

    pers_file = None
    return pers_file
