#
#  Module for appyling emission correction
#

import logging
from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class DataSet():
    """
    Input dataset for emission subtraction

    Parameters
    ----------

    """

    def __init__(self, input_dm):

        """
        Short Summary
        -------------
        Set file name of input file and its DM object

        Parameters
        ----------
        input_dm: data model object
            input Data Model object

        """

        try:
            model = datamodels.open(input_dm)
        except Exception:
            log.info('Error opening: %s ', input_dm)
            model = None

        self.input = model
        self.input_file = input_dm


    def do_all(self):
        """
        Short Summary
        -------------
        Execute all tasks for Emission Correction - which
            is a no-op before Build 1


        Parameters
        ----------

        Returns
        -------
        self.input_file: fits file
            emission-corrected input file data

        """
        self.apply_emission(get_em_file_name())

        return self.input


    def apply_emission(self, emission):
        """
        Short Summary
        -------------
        Emission Correction: may eventually subtract an emission image, but
            until Build 1 this will be a no-op


        Parameters
        ----------
        emission: emission object
            instance of object

        Returns
        -------

        """
        log.info("EmissionCorrection (no-op) applied")


def get_em_file_name():

    """
    Short Summary
    -------------
    Retrieve the particular emission reference file name. Will be done
        by call to CRDS eventually

    Parameters
    ----------
    input_obj: input
        File name or Data Model object

    Returns
    -------
    em_file: string
        name of emission file

    """

    em_file = None

    log.info('Will attempt to use emission file: %s', em_file)
    return em_file
