from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer
import logging
from pathlib import Path

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class DataTypes:
    """Class to handle reading input data to cube_build."""

    template = {
        "asn_rule": "",
        "target": "",
        "asn_pool": "",
        "asn_type": "",
        "products": [{"name": "", "members": [{"exptype": "", "expname": ""}]}],
    }

    def __init__(self, input_data, single, output_file, output_dir):
        """
        Assemble input data for processing.

        Open the input data using datamodels and determine if data is
        a single input model, an association, or a set of input models
        contained in a ModelContainer. The method populates the self.input_models
        which is a list of input models. An initial base name for the output file
        is constructed.

        Parameters
        ----------
        input_data : str, IFUImageModel or ModelContainer
           Input data to cube_build either a filename, single model,
           association table, or a ModelContainer
        single : bool
           If True then creating single mode IFUCubes for outlier detection
           or mrs_matching. If false then creating standard IFUcubes
        output_file : str
           Optional user provided output file name.
        output_dir : str
           Optional user provided output directory name

        Raises
        ------
        NotIFUImageModelError
           IFU data was not the input data
        TypeError
           Input data was not of correct form
        """
        self.input_models = []
        self.output_name = None

        # open the input_data with datamodels
        # if input_data is filename or model when it is opened it is a model
        # if input_data is an association name or ModelContainer then it is opened as a container

        input_models = datamodels.open(input_data)
        # if input_data is a filename, we will need to close the opened file
        self._opened = [input_models]

        if isinstance(input_models, datamodels.IFUImageModel):
            # It's a single image that's been passed in as a model
            # input_data is a model
            filename = input_models.meta.filename
            self.input_models.append(input_models)
            self.output_name = self.build_product_name(filename)

        elif isinstance(input_models, ModelContainer):
            self.output_name = "Temp"
            self.input_models = input_models
            if not single:  # find the name of the output file from the association
                self.output_name = input_models.asn_table["products"][0]["name"]
        else:
            # close files opened above
            self.close()
            raise TypeError(f"Failed to process file type {type(input_models)}")
        # If the user has set the output name, strip off *.fits.
        # Suffixes will be added to this name later, to designate the
        # channel+subchannel (MIRI MRS) or grating+filter (NRS IFU) the output cube covers.

        if output_file is not None:
            self.output_name = Path(output_file).stem

        if output_dir is not None:
            self.output_name = output_dir + "/" + self.output_name

    def close(self):
        """Close any files opened by this instance."""
        [f.close() for f in self._opened]

    # _______________________________________________________________________________
    def build_product_name(self, filename):
        """
        Determine the base of output name if an input data is a FITS filename.

        Parameters
        ----------
        filename : str
          If a string filename is given as the input data to cube_build, then
          determine the base name of the output IFU Cube filename.

        Returns
        -------
        single_product : str
          Output base filename.
        """
        indx = filename.rfind(".fits")
        indx_try = filename.rfind("_rate.fits")
        indx_try2 = filename.rfind("_cal.fits")

        if indx_try > 0:
            single_product = filename[:indx_try]
        elif indx_try2 > 0:
            single_product = filename[:indx_try2]
        else:
            single_product = filename[:indx]
        return single_product


# _______________________________________________________________________________


class NotIFUImageModelError(Exception):
    """Raise Exception if data is not of type IFUImageModel."""

    pass
