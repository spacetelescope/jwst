""" Class Data Type is used to read in the input data. It also determines
if the input data is a single science exposure, an association table, a
single datamodel or several data models stored in a ModelContainer.
"""
import os
from .. import datamodels
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# ******************************************************************************
class DataTypes():

    """ Class to handle reading input data to cube_build.
    """

    template = {"asn_rule": "",
                "target": "",
                "asn_pool": "",
                "asn_type": "",
                "products": [
                   {"name": "",
                    "members": [
                      {"exptype": "",
                       "expname": ""}
                      ]
                    }
                   ]
                }

    def __init__(self, input, single, output_file, output_dir):

        """ Read in input data and determine what type of input data.

        Open the input data using datamodels and determine if data is
        a single input model, an association, or a set of input models
        contained in a ModelContainer.

        Parameters
        ----------
        input : datamodel or  ModelContainter
           Input data to cube_build either a filename, single model,
           association table, or a ModelContainter
        single : boolean
           If True then creating single mode IFUCubes for outlier detection
           or mrs_matching. If false then creating standard IFUcubes
        output_file : str
           Optional user provided output file name.
        output_dir : str
           Optional user provided output directory name

        Notes
        -----
        The method populates the self.input_models which is list of input models.
        An initial base name for the output file is constructed.

        Raises
        ------
        NotIFUImageModel
           IFU data was not the input data
        TypeError
           Input data was not of correct form

        """

        self.input_models = []
        self.filenames = []
        self.output_name = None

        # open the input with datamodels
        # if input is filename or model when it is openned it is a model
        # if input if an assocation name or ModelContainer then it is openned as a container
        # print('***input type***',type(input))
        input_try = datamodels.open(input)

        if isinstance(input_try, datamodels.IFUImageModel):
            # print('this is a single file or Model ')
            # It's a single image that's been passed in as a model
            # input is a model
            self.filenames.append(input_try.meta.filename)
            self.input_models.append(input_try)
            self.output_name = self.build_product_name(self.filenames[0])

        elif isinstance(input_try, datamodels.ModelContainer):
            self.output_name = 'Temp'
            if not single:  # find the name of the output file from the association
                self.output_name = input_try.meta.asn_table.products[0].name
            for model in input_try:
                # check if input data is an IFUImageModel
                if not isinstance(model, datamodels.IFUImageModel):
                    raise NotIFUImageModel(
                        f"Input data is not a IFUImageModel, instead it is {model}")
                self.filenames.append(model.meta.filename)
            self.input_models = input_try

        else:
            raise TypeError("Failed to process file type {}".format(type(input_try)))

# if the user has set the output name - strip out *.fits
# later suffixes will be added to this name to designate the
# channel, subchannel or grating,filter the data is covers.

        if output_file is not None:
            basename, ext = os.path.splitext(os.path.basename(output_file))
            self.output_name = basename

        if output_dir is not None:
            self.output_name = output_dir + '/' + self.output_name

#        print('*****************',self.output_name)

# _______________________________________________________________________________
    def build_product_name(self, filename):
        """ Determine the base of output name if an input data is a fits filename.

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

        indx = filename.rfind('.fits')
        indx_try = filename.rfind('_rate.fits')
        indx_try2 = filename.rfind('_cal.fits')

        if indx_try > 0:
            single_product = filename[:indx_try]
        elif indx_try2 > 0:
            single_product = filename[:indx_try2]
        else:
            single_product = filename[:indx]
        return single_product

# _______________________________________________________________________________


class NotIFUImageModel(Exception):
    """ Raise Exception if data is not of type IFUImageModel
    """
    pass
