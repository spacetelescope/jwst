# Routines used for building cubes
import os
from .. import datamodels
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#********************************************************************************
class DataTypes():
#********************************************************************************

    """
    Class to handle reading the input to the processing, which
    can be a single science exposure or an IFU cube association table.
    The input and output member info is loaded into an ASN table model.
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

    def __init__(self, input,single,output_file,output_dir):

        self.input_models = []
        self.filenames = []
        self.output_name = None

        # open the input with datamodels
        # if input is filename or model when it is openned it is a model
        # if input if an assocation name or ModelContainer then it is openned as a container
#        print('***input type***',type(input))
        input_try = datamodels.open(input)
#        print('input_try',type(input_try))

        if isinstance(input_try, datamodels.IFUImageModel):
#            print('this is a single file or Model ')
            # It's a single image that's been passed in as a model
            # input is a model
            self.filenames.append(input_try.meta.filename)
            self.input_models.append(input_try)
            self.output_name = self.build_product_name(self.filenames[0])

        elif isinstance(input_try, datamodels.ModelContainer):
#            print('this is a model container type or association read in a ModelContainer')
            self.output_name  = 'Temp'
            if not single:  # find the name of the output file from the association
                with datamodels.ModelContainer(input_try) as input_model:
                    self.output_name =input_model.meta.asn_table.products[0].name
            for i in range(len(input_try)):
                # check if input data is an IFUImageModel
                if not  isinstance(input_try[i], datamodels.IFUImageModel):
                    serror = str(type(input_try[i]))
                    raise NotIFUImageModel("Input data is not a IFUImageModel, instead it is %s",serror)

                model = datamodels.IFUImageModel(input_try[i])
                self.input_models.append(model)
                self.filenames.append(model.meta.filename)

        else:
            raise TypeError

# if the user has set the output name - strip out *.fits
# later suffixes will be added to this name to designate the
# channel, subchannel or grating,filter the data is covers.

        if output_file !=None :
            basename,ext = os.path.splitext(os.path.basename(output_file))
            self.output_name = basename


        if output_dir !=None :
            self.output_name= output_dir + '/' + self.output_name

#        print('*****************',self.output_name)

    def build_product_name(self, filename):
        indx = filename.rfind('.fits')
        indx_try = filename.rfind('_rate.fits') # standard expected filename in CalSpec2
        indx_try2 = filename.rfind('_cal.fits') # standard expected filename


        if indx_try > 0:
            single_product = filename[:indx_try]
        elif indx_try2 > 0:
            single_product = filename[:indx_try2]
        else:
            single_product = filename[:indx]
        return single_product



# Raise Exceptions 
class NotIFUImageModel(Exception):
    pass
