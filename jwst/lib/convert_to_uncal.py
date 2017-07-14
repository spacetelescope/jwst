from astropy.io import fits
from jwst import datamodels
import numpy as np

input_file = 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'
output_file = 'f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_uncal.fits'

# Load the input model containing the level-2b-like data
input = datamodels.ImageModel(input_file)

# Create and load 4-D data array for level-1b-like data
data4 = np.zeros((1,1,2048,2048))
data4[0,0] = input.data * input.meta.exposure.group_time

# Add dark ref data to science array
dark_file = '/grp/crds/cache/references/jwst/jwst_nirspec_dark_0026.fits'
dark = datamodels.DarkModel(dark_file)
data4[0,0] += dark.data[0,0]
dark.close()

# Zero-out reference pixel values so that refpix step has no effect
data4[0,0,:4,:] = 0.0
data4[0,0,-4:,:] = 0.0
data4[0,0,:,:4] = 0.0
data4[0,0,:,-4:] = 0.0

# Add superbias into science array
bias_file = '/grp/crds/cache/references/jwst/jwst_nirspec_superbias_0030.fits'
bias = datamodels.SuperBiasModel(bias_file)
data4[0,0] += bias.data
bias.close()

# Stuff the data and meta data into an output model
output = datamodels.Level1bModel(data=data4)
output.update(input)

# Set important meta data values
output.meta.exposure.ngroups = 1
output.meta.exposure.frame_divisor = 1
output.meta.wcsinfo.cdelt1 = input.meta.wcsinfo.cdelt1 / 3600.
output.meta.wcsinfo.cdelt2 = input.meta.wcsinfo.cdelt2 / 3600.
output.meta.wcsinfo.pc1_1 = -1.0
output.meta.wcsinfo.v2_ref = 378.770400
output.meta.wcsinfo.v3_ref = -428.155200
output.meta.wcsinfo.v3yangle = 138.492300
output.meta.wcsinfo.vparity = -1
output.meta.wcsinfo.ra_ref = output.meta.wcsinfo.crval1
output.meta.wcsinfo.dec_ref = output.meta.wcsinfo.crval2
output.meta.wcsinfo.roll_ref = 0.0

# Add GROUP info
output.group = np.ndarray(
    (1, ),
    dtype=[
        ('integration_number', '<i2'),
        ('group_number', '<i2'),
        ('end_day', '<i2'),
        ('end_milliseconds', '<i4'),
        ('end_submilliseconds', '<i2'),
        ('group_end_time', 'S26'),
        ('number_of_columns', '<i2'),
        ('number_of_rows', '<i2'),
        ('number_of_gaps', '<i2'),
        ('completion_code_number', '<i2'),
        ('completion_code_text', 'S36'),
        ('bary_end_time', '<f8'),
        ('helio_end_time', '<f8')
    ]
)
output.group[0]['integration_number'] = 1
output.group[0]['group_number'] = 1
output.group[0]['end_day'] = 1
output.group[0]['end_milliseconds'] = 1
output.group[0]['end_submilliseconds'] = 1
output.group[0]['group_end_time'] = '20170713t163230'
output.group[0]['number_of_columns'] = 1
output.group[0]['number_of_rows'] = 1
output.group[0]['number_of_gaps'] = 1
output.group[0]['completion_code_number'] = 1
output.group[0]['completion_code_text'] = 'COMPLETE'
output.group[0]['bary_end_time'] = 1.0
output.group[0]['helio_end_time'] = 1.0


# Save to an output file
output.save(output_file)
fits.setval(output_file, 'datamodl', value='RampModel')

input.close()
output.close()
