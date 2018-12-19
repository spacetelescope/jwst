# Module containing constants that were previously in separate modules
# have been moved here
#

import numpy as np

# Constants used in nrm_model.py
# piston array, in waves
phi_nb = np.array([0.028838669455909766, -0.061516214504502634, \
                    0.12390958557781348, -0.020389361461019516, \
                    0.016557347248600723, -0.03960017912525625, \
                    -0.04779984719154552])

# define phi at the center of F430M band:
phi_nb = phi_nb * 4.3e-6 # phi_nb in m

ctrs = np.array([[0.00000000, -2.640000],
                 [-2.2863100, 0.0000000],
                 [2.2863100, -1.3200001],
                 [-2.2863100, 1.3200001],
                 [-1.1431500, 1.9800000],
                 [2.2863100, 1.3200001],
                 [1.1431500, 1.9800000]])

# Constants in (but not used) webb_psf.py; included here in case
#    they will be used
m_ = 1.0
mm_ = m_ / 1000.0
um_ = mm_ / 1000.0
nm_ = um_ / 1000.0
