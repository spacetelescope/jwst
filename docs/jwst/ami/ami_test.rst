.. _ami_unit_test:

AMI unit tests
==============

There are unit tests for AMI Analyze and AMI interface. 

test_ami_interface
------------------

 - Make sure ami_analyze fails if input is an input file of type  _calints
 - Make sure ami_analyze fails if input is CubeModel for _calints
 - Make sure that ami_analyze fails if no throughput reffile is available


test_ami_analyze
----------------

utils module tests:
+++++++++++++++++++

For the module *utils* we have several tests that compare the calculated value with a known value. The tests are:

 - Test of rebin() and krebin() 
 - Test of quadratic 
 - Test of findmax 
 - Test of makeA
 - Test of fringes2pistons
 - Test of rcrosscorrelate()
 - Test of crosscorrelate()

leastsqnrm module tests:
________________________

 - Test of rotatevectors()
        Positive x decreases under slight rotation, and positive y
        increases under slight rotation.
 - Test of flip()
        Change sign of 2nd coordinate of holes.
 - Test of mas2rad()
        Convert angle in milli arc-sec to radians.
 - Test of rad2mas()
        Convert input angle in radians to milli arc sec.
 - Test of sin2deltapistons()
        Each baseline has one sine and one cosine fringe with a coefficient
        that depends on the piston difference between the two holes that make
        the baseline.  For a 7-hole mask there are 21 baselines and therefore
        there are 42 sine and cosine terms that contribute to the fringe model.
        This function calculates the sine of this piston difference.
 - Test of cos2deltapistons()
        Each baseline has one sine and one cosine fringe with a coefficient
        that depends on the piston difference between the two holes that make
        the baseline.  For a 7-hole mask there are 21 baselines and therefore
        there are 42 sine and cosine terms that contribute to the fringe model.
        This function calculate the cosine of this piston difference.
 - Test of replacenan()
        Replace singularities encountered in the analytical hexagon Fourier
        transform with the analytically derived limits. (pi/4)
 - Test of hexpb()
        Calculate the primary beam for hexagonal holes.
 - Test of model_array
        Create a model using the specified wavelength.
 - Test of ffc
        Calculate cosine terms of analytic model.
 - Test of ffs
        Calculate sine terms of analytic model.
 - Test of return_CAs
        Calculate the closure amplitudes.
 -  Test of closurephase
         Calculate closure phases between each pair of holes.
 - Test of redundant_cps
        Calculate closure phases for each set of 3 holes.
 - Test of populate_symmamparray
        Populate the symmetric fringe amplitude array.
 - Test of populate_antisymmphasearray
        Populate the antisymmetric fringe phase array.
 - Test of tan2visibilities
        From the solution to the fit, calculate the fringe amplitude and phase.
 - Test of multiplyenv
        Multiply the envelope by each fringe 'image'.

hexee module tests:
+++++++++++++++++++
 -  Test of g_eeAG()
        Calculate the Fourier transform of one half of a hexagon that is
        bisected from one corner to its diametrically opposite corner.
 -  Test of glimit()
        Calculate the analytic limit of the Fourier transform of one half of the
        hexagon along eta=0.

analyticnrm2 module tests:
++++++++++++++++++++++++++

 -  Test of PSF() 
 -  Test of ASFhex() in the analyticnrm2 module FOR HEX 
 -  Test of interf() 
 -  Test of phasor()
 
webb_psf module test:
+++++++++++++++++++++

 -  Test of PSF()
        Create a Throughput datamodel, having a dummy filter bandpass data
        that peaks at 1.0 at the center and decreases in the wings.

