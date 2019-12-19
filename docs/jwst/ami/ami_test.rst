.. _ami_unit_test:

ami unit tests
==============

There are unit test for AMI Analyze and AMI interfarce. 

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

 - Test of rebin() and krebin() in utils module
 - Test of quadratic in utils module
 - Test of findmax in utils module
 - Test of makeA in utils module
 - Test of fringes2pistons in utils module
 - Test of rcrosscorrelate() in utils module
 - Test of crosscorrelate() in utils module

leastsqnrm module tests:
________________________

 - Test of rotatevectors() in leastsqnrm module.
        Positive x decreases under slight rotation, and positive y
        increases under slight rotation.
 - Test of flip() in leastsqnrm module.
        Change sign of 2nd coordinate of holes.
 - Test of mas2rad() in leastsqnrm module.
        Convert angle in milli arc-sec to radians.
 - Test of rad2mas() in leastsqnrm module.
        Convert input angle in radians to milli arc sec.
 - Test of sin2deltapistons() in leastsqnrm module.
        Each baseline has one sine and one cosine fringe with a coefficient
        that depends on the piston difference between the two holes that make
        the baseline.  For a 7-hole mask there are 21 baselines and therefore
        there are 42 sine and cosine terms that contribute to the fringe model.
        This function calculates the sine of this piston difference.
 - Test of cos2deltapistons() in leastsqnrm module.
        Each baseline has one sine and one cosine fringe with a coefficient
        that depends on the piston difference between the two holes that make
        the baseline.  For a 7-hole mask there are 21 baselines and therefore
        there are 42 sine and cosine terms that contribute to the fringe model.
        This function calculate the cosine of this piston difference.
 - Test of replacenan() in leastsqnrm module.
        Replace singularities encountered in the analytical hexagon Fourier
        transform with the analytically derived limits. (pi/4)
 - Test of hexpb() in leastsqnrm module.
        Calculate the primary beam for hexagonal holes.
 - Test of model_array in leastsqnrm module.
        Create a model using the specified wavelength.
 - Test of ffc in leastsqnrm module.
        Calculate cosine terms of analytic model.
 - Test of ffs in leastsqnrm module.
        Calculate sine terms of analytic model.
 - Test of return_CAs in leastsqnrm module.
        Calculate the closure amplitudes.
 -  Test of closurephase in leastsqnrm module.
         Calculate closure phases between each pair of holes.
 - Test of redundant_cps in leastsqnrm module.
        Calculate closure phases for each set of 3 holes.
 - Test of populate_symmamparray in leastsqnrm module.
        Populate the symmetric fringe amplitude array.
 - Test of populate_antisymmphasearray in leastsqnrm module.
        Populate the antisymmetric fringe phase array.
 - Test of tan2visibilities in leastsqnrm module.
        From the solution to the fit, calculate the fringe amplitude and phase.
 - Test of multiplyenv in leastsqnrm module.
        Multiply the envelope by each fringe 'image'.

hexee module tests:
+++++++++++++++++++
 -  Test of g_eeAG() in the hexee module.
        Calculate the Fourier transform of one half of a hexagon that is
        bisected from one corner to its diametrically opposite corner.
 -  Test of glimit() in the hexee module.
        Calculate the analytic limit of the Fourier transform of one half of the
        hexagon along eta=0.

analyticnrm2 module tests:
++++++++++++++++++++++++++

 -  Test of PSF() in the analyticnrm2 module 
 -  Test of ASFhex() in the analyticnrm2 module FOR HEX 
 -  Test of interf() in the analyticnrm2 module
 -  Test of phasor() in the analyticnrm2 module

webb_psf module test:
+++++++++++++++++++++

 -  Test of PSF() in the webb_psf module:
        Create a Throughput datamodel, having a dummy filter bandpass data
        that peaks at 1.0 at the center and decreases in the wings.

