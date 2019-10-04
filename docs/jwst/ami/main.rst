Steps in the Package
====================
The Aperture Masking Interferometry (AMI) package currently consists
of three steps:

1) :ref:`ami_analyze <ami_analyze_step>`: apply the LG algorithm to a NIRISS AMI exposure
2) :ref:`ami_average <ami_average_step>`: average the results of LG processing for multiple
   exposures
3) :ref:`ami_normalize <ami_normalize_step>`: normalize the LG results for a science target
   using LG results from a reference PSF target

These three steps are applied to an association of AMI exposures using the
pipeline module :ref:`calwebb_ami3 <calwebb_ami3>`.

