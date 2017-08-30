Reference File
==============

The gain_scale correction step uses the gain reference file. The only purpose
of the reference file is to retrieve the `GAINFACT` keyword value from its
header (the reference file data are not used in any way). If the `ramp_fit`
step, which also uses the gain reference file, succeeded in finding the
`GAINFACT` keyword in this reference file, it will store the value in the
`GAINFACT` keyword in the science data, in which case the `gain_scale` step
will not reload the gain reference file.
