
Description
===========

This routine takes  MIRI  IFU calibrated 2D slope images and produces
3-D spectral cubes. A future version will also create IFU cubes from NIRSPEC data. 
In this cube_build routine the IFU slice distorted and disjointed 2D spectra are corrected
for distortion and put back together into a rectangular cube with three orthogonal axes, two 
spatial and one spectral with regular sampling in the three axes. The disortion information 
should have been incorporated into the calibrated images using the latest assign_wcs pipeline step.

At this time the cube_build method  only creates cubes from single exposures. Currently the information
needed from dithers and pointing information is not available to the cube_build program. Until the
needed telemetry information is available, a future version of the cube building routine 
will allow the user to provide a list offsets to be applied to a list of exposures.


The cube_build package can take as input either: 

  * a single 2D input image updated by assign_wcs

  * an association table (in json format)

In future versions if the user wants to create a cube from a set of exposures, the list of files must be provided
in an association table. There are a number of arguments the user can provide either in a configuration file or
on the command line that control the sampling size of the cube as well as the type of data that is combined to
create a cube. See the **Arguments** section in the documentation for more details.  



Assumption
----------
It is assumed the assign_wcs step has been run on the data, attaching the distortion and pointing
information to the image. It is also assumed that the input data is MIRI  data. 


Background Information
----------------------
The JWST integral field spectrographs obtain simultaneous spectral and spatial data on a relatively compact
region of the sky. The MIRI Medium Resolution Spectrograph (MRS) consists of four integral field units
providing four simultaneous and overlapping fields of view ranging from 3.7" X 3.7" to ~7.7" X 7.7" covering a
wavelength range from 5 to 28 microns. The optics system for the four IFUs is split into two paths. One path
is dedicated to the two short wavelength IFUs and the other one handles the two longer wavelength IFUs.
There is one 1024 X 1024 SCA for each path. Light entering the MRS is spectrally separated into four
channels by dichroic mirrors. Each of these channels has its own IFU that divides the image into several
slices. Each slice is then dispersed using a grating spectrograph and imaged on one half of a SCA. While
four channels are observed simultaneously, each exposure only records the spectral coverage of
approximately one third of the full wavelength range of each channel. The full 5 to 28 micron spectrum is then
obtained by making three exposures using three different gratings and three different dichroic sets. 
We refer to a sub-channel as one of the three possible configurations of the channel where each
sub-channel covers ~1/3 of the full spectrum for the channel. Each of the four channels have a different sampling 
of the field, so the FOV, slice width, number of slices and plate scales are different for each channel. 


Terminology
-----------

**We use the following terminology to define the Spectral range divisions of MIRI:**

- *Channel* the spectral range covered by each MIRI IFU. The channels are labeled as 1, 2, 3 or 4.
- *Sub-Channel* each of the 3 sub-ranges that a channel is divided into. We  will designate these as *Short*, *Medium*, or *Long*.
- *Band*  is one of the 12 sub-ranges the spectral range of the MRS can be divided, each band has unique channel/sub-channel combination, i.e., 
  the shortest wavelength range on MIRI is covered by Band 1-SHORT and the longest is covered by Band 4-LONG.  

**Coordinate Systems:**

The MRS, being an integral field spectrograh, measures the intensity of in a region of the sky as a function of 
wavelength. There are a number of different coordinate systems used in the cube building process. Here is an 
overview of these coordinate systems:

- *Detector System* is defined by the hardware and presents raw detector pixel values. Each detector or SCA 
  will have its own pixel-based coordinate system. In the case of MIRI we have two detector systems because
  the the MIRI IFUs disperse data onto two SCAs.
- *MRS-FOV* this the angular coordinate system attached to the FOV of each MRS band. There are twelve MRS-FOV systems
  for MIRI, since there are twelve bands (1A, 1B, 1C,... 4C). Each system has two orthogonal axes, one parallel 
  (**alpha**) and the other perpendicular (**beta**) to the projection of the long axes of the slices in the FOV. 
- *Telescope (V2,V3)* : the V2,V3 coordinates locate points on  a spherical coordinate system. The frame is tied
  to JWST and applies to the whole field of view, encompassing all the instruments. The coordinate (V2,V3) are Euler
  angles in a spherical frame rather than Cartesian. The transformation between the V2-V3 and MRS-FOV system is fixed 
  mission and is determined during ground testing. 
- *Absolute* is the standard astronomical equatorial system of Ra-Dec. 
- *Cube* is a three dimensional system with two spatial axes and one spectral axis. 


Output
---------

If a single 2D calibrated file is used as the input then the resulting cube for MIRI data will be two cubes (one
for each channel contained in the input file).  In this case the name of the output files are derived from the input
file name plus a suffix  containing the channel number of band of the cube. For example a MIRIFULONG
calibrated image for the shortest wavelength would by default produce two cubes 1) *filename +  CH3_SHORT.fits*
and 2) *filename + CH4_SHORT.fits*. The output cooridinate system for single exposures are in the **alpha** - **beta** coordinate
system. If the cube is created from combining dither data or data from different bands then the output coordinate
system will ultimately  be in Right ascension, declination, but for the foreseeable future these cubes will be
in the V2-V3 system. 

The user can choose to limit the data contained in the cube to a specified list of channels and/or subchannels.
For example, the user can select to only  create a cube by combining data  from  a specific MIRI  channel 
or specific subchannel by using the following arguments:  

* ``--channel #``, where the only valid values for # are 1,2,3,4.
  This argument is only valid for MIRI data. If the ``--channel`` argument is given, then only data corresponding 
  to that channel  will be used in constructing the cube.  
  If the user wants more than one  channel to make cube, then all the values are contained in the string with a space 
  between each channel number. For example, to create a cube with channel 1 and 2 the argument list is 
  ``--channel='1 2'``. If this value is not specified then all the  channels contained in the input list of files  will be used 
  in constructing the cube. 

* ``--subchannel type``, where the only valid values for type are SHORT,MEDIUM,LONG.
  This argument is only valid for MIRI data. If the ``--subchannel`` argument is given, then only data corresponding 
  to that subchannel will be used in  constructing the cube. If the user wants more than one subchannel, then all 
  the values are contained in the string
  with a space between each band type. For example, to create a  cube with band  LONG and MEDIUM the argument 
  list is ``--subchannel='LONG MEDIUM'``. If this value is not specified then all the subchannels
  contained in the input list of files will be used in constructing the cube.

If an association table is used to supply a list of files to create the cubes from, then the output
**basename** is given in the table. The following arguments control the type of cubes to be created:

* ``--CombinedCube`` uses all the files given in the assoication table regardless of channel or band. 
  The output  has the name is *basename_Combined.fits*.
* ``--ChannelCube`` creates a cube for each channel the data covers. The output name is *basename_CH#.fits*, 
  where # is replaced by 1,2,3 or 4.
* ``-- BandCube`` creates a cube for each band the data covers. The output name is *basename_CH#X.fits*. 
  The # represents the channel # (1,2,3 or 4) and the X is replaced by  SHORT, MEDIUM or LONG.
* ``--SingleChannelExposureCube`` creates a cube for each channel in an input calibrated image. 
  The output name is  *calibrated_filename_CH#X.fits*, where # is replaced by 1,2,3, or 4 and X is 
  replaced by SHORT, MEDIUM or LONG.

In order to better explain how all these options can interact lets take the cases where we have
an association table with 12 files: 4 dither positions of MIRIFUSHORT data at each sub-channel (SHORT,MEDIUM,LONG).
* Example 1: ``--ChannelCube`` would create two cubes for each channel containing all the data in the sub-channels.
* Example 2: ``--Channel=1`` would create a single cube for channel 1 using all the data at the three sub-channels;
* Example 3: ``--Channel=1`` ``--BandCube`` would create three cubes, basename_CH1Short.fits, basename_CH1Medium.fits, 
basename_CH1LONG.fits. 



Algorithm
---------
Based on the arguments defining the type of cubes to create, the program loops over each cube type and selects the data from
each exposure that should be included in the cube. The output cube is defined using the WCS information of all the input data.
This output cube WCS defines a field-of-view that encompasses the undistorted footprints on 
the sky of all the input images. The cube sample size in the three dimensions (plate scale) is either set by the user or 
determined from defaults. Each channel has a predefined scale to use for each dimension. If the data consists of more  than one 
channel of data the output scale corresponds to the channel with the smallest scale. For single exposure cubes the output
WCS system will be in **alpha** - **beta**, for dithered expsures or data from different bands the output WCS system will
be in right ascension-declation (or V2-V3 until the needed telemetry information is available to give ra-dec). 

All the pixels on each exposure that are included in output cube are mapped to the cube coordinate system. This input-to-output 
pixel mapping is determined via a mapping function derived from the WCS of each input image  and the WCS of output cube. This 
mapping process corrects for the optical distortions and uses the spacecraft telemetry information in one rebinning step to map 
a pixel from the  the detector to the cube coordinate system. The mapping is actually a series of chained transformations 
(detector -> alpha-beta-lambda), (alpha-beta-lambda -> V2, V3 lambda), (V2-V3-Lambda - > right ascension-declination-lambda),
and (right ascension-declination-lambda -> Cube coordinate1,-Cube Coordinate2-lambda).  The reverse of each transformation 
is also possible. 

The mapping process results in an irregulary spaced "cloud of points" in the cube coordinate system. Each point in the cloud  
contains information of the flux of the original detector pixel and error of this flux. The final flux that is derived for each 
cube pixel (**spaxel**) is a combination of all the "*point cloud*" values with a specified *region of interest* from the center of 
the spaxel. How to best combine the point cloud values into a final flux is an  on-going process. The current method uses a 
weighting function based on the distance between the center of spaxel center and point cloud member as well as the width 
of the PSF and LSF. The width of the MIRI PSF varies with wavelength, broader for longer wavelengths. The resolving power of 
the MRS  varies with wavelength and band.  Adjacent point-cloud elements may in fact originate from 
different exposures rotated from one another and even from different spectral bands. In order to properly weight the 
distances  between the point cloud element and spaxel the distances are determined in the alpha-beta coordinate system and 
then normalized by the width of the PSF and the LSF.

The algorithm for weighting the point cloud values depends on the WCS of the final cube.
If the cube coordinate system is in alpha-beta-lambda then the alpha dimension is contained in naxis 1, the beta dimension
in naxis 2 and wavelenght in naxis 3.  For a V2-V3 cube the V2 dimesnion is in naxis 1, V3 in naxis 3 and wavelength in naxis 3.
For ra-dec cubes, the right ascension is in naxis 1, declination in naxis 2 and wavelength in naxis 3. 
In order to explain how the point cloud values are weigthed we will make the following definations:

* :math:`Radius_x` -is the size of the region of interest in the naxis1 dimension of the cube
* :math:`Radius_y` -is the size of the region of interest in the naxis2 dimension of the cube
* :math:`Radius_z` -is the size of the region of interest in the naxis3 dimension of the cube


* xdistance = distance between point in the cloud and spaxel center in the alpha dimension in units of arc seconds
* ydistance = distance between point in the cloud and spaxel center in the beta dimension in units of arc seconds
* zdistance = distance between point cloud and spaxel center in the lambda dimension in units of microns

These distances are determined in the **alpha** - **beta** system from where the point cloud value orginated. We want to combine
many points -possibly coming from a variety of bands- together. To apply the correct weighting to these points we
normalize the distance between the cube spaxel and point cloud value by the PSF and the LSF which where defined 
in the **alpha**-**beta** coordinate system.  We therefore, transform the cube spaxel coordinates to each **alpha-beta** system
that is found within the region of interset. 


* xnorm  width of the PSF in the alpha dimension in units of arc seconds
* ynorm  width of the PSF in the beta dimension  in units of arc seconds
* znorm width of LSF in lambda dimension in units of microns

* xn = xdistance/xnorm
* yn = ydistance/ynorm
* zn = zdistance/znorm

* n = the number of point cloud points within the region of interest of a given spaxel 

For each spaxel find the n points in the cloud what fall within Radius_X, Radius_Y and Radius_Z. Using these
n points calculated the 

The spaxel flux K =  
:math:`\frac{ \sum_{i=1}^n Flux_i w_i}{\sum_{i=1}^n w_i}`

Where 

:math:`w_i = (\frac{Radius_x – xn}{Radius_x xn})^p + (\frac{Radius_y – yn}{Radius_y yn})^q + (\frac{Radius_z – zn}{Radius_z zn})^r`

The default values for the p,q,r and 2, 2 and 2 respectively. The optiminal choice of these values is still TBD, but 
one should consider the degree of smoothing desired in the interpolation, the density of the point cloud elements,
and the region of interest when chosing these values. 

