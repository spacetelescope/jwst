Description
===========

Processing multiple datasets together allows for the identification of bad pixels
or cosmic-rays that remain in each of the input images, many times at levels which
were not detectable by the :ref:`jump <jump_step>` step. The ``outlier_detection`` step 
implements the following algorithm to identify and flag any remaining cosmic-rays or
other artifacts left over from previous calibrations:

  - build a stack of input data
  
    - all inputs will need to have the same WCS since outlier detection assumes 
      the same flux for each point on the sky, and variations from one image to
      the next would indicate a problem with the detector during readout of that
      pixel
    - if needed, each input will be resampled to a common output WCS
  
  - create a median image from the stack of input data

    - this median operation will ignore any input pixels which have a weight 
      which is too low (<70% max weight)
  
  - create "blotted" data from the median image to exactly match each original
    input dataset
  
  - perform a statistical comparison (pixel-by-pixel) between the median blotted
    data with the original input data to look for pixels with values that are  
    different from the mean value by more than some specified sigma 
    based on the noise model

    - the noise model used relies on the error array computed by previous 
      calibration steps based on the readnoise and calibration errors
  
  - flag the DQ array for the input data for any pixel (or affected neighboring
    pixels) identified as a statistical outlier

The outlier detection step serves as a single interface to apply this general 
process to any JWST data, with specific variations of this algorithm for each 
type of data.  Sub-classes of the outlier detection algorithm have been developed
specifically for 

  - Imaging data
  - IFU spectroscopic data
  - TSO data
  - coronagraphic data
  - spectroscopic data
  
This allows the outlier_detection step to be tuned to the variations in each type 
of JWST data.
           
Reference Files
===============

The ``outlier_detection`` step does not use any reference files.
