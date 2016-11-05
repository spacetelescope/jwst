============
Dark current
============


The dark current step removes dark current from a JWST exposure by subtracting
dark current data stored in a dark reference file.  The reference file records a high signal-to-noise 
ramp of the detector dark signal (i.e., the signal detected in the absence of photons from the sky).  
It is constructed by averaging the individual frames of many long, dark exposures.


.. toctree::
   :maxdepth: 4

   dark_current_reference_files
