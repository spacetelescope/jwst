*****************************************************
JWST skymatch pipeline step 
*****************************************************

Description
===========

The skymatch step compares the signal levels in the overlap regions
of a set of images and computes the signal offsets for each image
that will minimize the residuals across the entire set. By default
the sky value computed for each image is recorded, but not actually
subtracted from the images.

.. moduleauthor:: Mihai Cara <help@stsci.edu>

.. currentmodule:: jwst.skymatch.skymatch_step

.. automodule:: jwst.skymatch.skymatch_step
   :members:
..   :undoc-members:
