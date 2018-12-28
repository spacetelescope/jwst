Tasks in the Package
====================
The coronagraphic package currently consists of the following tasks:

1) :ref:`stack_refs <stack_refs_step>`: Accumulate PSF target images from multiple exposures
   into a single product, for use in other steps
2) :ref:`align_refs <align_refs_step>`: Align PSF images to science target images
3) :ref:`klip <klip_step>`: Apply the KLIP algorithm to fit and subtract an optimized PSF
   from science target images
4) :ref:`hlsp <hlsp_step>`: Create high-level science products from a PSF-subtracted image

These tasks, and others, are applied to an association of coronagraphic exposures using the
pipeline module :ref:`calwebb_coron3 <calwebb_coron3>`.
