Step Arguments
==============
The ``ami_analyze`` step has several optional arguments. In most cases the
default arguments will be suitable but more advanced users may wish to test
other options.

``--oversample`` (int, default=3)
  The oversampling factor to be used in the model fit.

``--rotation`` (float, default=0.0)
  Initial guess for the rotation of the PSF in the input image, in
  units of degrees.

``--psf_offset`` (str, default='0.0 0.0')
  List of PSF offset values to use when creating the model array.

``--rotation_search`` (str, default='-3 3 1')
  List of start, stop, and step values that define the list of
  rotation search values. The default setting of '-3 3 1'
  results in search values of [-3, -2, -1, 0, 1, 2, 3].

``--bandpass`` (str, default=None)
  ASDF file containing suitable array to override filter/source.

``--usebp`` (bool, default=True)
  If `True`, exclude pixels marked DO_NOT_USE from fringe fitting.

``--firstfew`` (int, default=None)
  If not None, process only the first few integrations.

``--chooseholes`` (str, default=None)
  If not None, fit only certain fringes, e.g., ['B4','B5','B6','C2'].

``--affine2d`` (str, default='commissioning')
  ASDF file containing user-defined affine parameters.
  Note that the default argument is a special case;
  'commissioning' is currently the only string other than an ASDF filename
  that is accepted. If `None` is passed, it will perform a rotation search
  (least-squares fit to a PSF model) and use that for the affine transform.

``--run_bpfix`` (bool, default=True)
  Run Fourier bad pixel fix on cropped data.

Creating ASDF files
^^^^^^^^^^^^^^^^^^^
The optional arguments ``bandpass`` and ``affine2d`` must be written to `ASDF <https://asdf-standard.readthedocs.io/>`_
files to be used by the step. The step expects the contents to be stored with particular keys but the format is not currently
enforced by a schema; incorrect ASDF file contents will cause the step to revert back to the defaults for each argument.

Examples of how to create ASDF files containing the properly formatted information for each of the arguments follows.

.. code-block:: python

   # Create a F480M filter + Vega bandpass ASDF file

   import asdf
   from jwst.ami import utils
   from stdatamodels.jwst import datamodels
   from synphot import SourceSpectrum

   # F480M throughput reference file from JWST CRDS
   throughput_file = 'jwst_niriss_throughput_0012.fits'
   nspecbin=19
   throughput_model = datamodels.open(throughput_file)

   filt_spec = utils.get_filt_spec(throughput_model)
   src_spec = SourceSpectrum.from_vega()
   bandpass = utils.combine_src_filt(filt_spec,
                                    src_spec,
                                    trim=0.01,
                                    nlambda=nspecbin)

   # This bandpass has shape (19, 2); each row is [throughput, wavelength]
   asdf_name = 'bandpass_f480m_vega.asdf'
   tree = {"bandpass": bandpass}
   with open(asdf_name, 'wb') as fh:
        af = asdf.AsdfFile(tree)
        af.write_to(fh)
   af.close()
   throughput_model.close()


.. code-block:: python

   # Create an affine transform ASDF file to use for the model

   import asdf
   tree = {
       'mx': 1., # dimensionless x-magnification
       'my': 1., # dimensionless y-magnification
       'sx': 0., # dimensionless x shear
       'sy': 0., # dimensionless y shear
       'xo': 0., # x-offset in pupil space
       'yo': 0., # y-offset in pupil space
       'rotradccw': None
   }

   affineasdf = 'affine.asdf'

   with open(affineasdf, 'wb') as fh:
       af = asdf.AsdfFile(tree)
       af.write_to(fh)
   af.close()
