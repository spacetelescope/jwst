Description
===========

:Class: `jwst.spectral_leak.SpectralLeakStep`
:Alias: spectral_leak
	
The MIRI MRS filters are designed to keep out-of-band light from interfering with the desired first
order wavelengths dispersed in a given band. However, around 12.2 µm (channel 3A) a few-percent spectral leak
admits second-order light from 6 µm (channel 1B) into the bandpass. This results in 
spectra produced by the pipeline containing additional flux around 12.2 µm that is only proportional to the object flux at 6 µm.


Applying the optimal spectral leak correction to MIRI MRS data in the  :ref:`calwebb_spec3 <calwebb_spec3>` pipeline corrects for
this feature in  extracted channel 3A spectrum
for a given target using the channel 1B spectrum of that target (if available). Note that since the channel 1B FOV is smaller
than that for Ch3A no such correction is possible in general for extended sources that fill the entire FOV. An example of an
uncorrected and corrected spectrum is given in the figure below for a G dwarf star.

.. figure:: Example_Spectral_Leak_Corrected.png
   :scale: 50%
   :align: center

Figure: MRS spectral leak as seen in G Dwarf star. The red extracted spectrum does not have the spectral_leak step applied,
while the the black extracted spectrum has the spectral leak correction applied. 

