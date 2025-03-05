.. _emicorr_reffile:

EMICORR Reference File
----------------------

:REFTYPE: EMICORR
:Data model: `~jwst.datamodels.EmiModel`

The EMICORR reference file contains data necessary for removing
noise from contaminating EMI frequencies in MIRI ramp data.

.. include:: ../references_general/emi_selection.inc

.. include:: ../includes/standard_keywords.inc

EMICORR Reference File Format
+++++++++++++++++++++++++++++

MIRI EMICORR reference files are in ASDF format.  The EMICORR
reference file contains the frequencies for which the image will
be corrected.

Example file contents::

  root (AsdfObject)
    ├─frequencies (dict)
      │ ├─Hz10 (dict)
      │ │ ├─frequency (float): 10.039216
      │ │ └─phase_amplitudes (NDArrayType): shape=(500,), dtype=float32
      │ ├─Hz10_slow_MIRIFULONG (dict)
      │ │ ├─frequency (float): 10.039216
      │ │ └─phase_amplitudes (NDArrayType): shape=(250,), dtype=float32
    └─subarray_cases (dict)
      ├─BRIGHTSKY (dict)
      │ ├─frameclocks (int): 23968
      │ ├─freqs (dict)
      │ │ ├─FAST (list)
      │ │ │ └─[0] (str): Hz10
      │ │ └─SLOW (dict)
      │ │   ├─MIRIFULONG (list)
      │ │   │ └─[0] (str): Hz10_slow_MIRIFULONG
      │ │   ├─MIRIFUSHORT (list)
      │ │   │ └─[0] (str): Hz10_slow_MIRIFUSHORT
      │ │   └─MIRIMAGE (list)
      │ │     └─[0] (str): Hz10_slow_MIRIMAGE
      │ └─rowclocks (int): 82
      ├─MASK1140 (dict)
      │ ├─frameclocks (int): 23968
      │ ├─freqs (dict)
      │ │ ├─FAST (list)
      │ │ │ ├─[0] (str): Hz390
      │ │ │ └─[1] (str): Hz10
      │ │ └─SLOW (dict)
      │ │   ├─MIRIFULONG (list)
      │ │   │ ├─[0] (str): Hz390
      │ │   │ └─[1] (str): Hz10_slow_MIRIFULONG
      │ │   ├─MIRIFUSHORT (list)
      │ │   │ ├─[0] (str): Hz390
      │ │   │ └─[1] (str): Hz10_slow_MIRIFUSHORT
      │ │   └─MIRIMAGE (list)
      │ │     ├─[0] (str): Hz390
      │ └─rowclocks (int): 82


Frequency Selection
^^^^^^^^^^^^^^^^^^^

The frequency to be corrected for will be selected according to the
dictionary contained in the key ``subarray_cases`` in the reference
file. This contains the subarray name, readout pattern, and detector
to match to the names of the frequencies to be used in the correction.
