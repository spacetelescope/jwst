.. _level3-associations:

###################
Level3 Associations
###################

.. toctree::
   :maxdepth: 2

   self
   level3_asn_technical
   level3_asn_rules

.. _level3-asn-jwst-overview:

JWST Associations
=================

The basic chunk in which science data arrives from the observatory is
termed an `exposure`. An exposure contains the data from a single set
of integrations per detector per instrument. In general, it takes many
exposures to make up a single observation, and a whole program is made
up of a large number of observations.

On first arrival, an exposure is termed to be at `Level1b`: The only
transformation that has occured is the extraction of the science data
from the telescope telemetry into a FITS file. At this point, the
science exposures enter the calibration pipeline.

The pipeline consists of two stages: Level2 processing and Level3
processing. Level2 processing is the calibration necessary to remove
instrumental effects from the data. The resulting files contain flux
and spatially calibrated data, called `Level2b` data. The information
is still in individual exposures.

To be truly useful, the exposures need to be combined and, in the case
of multi-object spectrometry, separated, into data that is
source-oriented. This type of calibration is called `Level3`
processing. Due to the nature of the individual instruments, observing
modes, and the interruptability of the observatory itself, how to
group the right exposures together is not straight-forward.

Enter the :ref:`association-generator`. Given a set of exposures,
called the :ref:`Association Pool <asn-pool>`, and a set of rules found in an
:ref:`Association Registry <asn-registry>`, the generator groups the exposures into
individual :ref:`associations <association>`. These associations are
then used as input to the Level3 calibration steps to perform the
transformation from exposure-based data to source-based, high(er)
signal-to-noise data.

In short, Level 3 associations are created running the
:ref:`asn_generate <asn-generate>` task on an :ref:`Association Pool
<asn-pool>` using the default :ref:`Level 3 Association Rules
<level3-asn-rules>` to produce :ref:`level3-associations`.
