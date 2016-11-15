.. _level3-associations:

###################
Level3 Associations
###################

.. toctree::
   :maxdepth: 2

   self
   level3_asn_technical
   level3_asn_rules
   level3_asn_reference

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

User Usage
==========

In general, the only use of associations is to create them from an
association pool using :ref:`asn_generate` and then use them as input
to the appropriate JWST Level3 pipeline step or steps.

When there is need to modify an association, there are two ways of
doing so, depending on the reason for modification. If there is an
issue at the exposure level, either to add or remove an exposure, then
the association pool itself should be modified. After pool
modification, :ref:`asn_generate` is then re-run to create new
associations. This ensures that all associations related to the
exposure are modified.

If there is an issue with a particular association, then one can
simply edit that association as need be. The :ref:`technical
specifications <asn-level3-techspecs>` describe the contents of Level3
associations and how to modify them.

For scripting usage, the high-level routines to use are described in
the :ref:`design` section.
