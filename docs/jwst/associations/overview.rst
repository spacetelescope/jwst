.. _asn-overview:

====================
Association Overview
====================

.. _asn-what-are-associations:

What are Associations?
======================

Associations are basically just lists of things, mostly exposures,
that are somehow related. With respect to JWST and the Data Management
System (DMS), associations have the following characteristics:

- Relationships between multiple exposures are captured in an association.
- An association is a means of identifying a set of exposures that belong together and may be dependent upon one another.
- The association concept permits exposures to be calibrated, archived, retrieved, and reprocessed as a set rather than as individual objects.
- For each association, DMS will generate the most combined and least combined data products.

.. _asn-associations-and-jwst:

Associations and JWST
=====================

The basic chunk in which science data arrives from the observatory is
termed an *exposure*. An exposure contains the data from a single set
of integrations per detector per instrument. In general, it takes many
exposures to make up a single observation, and a whole program is made
up of a large number of observations.

On first arrival, an exposure is termed to be at *Level1b*: The only
transformation that has occurred is the extraction of the science data
from the observatory telemetry into a FITS file. At this point, the
science exposures enter the calibration pipeline.

The pipeline consists of three stages: Stage 1, Stage 2, and Stage 3
processing. Stage 2 processing is the calibration necessary to remove
instrumental effects from the data. The resulting files contain flux
and spatially calibrated data, called *Stage 2b* data. The information
is still in individual exposures.

.. note::

   Older documentation and code may refer to the stages as **levels**. They
   are synonymous.

To be truly useful, the exposures need to be combined and, in the case
of multi-object spectrometry, separated, into data that is
source-oriented. This type of calibration is called *Stage 3*
processing. Due to the nature of the individual instruments, observing
modes, and the interruptibility of the observatory itself, how to
group the right exposures together is not straight-forward.

Enter the :ref:`Association Generator <design-generator>`. Given a set of exposures,
called the :ref:`Association Pool <design-pool>`, and a set of rules found in an
:ref:`Association Registry <design-registry>`, the generator groups the exposures into
individual :ref:`associations <design-association>`. These associations are
then used as input to the Stage 3 calibration steps to perform the
transformation from exposure-based data to source-based, high(er)
signal-to-noise data.

In short, Stage 2 and Stage 3 associations are created running the
:ref:`asn-generate` task on an :py:class:`~jwst.associations.AssociationPool`
using the default :ref:`Level2<asn-level2-techspecs>` and :ref:`Level
3<asn-level3-techspecs>` association rules to produce Stage 2 and Stage 3
associations. When retrieving the data from the archive, users will find the
list of associated data in JSON files that are submitted together with the
requested Stage 2 or Stage 3 data.

Association Pools
-----------------

The information about what data will be associated is constructed with the
information derived from the Astronomer Proposal Tool and the rules on how data
should be associated that are defined by the instrument teams. All the
information from a single proposal is captured in a single file known as the
:ref:`Association Pool<design-pool>`.

.. _asn-usage:

Usage
=====

Users should not need to run the generator. Instead, it is expected that one
edits an already existing association that accompanies the user's JWST data.
Care should be taken if editing an association file.  Keep in mind all input
files listed in the association file are in the same directory as the
association file and no path information can be put in ``expname``, only the
file name.  Or, if need be, an association can be created based on the existing
:ref:`Stage 2 <asn-level2-example>` or :ref:`Stage 3 <asn-level3-example>` examples.
If, however, the user *does* need to run the generator, :ref:`Association Generator
<design-generator>` documentation will be helpful.

Once an association is in-hand, one can pass it as input to a pipeline
routine. For example::

  % strun calwebb_image3  jw12345-o001_20210311t170002_image3_001_asn.json

Programmatically, to read in an Association, one uses the
:py:func:`~jwst.associations.load_asn` function:

.. code-block:: python

   from jwst.associations import load_asn

   with open('jw12345-o001_20210311t170002_image3_001_asn.json') as fp:
       asn = load_asn(fp)

What exactly is returned depends on what the association is. However,
for all Stage 2 and Stage 3 associations, a Python ``dict`` is returned,
whose structure matches that of the JSON or YAML file. Continuing
from the above example, the following shows how to access the first
exposure file name of a Stage 3 associations:

.. code-block:: python

   exposure = asn['products'][0]['members'][0]['expname']

Since most JWST data are some form of a 
:ref:`JWST Data Model<jwst-data-models>`
an association can be opened with
:ref:`datamodels.open<stdatamodels:datamodels-open>`
which returns a
:py:class:`~jwst.datamodels.ModelContainer`. All members of the association that can
be represented as a ``DataModel``, will be available in the ``ModelContainer``
as their respective DataModels.

.. code-block:: python

  from stdatamodels.jwst.datamodels import open as dm_open
  container_model = dm_open('jw12345-o001_20210311t170002_image3_001_asn.json')

Utilities
=========

There are a number of utilities to create user-specific associations that are
documented under :ref:`Association Commands<association-commands>`.
