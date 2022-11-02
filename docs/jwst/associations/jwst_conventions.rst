.. _asn-jwst-conventions:

================
JWST Conventions
================

.. _asn-jwst-naming:

Association Naming
==================

When produced through the ground processing, all association files are
named according to the following scheme::

  jwPPPPP-TNNNN_YYYYMMDDtHHMMSS_ATYPE_MMMMM_asn.json

where:

  * ``jw``: All JWST-related products begin with ``jw``
  * ``PPPPP``: 5 digit proposal number
  * ``TNNNN``: Candidate Identifier. Can be one of the following:

    * ``oNNN``: Observation candidate specified by the letter ``o`` followed
      by a 3 digit number.
    * ``c1NNN``: Association candidate, specified by the letter 'c',
      followed by a
      number starting at 1001.
    * ``a3NNN``: Discovered whole program associations, specified by the
      letter 'a', followed by a number starting at 3001
    * ``rNNNN``: Reserved for future use. If you see this in practice,
      file an issue to have this document updated.

  * ``YYYYMMDDtHHMMSS``: This is generically referred to as the ``version_id``.
    DMS specifies this as a  timestamp. Note:
    When used outside the workflow, this field is user-specifiable.
  * ``ATYPE``: The type of association. See
    :ref:`asn-jwst-association-types`
  * ``MMMMM``: A counter for each type of association created.

.. _asn-jwst-association-types:

Association Types
=================

Each association is intended to make a specific science
product. The type of science product is indicated by the ``ATYPE`` field
in the association file name (see :ref:`asn-jwst-naming`), and in the ``asn_type`` meta
keyword of the association itself (see :ref:`asn-level3-meta-keywords`).

The pipeline uses this type as the key to indicate which Level 2 or
Level 3 pipeline module to use to process this association.

The current association types are:

  * ``ami3``: Intended for :ref:`calwebb_ami3 <calwebb_ami3>` processing
  * ``coron3``: Intended for :ref:`calwebb_coron3 <calwebb_coron3>` processing
  * ``image2``: Intended for :ref:`calwebb_image2 <calwebb_image2>` processing
  * ``image3``: Intended for :ref:`calwebb_image3 <calwebb_image3>` processing
  * ``nrslamp-spec2``: Intended for :ref:`calwebb_spec2 <calwebb_spec2>` processing
  * ``spec2``: Intended for :ref:`calwebb_spec2 <calwebb_spec2>` processing
  * ``spec3``: Intended for :ref:`calwebb_spec3 <calwebb_spec3>` processing
  * ``tso3``: Intended for :ref:`calwebb_tso3 <calwebb_tso3>` processing
  * ``tso-image2``: Intended for :ref:`calwebb_image2 <calwebb_image2>` processing
  * ``tso-spec2``: Intended for :ref:`calwebb_spec2 <calwebb_spec2>` processing
  * ``wfs-image2``: Intended for :ref:`calwebb_image2 <calwebb_image2>` processing
  * ``wfs-image3``: Intended for :ref:`calwebb_wfs-image3 <calwebb_wfs-image3>` processing

Field Guide to File Names
=========================

The high-level distinctions between stage 2, stage 3, exposure-centric, and
target-centric files can be determined by the following file patterns. These
patterns are not intended to fully define all the specific types of files there
are. However, these are the main classifications, from which the documentation
for the individual calibrations steps and pipelines will describe any further
details.

The most general regex matches all files that have been produced by Stage 3 processing::

    .+[aocr][0-9]{3:4}.+

The following regexes differentiate between exposure-centric and target-centric files.

- Files containing exposure-centric data

  The following regex matches files names produced by either Stage 2 or 3
  calibration and containing exposure-centric data::

    jw[0-9]{11}_[0-9]{5}_[0-9]{5}_.+\.fits

- Files containing target-centric data

  The following regex matches file names produced by Stage 3 calibration and containing target-centric data::

    jw[0-9]{5}-[aocr][0-9]{3:4}_.+
