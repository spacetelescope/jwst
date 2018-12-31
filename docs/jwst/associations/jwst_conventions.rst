.. _asn-jwst-conventions:

================
JWST Conventions
================

.. _asn-jwst-naming:

Naming Conventions
==================

When produced through the ground processing, all association files are
named according to the following scheme::

  jwPPPPP-TNNNN_YYYYMMDDtHHMMSS_ATYPE_MMM_asn.json

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
  * ``MMM``: A counter for each type of association created.

.. _asn-jwst-association-types:

Association Types
=================

Each association is intended to make a specific science
product. The type of science product is indicated by the ``ATYPE`` field
in the association file name (see :ref:`asn-jwst-naming`), and in the ``asn_type`` meta
keyword of the association itself (see :ref:`asn-association-meta-keywords`).

The pipeline uses this type as the key to indicate which Level 2 or
Level 3 pipeline module to use to process this association.

The current association types are:

  * ``ami3``: Intended for ``calwebb_ami3.cfg`` processing
  * ``coron3``: Intended for ``calwebb_coron3.cfg`` processing
  * ``image2``: Intended for ``calwebb_image2.cfg`` processing
  * ``image3``: Intended for ``calwebb_image3.cfg`` processing
  * ``nrslamp-spec2``: Intended for ``calwebb_nrslamp-spec2.cfg`` processing
  * ``spec2``: Intended for ``calwebb_spec2.cfg`` processing
  * ``spec3``: Intended for ``calwebb_spec3.cfg`` processing
  * ``tso3``: Intended for ``calwebb_tso3.cfg`` processing
  * ``tso-image2``: Intended for ``calwebb_tso-image2.cfg`` processing
  * ``tso-spec2``: Intended for ``calwebb_tso-spec2.cfg`` processing
  * ``wfs-image2``: Intended for ``calwebb_wfs-image2.cfg`` processing
  * ``wfs-image3``: Intended for ``calwebb_wfs-image3.cfg`` processing
