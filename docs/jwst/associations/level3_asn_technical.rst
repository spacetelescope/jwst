.. _asn-level3-techspecs:

Level 3 Associations: Technical Specifications
==============================================

Logical Structure
-----------------

Independent of the actual format, all Level 3 associations have the
following structure. Again, the structure is defined and enforced by
the Level 3 schema

  * Top level, or meta, key/values
  * List of products, each consisting of

    * Output product name template
    * List of exposure members, each consisting of

      * filename of the input exposure
      * Type of exposure
      * Errors from the observatory log
      * Association Candidates this exposure belongs to

.. _asn-level3-example:

Example Association
-------------------

The following example will be used to explain the contents of an association::

    {
        "degraded_status": "No known degraded exposures in association.",
        "version_id": "20160826t131159",
        "asn_type": "image3",
        "asn_id": "c3001",
        "constraints": "Constraints:\n    opt_elem2: CLEAR\n    detector: (?!NULL).+\n    target_name: 1\n    exp_type: NRC_IMAGE\n    wfsvisit: NULL\n    instrument: NIRCAM\n    opt_elem: F090W\n    program: 99009",
        "asn_pool": "mega_pool",
        "asn_rule": "Asn_Image",
        "target": "1",
        "program": "99009",
        "products": [
            {
                "name": "jw99009-a3001_t001_nircam_f090w",
                "members": [
                    {
                        "exposerr": null,
                        "expname": "jw_00001_cal.fits",
                        "asn_candidate": "[('o001', 'observation')]",
                        "exptype": "science"
                    },
                    {
                        "exposerr": null,
                        "expname": "jw_00002_cal.fits",
                        "asn_candidate": "[('o001', 'observation')]",
                        "exptype": "science"
                    }
                ]
            }
        ]
    }

.. _asn-association-meta-keywords:

Association Meta Keywords
-------------------------

The following are the top-level, or meta, keywords of an association.

program *optional*
  Program number for which this association was created.

target *optional*
  Target ID for which this association refers to. DMS currently uses
  the TARGETID header keyword in the Level2 exposure files, but there
  is no formal restrictions on value.

asn_type *optional*
  The type of association represented. See :ref:`asn-jwst-association-types`

asn_id *optional*
  The association id. The id is what appears in the :ref:`asn-jwst-naming`

asn_pool *optional*
  Association pool from which this association was created.

asn_rule *optional*
  Name of the association rule which created this association.

degraded_status *optional*
  Error status from the observation logs. If none the phrase "No
  known degraded exposures in association." is used.

version_id *optional*
  Version identifier. DMS uses a time stamp with the format
  ``yyyymmddthhmmss``
  Can be None or NULL

constraints *optional*
  List of constraints used by the association generator to create this
  association. Format and contents are determined by the defining
  rule.


``products`` Keyword
^^^^^^^^^^^^^^^^^^^^

Association products have two components:

name *optional*
  The string template to be used by Level 3 processing tasks to create
  the output file names. The product name, in general, is a prefix on
  which the individual pipeline and step modules will append whatever
  suffix information is needed.

  If not specified, the Level3 processing modules will create a name root.

members *required*
  This is a list of the exposures to be used by the Level 3 processing
  tasks. This keyword is explained in detail in the next section.

``members`` Keyword
^^^^^^^^^^^^^^^^^^^

``members`` is a list of objects, each consisting of the following
keywords

expname *required*
  The exposure file name

exptype *required*
  Type of information represented by the exposure. Possible values are

  * ``science`` *required*

    The primary science exposures. There is usually more than one
    since Level3 calibration involves combining multiple science
    exposures. However, at least one exposure in an association needs
    to be ``science``.
    
  * ``psf`` *optional*

    Exposures that should be considered PSF references for
    coronagraphic and AMI calibration.

exposerr *optional*
  If there was some issue the occured on the observatory that may have
  affected this exposure, that condition is listed here. Otherwise the
  value is ``null``

asn_candidate *optional*
  Contains the list of association candidates this exposure belongs
  to.

Editing the member list
^^^^^^^^^^^^^^^^^^^^^^^

As discussed previously, a member is made up of a number of keywords,
formatted as follows::

  {
      "expname": "jw_00003_cal.fits",
      "exptype": "science",
      "exposerr": null,
      "asn_candidate": "[('o001', 'observation')]"
  },

To remove a member, simply delete its corresponding set.

To add a member, one need only specify the two required keywords::

  {
      "expname": "jw_00003_cal.fits",
      "exptype": "science"
  },
