.. _asn-level3-techspecs:

Stage 3 Associations: Technical Specifications
==============================================

Logical Structure
-----------------

Independent of the actual format, all stage 3 associations have the
following structure. Again, the structure is defined and enforced by
the stage 3 schema

  * :ref:`Informational Meta Keywords<asn-level3-meta-keywords>`
  * List of :ref:`products<asn-level3-products>`, each consisting of
    
    * Output product name
    * List of :ref:`exposure members<asn-level3-members>`, each consisting of
      
      * Filename of the exposure that is a member of this association
      * Type of exposure
      * If present, information about errors from the observatory log
      * If present, the Association Candidates this exposure belongs to

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

.. _asn-level3-meta-keywords:

Association Meta Keywords
-------------------------

The following are the informational, meta keywords of an association.

asn_id *required*
  The association id. The id is what appears in the :ref:`asn-jwst-naming`

asn_pool *required*
  Association pool from which this association was created.

asn_rule *optional*
  Name of the association rule which created this association.

asn_type *optional*
  The type of association represented. See :ref:`asn-jwst-association-types`

code_version *optional*
  The version of the generator that created this association. Typically this is the version
  of the jwst python package.

constraints *optional*
  List of constraints used by the association generator to create this
  association. Format and contents are determined by the defining
  rule.

degraded_status *optional*
  If any of the included members have an actual issue,
  as reported by the ``exposerr`` keyword, ``degraded_status`` will have the
  value ``One or more members have an error associated with them.`` If no errors
  are indicated, the value will be ``No known degraded exposures in
  association.``

program *optional*
  Program number for which this association was created.

target *optional*
  Target ID to which this association refers. JWST currently uses
  the TARGETID header keyword in the stage 2 exposure files, but there
  are no formal restrictions on value.
  
version_id *optional*
  Version identifier. DMS uses a time stamp with the format
  ``yyyymmddthhmmss``
  Can be None or NULL

.. _asn-level3-products:

``products`` Keyword
^^^^^^^^^^^^^^^^^^^^

Association products have two components:

name *optional*
  The string template to be used by stage 3 processing tasks to create
  the output file names. The product name, in general, is a prefix on
  which the individual pipeline and step modules will append whatever
  suffix information is needed.

  If not specified, the stage 3 processing modules will create a name root.

members *required*
  This is a list of the exposures to be used by the stage 3 processing
  tasks. This keyword is explained in detail in the next section.

.. _asn-level3-members:

``members`` Keyword
^^^^^^^^^^^^^^^^^^^

``members`` is a list of dictionaries, one for each member exposure in the
association. Each member has the following keywords.

expname *required*
  The exposure file name

exptype *required*
  Type of information represented by the exposure. Possible values are

  * ``science``: The primary science exposures. There is usually more than one
    since stage 3 calibration involves combining multiple science
    exposures. However, at least one exposure in an association needs
    to be ``science``.

  * ``background``: Exposures used for background subtraction.
    
  * ``psf``: Exposures that should be considered PSF references for
    coronagraphic and AMI calibration.

exposerr *optional*
  If there was some issue the occurred on the observatory that may have
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
