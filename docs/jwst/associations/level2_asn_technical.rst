.. _asn-level2-techspecs:

Stage 2 Associations: Technical Specifications
==============================================

Logical Structure
-----------------

All stage 2 associations have the following structure. The structure is defined
and enforced by the stage 2 schema.

  * :ref:`Informational Meta Keywords<asn-level2-meta-keywords>`
  * List of :ref:`products<asn-level2-products>`, each consisting of
    
    * Output product name
    * List of :ref:`exposure members<asn-level2-members>`, each consisting of
      
      * Filename of the exposure that is a member of this association
      * Type of exposure
      * If present, information about errors from the observatory log

.. _asn-level2-example:
   
Example Association
-------------------

The following example will be used to explain the contents of an association::

    {
        "asn_type": "image2",
        "asn_rule": "candidate_Asn_Lv2Image",
        "version_id": "20210610t121508",
        "code_version": "1.2.3",
        "degraded_status": "No known degraded exposures in association.",
        "program": "00623",
        "constraints": "DMSAttrConstraint({'name': 'program', 'sources': ['program'], 'value': '623'})\nDMSAttrConstraint({'name': 'is_tso', 'sources': ['tsovisit'], 'value': None})\nDMSAttrConstraint({'name': 'instrument', 'sources': ['instrume'], 'value': 'miri'})\nDMSAttrConstraint({'name': 'detector', 'sources': ['detector'], 'value': 'mirimage'})\nDMSAttrConstraint({'name': 'opt_elem', 'sources': ['filter'], 'value': 'f1130w'})\nDMSAttrConstraint({'name': 'opt_elem2', 'sources': ['pupil', 'grating'], 'value': None})\nDMSAttrConstraint({'name': 'opt_elem3', 'sources': ['fxd_slit'], 'value': None})\nDMSAttrConstraint({'name': 'subarray', 'sources': ['subarray'], 'value': 'brightsky'})\nDMSAttrConstraint({'name': 'channel', 'sources': ['channel'], 'value': None})\nDMSAttrConstraint({'name': 'slit', 'sources': ['fxd_slit'], 'value': None})\nConstraint_Image_Science({'name': 'exp_type', 'sources': ['exp_type'], 'value': 'mir_image'})\nConstraint_Single_Science({'name': 'single_science', 'value': False})\nDMSAttrConstraint({'name': 'asn_candidate', 'sources': ['asn_candidate'], 'value': \"\\\\('o037',\\\\ 'observation'\\\\)\"})",
        "asn_id": "o037",
        "asn_pool": "jw00623_20210610t121508_pool",
        "target": "9",
        "products": [
            {
                "name": "jw00623037001_02101_00001_mirimage",
                "members": [
                    {
                        "expname": "jw00623037001_02101_00001_mirimage_rate.fits",
                        "exptype": "science",
                        "exposerr": "null"
                    }
                ]
            }
        ]
    }

.. _asn-level2-meta-keywords:

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
  Target ID for which this association refers to. JWST currently uses
  the TARGETID header keyword in the stage 2 exposure files, but there
  is no formal restrictions on value.
  
version_id *optional*
  Version identifier. DMS uses a time stamp with the format
  ``yyyymmddthhmmss``
  Can be None or NULL

.. _asn-level2-products:

``products`` Keyword
^^^^^^^^^^^^^^^^^^^^

A list of products that would be produced by this association. For
stage 2, each product is an exposure. Each product should have one
``science`` member, the exposure on which the stage 2 processing will
occur.

Association products have two components: 

name *optional*
  The string template to be used by stage 2 processing tasks to create
  the output file names. The product name, in general, is a prefix on
  which the individual pipeline and step modules will append whatever
  suffix information is needed.

  If not specified, the stage 2 processing modules will create a name
  based off the name of the ``science`` member.

members *required*
  This is a list of the exposures to be used by the stage 2 processing
  tasks. This keyword is explained in detail in the next section.

.. _asn-level2-members:

``members`` Keyword
^^^^^^^^^^^^^^^^^^^

``members`` is a list of dictionaries, one for each member exposure in the
association. Each member has the following keywords.

expname *required*
  The exposure file name.  This must be a filename only, with no path.  This
  file must be in the same directory as the association file, so path is not
  needed.  If path data is included, an exception is raised when loading the
  association file.

exptype *required*
  Type of information represented by the exposure. Possible
  values are as follows. *Note that this is not the same as the exposure
  ``EXP_TYPE`` header keyword.*

  * ``science``: Primary science exposure. For each product, only one exposure can
    be ``science``.
    
  * ``background``: Background exposure to subtract.
    
  * ``imprint``: Imprint exposure to subtract.
    
  * ``sourcecat``: The catalog of sources to extract spectra for. Usually produced by
    :ref:`calwebb_image3 <calwebb_image3>` for wide-field slitless spectroscopy.

  * ``segmap``: The 2D segmentation map used to produce the source catalog. Usually produced by
    :ref:`calwebb_image3 <calwebb_image3>` for wide-field slitless spectroscopy.

  * ``direct_image``: The direct image used to produce the source catalog.

Editing the member list
^^^^^^^^^^^^^^^^^^^^^^^

As discussed previously, a member is made up of a number of keywords,
formatted as follows::

  {
      "expname": "jw_00003_cal.fits",
      "exptype": "science",
  },

To remove a member, simply delete its corresponding set.

To add a member, one need only specify the two required keywords::

  {
      "expname": "jw_00003_cal.fits",
      "exptype": "science"
  },
