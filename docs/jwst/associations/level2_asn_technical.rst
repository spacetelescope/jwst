.. _asn-level2-techspecs:

Level 2 Associations: Technical Specifications
==============================================

Logical Structure
-----------------

Independent of the actual format, all Level 2 associations have the
following structure. Again, the structure is defined and enforced by
the Level 2 schema

  * Top level, or meta, key/values
  * List of products, each consisting of
    
    * Output product name template
    * List of exposure members, each consisting of
      
      * filename of the input exposure
      * Type of exposure
      * Errors from the observatory log
      * Association Candidates this exposure belongs to

.. _asn-level2-example:
   
Example Association
-------------------

The following example will be used to explain the contents of an association::

    {
        "asn_rule": "Asn_Lv2Spec",
        "asn_pool": "jw82600_001_20160304T145416_pool",
        "program": "82600",
        "asn_type": "spec2",
        "products": [
            {
                "name": "test_lrs1",
                "members": [
                    {
                        "expname": "test_lrs1_rate.fits",
                        "exptype": "science"
                    }
                ]
            },
            {
                "name": "test_lrs2bkg",
                "members": [
                    {
                        "expname": "test_lrs2bkg_rate.fits",
                        "exptype": "science"
                    }
                ]
            },
            {
                "name": "test_lrs2",
                "members": [
                    {
                        "expname": "test_lrs2_rate.fits",
                        "exptype": "science"
                    },
                    {
                        "expname": "test_lrs2bkg_rate.fits",
                        "exptype": "background"
                    }
                ]
            }
        ]
    }

Association Meta Keywords
-------------------------

The following are the top-level, or meta, keywords of an association.

program *optional*
  Program number for which this association was created.
  
asn_type *optional*
  The type of association represented. See :ref:`asn-jwst-association-types`

asn_id *optional*
  The association id. The id is what appears in the :ref:`asn-jwst-naming`
  
asn_pool *optional*
  Association pool from which this association was created.

asn_rule *optional*
  Name of the association rule which created this association.
  
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

A list of products that would be produced by this association. For
Level2, each product is an exposure. Each product should have one
``science`` member, the exposure on which the Level2b processing will
occur.

Association products have two components: 

name *optional*
  The string template to be used by Level 2b processing tasks to create
  the output file names. The product name, in general, is a prefix on
  which the individual pipeline and step modules will append whatever
  suffix information is needed.

  If not specified, the Level2b processing modules will create a name
  based off the name of the ``science`` member.

members *required*
  This is a list of the exposures to be used by the Level 2b processing
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
    
    Primary science exposure. For each product, only one exposure can
    be ``science``.
    
  * ``background`` *optional*
    
    Off-target background exposure to subtract.
    
  * ``imprint`` *optional*
    
    Imprint exposure to subtract.
    
  * ``sourcecat`` *optional*
    
    The catalog of sources to extract spectra for. Usually produced by
    ``calwebb_image3`` for wide-field slitless spectroscopy.

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
