Level 3 Associations
````````````````````
Level 3 associations associate which level2b files to process to
create higher ordered Level 3 product3. Associations are used as the
primary input to the Level 3 processing tasks. The associations
themselves also define the primary output names to be used by the
level 3 processing tasks when creating the output files.

For Level 3, there is a one-to-one correspondence between associations
and the Level 3 products. However, there are cases where one
association can define multiple Level 3 outputs. Also, Level 3 tasks
are free to create as many auxiliary products as seen fit. Regardless,
any outputs created must use the product names, as defined in the
association, as a template for name creation.

Level 3 associations are created running the `asn_generate` task on an
`association pool` using the default `Level 3 Association Rules`.

The structure of an association is defined by a schema. Level 3
associations follow the Level 3 Association schema. Associations are
currently in JSON format.

Logical Structure
=================

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

Example Association
===================

The following example will be used to explain the contents of an association::
  
    {
        "degraded_status": "No known degraded exposures in association.",
        "creation_time": "20160826T131159",
        "asn_type": "image",
        "constraints": "Constraints:\n    opt_elem2: CLEAR\n    pointing_type: SCIENCE\n    detector: (?!NULL).+\n    target_name: 1\n    exp_type: NRC_IMAGE\n    wfsvisit: NULL\n    instrument: NIRCAM\n    opt_elem: F090W\n    program: 99009",
        "asn_pool": "mega_pool",
        "asn_rule": "Asn_Image",
        "targname": "1",
        "program": "99009",
        "products": [
            {
                "name": "jw99009-a3001_t001_nircam_f090w_{product_type}.fits",
                "members": [
                    {
                        "exposerr": null,
                        "expname": "jw_00001_cal.fits",
                        "asn_candidate": "[('o001', 'OBSERVATION')]",
                        "exptype": "SCIENCE"
                    },
                    {
                        "exposerr": null,
                        "expname": "jw_00002_cal.fits",
                        "asn_candidate": "[('o001', 'OBSERVATION')]",
                        "exptype": "SCIENCE"
                    },
                ]
            }
        ]
    }

JSON Format
-----------

As with any JSON file, the basic unit of information is the key/value
pair::
  
  "key": "value"
  
Values can be nested, either by specifying a list, using `[]`
notation, or a dictionary (object in JSON nomenclature), using the `{}` notations::

  "key_with_an_array": [
    value,
    value,
    value
  ],
  "key_with_dictionary": {
    "key": value,
    "key": value
  }

The `json.org <http://www.json.org/>`_ website has a full description
of the format.

Ordering and Indention
^^^^^^^^^^^^^^^^^^^^^^

Neither JSON or the logical layout of the associations themselves
define or require that key/values appear in a particular order. It is
entirely possible that the `products` keyword may appear in the middle
of the file, surround by the other top-level keywords. This is
perfectly acceptable. However, it may be disconcerting at first if one is
editing an association.

What is important is the indention of the nested values. Indention
should be done only with spaces, to ensure that visual inspection is
correct. How much indentation to use is arbitrary, but must be
consistent: All nested information for a key must lie at the same
indentation.

Association Meta Keywords
-------------------------

The following are the top-level, or meta, keywords of an association.

program
  Program number for which this association was created.
  
targname
  Target ID for which this association refers to.

asn_type
  The type of association represented. This is used to indicate what
  Level 3 task should be used to process this association. Possible
  values are:

  * `image`: suitable for CALIMAGE3 processing
  * `spec`: suitable for CALSPECE3 processing
  * `wfs`: Wave front sensing data, used by `wfs_combine`
  * `ami`: Aperture Mask Interferometry
  * `coron`: Coronography
  * `tso`: Time-series Observations

asn_pool
  Association pool from which this association was created.

asn_rule
  Name of the association rule which created this association.
  
degraded_status
  Error status from the observation logs. If none the phrase "No
  known degraded exposures in association." is used.

creation_time
  Time of creation or version id. If time, it will have the format `yyyymmddThhmmss`

constraints
  List of constraints used by the association generator to create this
  association. Format and contents are determined by the defining
  rule.


`products` Keyword
------------------

Association products have to components:

name
  The string template to be used by Level 3 processing tasks to create
  the output file names. The template has one, or more, replacement
  fields to be used by downstream tasks to fill in task-specific
  information. All templates have one replacement field,
  `product_type`. For example, CALIMAGE3 will fill this with the value
  `i2d`.

  Associations of type `spec` will have an extra replacement field,
  `source_id`. This is meant for the number of multi-object modes that
  exist, since target/source information is not known until Level 3
  processing.

members
  This is a list of the exposures to be used by the Level 3 processing
  tasks. This keyword is explained in detail in the next section.

`members` Keyword
-----------------

`members` is a list of objects, each consisting of the following
keywords

expname *required*
  The exposure file name

exptype *required*
  Type of information represented by the exposure. Possible values are

  * `SCIENCE`
  * `TARGET_AQUISITION`

exposerr *optional*
  If there was some issue the occured on the observatory that may have
  affected this exposure, that condition is listed here. Otherwise the
  value is `null`

asn_candidate *optional*
  Contains the list of association candidates this exposure belongs
  to.

Editing the member list
-----------------------

As discussed previously, a member is made up of a number of keywords,
formatted as follows::

  {
      "expname": "jw_00003_cal.fits",
      "exptype": "SCIENCE",
      "exposerr": null,
      "asn_candidate": "[('o001', 'OBSERVATION')]"
  },

To remove a member, simply delete its corresponding set.

To add a member, one need only specify the two required keywords::

  {
      "expname": "jw_00003_cal.fits",
      "exptype": "SCIENCE"
  },
