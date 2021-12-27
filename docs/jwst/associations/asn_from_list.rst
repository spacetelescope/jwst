.. _asn-from-list:

asn_from_list
=============

Create an association using either the command line tool
``asn_from_list`` or through the Python API using either
:class:`jwst.associations.asn_from_list.Main` or
:func:`jwst.associations.asn_from_list.asn_from_list`

Command Line
------------

.. code-block:: shell

    asn_from_list --help

Usage
-----

Level2 Associations
^^^^^^^^^^^^^^^^^^^

Refer to :ref:`asn-level2-techspecs` for a full description of Level2
associations.

To create a Level2 association, use the following command:

.. code-block:: shell

   asn_from_list -o l2_asn.json -r DMSLevel2bBase *.fits

The ``-o`` option defines the name of the association file to create.

The ``-r DMSLevel2bBase`` option indicates that a Level2 association is
to be created.

Each file in the list will have its own ``product`` in the association file.
When used as input to :ref:`calwebb_image2 <calwebb_image2>` or
:ref:`calwebb_spec2 <calwebb_spec2>`, each product is processed independently,
producing the Level2b result for each product.

For those exposures that require an off-target background or imprint
image, modify the ``members`` list for those exposure, adding a new
member with an ``exptype`` of ``background`` or ``imprint`` as
appropriate. The ``expname`` for these members are the Level2a exposures
that are the background/imprint to use.

An example product that has both a background and imprint exposure
would look like the following::

  "products": [
      {
          "name": "jw99999001001_011001_00001_nirspec",
          "members": [
              {
                  "expname": "jw99999001001_011001_00001_nirspec_rate.fits",
                  "exptype": "science"
              },
              {
                  "expname": "jw99999001001_011001_00002_nirspec_rate.fits",
                  "exptype": "background"
              },
              {
                  "expname": "jw99999001001_011001_00003_nirspec_rate.fits",
                  "exptype": "imprint"
              }
          ]
      }
  ]

Level3 Associations
^^^^^^^^^^^^^^^^^^^

Refer to :ref:`asn-level3-techspecs` for a full description of Level3
associations.

To create a Level3 association, use the following command:

.. code-block:: shell

   asn_from_list -o l3_asn.json --product-name l3_results *.fits

The ``-o`` option defines the name of the association file to create.

The ``--product-name`` will set the ``name`` field that the Level3 calibration
code will use as the output name. For the above example, the output files
created by :ref:`calwebb_image3 <calwebb_image3>`, or other Level3 pipelines,
will all begin with **l3_results**.

The list of files will all become ``science`` members of the
association, with the presumption that all files will be combined.

For coronagraphic or AMI processing, set the ``exptype`` of the
exposures that are the PSF reference exposures to **psf**.  If the
PSF files are not in the ``members`` list, edit the association and add
them as members. An example product with a psf exposure would look
like::

  "products": [
      {
          "name": "jw99999-o001_t14_nircam_f182m-mask210r",
          "members": [
              {
                  "expname": "jw99999001001_011001_00001_nircam_cal.fits",
                  "exptype": "science"
              },
              {
                  "expname": "jw99999001001_011001_00002_nircam_cal.fits",
                  "exptype": "science"
              },
              {
                  "expname": "jw99999001001_011001_00003_nircam_cal.fits",
                  "exptype": "psf"
              }
          ]
      }
  ]


API
---

There are two programmatic entry points: The
:py:class:`~jwst.associations.asn_from_list.Main` is the highest level
entry and is what is instantiated when the command line
``asn_from_list`` is used. ``Main`` handles the command line interface.

:py:func:`~jwst.associations.asn_from_list.asn_from_list` is the main
mid-level entry point.


.. automodapi:: jwst.associations.asn_from_list
