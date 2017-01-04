.. _asn-generate:

asn_generate
============

Command Line
------------

.. program-output:: asn_generate --help

API
---

There are two programmatic entry points: the :ref:`Main class
<main>` and the :ref:`generate <generate_function>` function. `Main`
is the highest level entry and is what is instantiated when the
command line `asn_generate` is used. `Main` basically parses the
command line options, creates the :ref:`AssociationPool <asn-pool>`
and :ref:`AssociationRegistry <asn-registry>` instances, calls
`generate`, and saves the resulting associations.

`generate` is the main mid-level entry point. Given an
`AssociationPool` and an `AssociationRegistry`, `generate`
returns a list of associations and the orphaned exposure table.

.. _main:

Main
^^^^

.. autoclass:: jwst.associations.main.Main
   :noindex:
      
.. _generate_function:

generate()
^^^^^^^^^^

.. autofunction:: jwst.associations.generate.generate
   :noindex:
