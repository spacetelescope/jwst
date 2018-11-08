.. _asn-generate:

asn_generate
============

Association generation is done either using the command line tool
``asn_generate`` or through the Python API using either
``Main`` or ``generate``

Command Line
------------

.. code-block:: shell

    asn_generate --help

Association Candidates
^^^^^^^^^^^^^^^^^^^^^^

A full explanation of association candidates be found under the
:ref:`design <design-candidates>` section.

Default Rules
^^^^^^^^^^^^^

The default rules are the :ref:`Level2 <level2-associations>` and
:ref:`Level3 <level3-associations>`. Unless the ``--ignore-default``
option is specified, these rules are included regardless of any other
rules also specified by the ``-r`` options.

DMS Workflow
^^^^^^^^^^^^
The JWST pipeline environment has specific requirements that must be
met by any task running in that environment. The ``--DMS`` option
ensures that ``asn_generate`` conforms to those specifications.

API
---

There are two programmatic entry points: the :ref:`Main class <main>`
and the :ref:`generate <generate_function>` function. ``Main`` is the
highest level entry and is what is instantiated when the command line
``asn_generate`` is used. ``Main`` parses the command line options,
creates the :ref:`AssociationPool <asn-pool>` and
:ref:`AssociationRegistry <asn-registry>` instances, calls ``generate``,
and saves the resulting associations.

``generate`` is the main mid-level entry point. Given an
``AssociationPool`` and an ``AssociationRegistry``, ``generate``
returns a list of associations and the orphaned exposure table.
