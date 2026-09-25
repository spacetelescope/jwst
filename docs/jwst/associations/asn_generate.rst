.. _asn-generate:

asn_generate
============

Association generation is done either using the command line tool
``asn_generate`` or through the Python API using either:

* :class:`~jwst.associations.main.Main`, which is the highest level
  entry and is what is instantiated when the command line is used.
  It parses the command line options, creates the
  :class:`~jwst.associations.pool.AssociationPool` and
  :class:`~jwst.associations.registry.AssociationRegistry` instances,
  calls :func:`~jwst.associations.generator.generate.generate`,
  and saves the resulting associations.
* :func:`~jwst.associations.generator.generate.generate`,
  which is the main mid-level entry point. Given
  an :class:`~jwst.associations.pool.AssociationPool` and an
  :class:`~jwst.associations.registry.AssociationRegistry`,
  it returns a list of associations.

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

The default rules are the :ref:`Stage 2 <asn-level2-techspecs>` and
:ref:`Stage 3 <asn-level3-techspecs>`. Unless the ``--ignore-default``
option is specified, these rules are included regardless of any other
rules also specified by the ``-r`` options.

DMS Workflow
^^^^^^^^^^^^
The JWST pipeline environment has specific requirements that must be
met by any task running in that environment. The ``--DMS`` option
ensures that ``asn_generate`` conforms to those specifications.
