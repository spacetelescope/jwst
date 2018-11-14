.. _asn-generate:

asn_generate
============

Association generation is done either using the command line tool
``asn_generate`` or through the Python API using either
:py:class:`~jwst.associations.Main` or :py:func:`~jwst.associations.generate`.

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

The default rules are the :ref:`Level2 <asn-level2-techspecs>` and
:ref:`Level3 <asn-level3-techspecs>`. Unless the ``--ignore-default``
option is specified, these rules are included regardless of any other
rules also specified by the ``-r`` options.

DMS Workflow
^^^^^^^^^^^^
The JWST pipeline environment has specific requirements that must be
met by any task running in that environment. The ``--DMS`` option
ensures that ``asn_generate`` conforms to those specifications.

API
---

There are two programmatic entry points: the :py:class:`~jwst.associations.Main`
class and the :py:func:`~jwst.associations.generate` function.
:py:class:`~jwst.associations.Main` is the highest level entry and is what is
instantiated when the command line ``asn_generate`` is used.
:py:class:`~jwst.associations.Main` parses the command line options, creates the
:py:class:`~jwst.associations.AssociationPool` and :py:class:`~jwst.associations.AssociationRegistry`
instances, calls ``generate``, and saves the resulting associations.

:py:func:`~jwst.associations.generate` is the main mid-level entry point. Given
an :py:class:`~jwst.associations.AssociationPool` and an
:py:class:`~jwst.associations.AssociationRegistry`,
:py:func:`~jwst.associations.generate` returns a list of associations.
