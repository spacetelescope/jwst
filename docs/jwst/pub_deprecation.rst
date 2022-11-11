.. _pub-deprecation:

CRDS PUB Server Freeze and Deprecation
======================================

Why and When
------------

As of November 10, 2022, all observers should use the standard CRDS OPS server
for JWST calibration reference files:

`https://jwst-crds.stsci.edu <https://jwst-crds.stsci.edu>`_

The PUB server:

`https://jwst-crds-pub.stsci.edu <https://jwst-crds-pub.stsci.edu>`_

was set up in anticipation of rapid reference file updates during commissioning
and Cycle 1. However, due to the trouble-free commissioning process,
the smooth transition to science operations, and the subsequent
confusion that has resulted from having two servers, it has been
decided that the PUB server is no longer needed and will be
decommissioned. To make this transition as smooth as possible, this
update will take place in stages.

On November 10, 2022, all observers should begin to transition to using only the CRDS
OPS server, `https://jwst-crds.stsci.edu <https://jwst-crds.stsci.edu>`_. See
the :ref:`software documentation <pub-deprecation-transition>` for
instructions about how to configure CRDS.

In mid-November, the PUB server will be frozen and a notification will be sent via the
reference file mailing lists and the JWST Observer News mailing list.

On December 2nd, access to the PUB server will no longer be available externally.
The frozen PUB database will be maintained internally for 3 months. On March 1, the PUB server will
be fully decommissioned and the institute will
retain an internal archive of the maps and calibration reference files.
Observers who wish to use historical files from the PUB server in the future will need to file a
JWST Pipeline help desk ticket to access the information.

Part of the decommissioning process will include establishing guidance for how best to
maintain reproducibility for new papers and for already-published papers that used the PUB server. This information
will be included in the new JDox page:

`How to Cite Your Data Reduction and Reference Files <https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline-overview/jwst-data-calibration-reference-files>`_


.. _pub-deprecation-transition:

Transition Procedure
--------------------

If using the PUB server, there are two simple tasks that need to be done to
ensure a successful transition from using the PUB server to the JWST OPS server.

First, the folder containing the local CRDS cache, pointed to by the environment
variable CRDS_PATH, should be cleared. Under Linux-like systems, the command is

::

   $ rm -rf $CRDS_PATH/*

Second, ensure that the environment variable CRDS_SERVER_URL is pointing to the
JWST OPS server, https://jwst-crds.stsci.edu:

::

   $ export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

Following these two steps ensures that further calibration processing will use
references from the standard CRDS server.
