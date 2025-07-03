.. _reference_files_crds:

=========================================
Reference Files, Parameter Files and CRDS
=========================================

The JWST pipeline uses version-controlled :ref:`reference files <crds_reference_files>` and
:ref:`parameter files <crds_parameter_files>` to supply pipeline steps with necessary data
and set pipeline/step parameters, respectively. These files both use the ASDF format,
and are managed by the Calibration References Data System (:ref:`CRDS <crds>`) system.

.. _crds_reference_files:

Reference Files
================

Most pipeline steps rely on the use of reference files that contain different
types of calibration data or information necessary for processing the data. The
reference files are instrument-specific and are periodically updated as the data
processing evolves and the understanding of the instruments improves. They are
created, tested, and validated by the JWST Instrument Teams. The teams ensure
all the files are in the correct format and have all required header keywords.
The files are then delivered to the Reference Data for Calibration and Tools
(ReDCaT) Management Team. The result of this process is the files being ingested
into the JWST Calibration Reference Data System (CRDS), and made available to
users, the pipeline team and any other ground subsystem that needs access to
them.

Information about all the reference files used by the Calibration Pipeline can
be found at :ref:`reference_file_information`, as well as in the documentation
for each Calibration Step that uses a reference file. Information on reference
file types and their correspondence to calibration steps is described within the
table at :ref:`reference_file_types`.

.. _crds_parameter_files:

Parameter Files
===============

Parameter files, which like reference files are encoded in ASDF and
version-controlled by CRDS, define the 'best' set of parameters for pipeline
steps as determined by the JWST instrument teams, based on instrument, observing
model, filter, etc. They also may evolve over time as understanding of caibration
improves.

By default, when running the pipeline via ``strun`` or using the ``pipeline/step.call()``
method when using the Python interface, the appropriate parameter file will be determined
and retrieved by CRDS to set step parameters.

.. _crds:

CRDS
====

Calibration References Data System (CRDS) is the system that manages the
reference files that the pipeline uses. For the JWST pipeline, CRDS manages both
data reference files as well as parameter reference files which contain step
parameters.

CRDS consists of external servers that hold all available reference files, and
the machinery to map the correct reference files to datasets and download them
to a local cache directory.

When the Pipeline is run, CRDS uses the metadata in the input file to determine
the correct reference files to use for that dataset, and downloads them to a
local cache directory if they haven't already been downloaded so they're
available on your filesystem for the pipeline to use.

**The environment variables `crds_context` and `crds_server` must be set before running the pipeline**

 
.. _crds_context:

Reference Files Mappings (CRDS Context)
---------------------------------------
One of the main functions of CRDS is to associate a dataset with its best
reference files - this mapping is referred to as the 'CRDS context' and is
defined in a `.pmap` file, which itself is version-controlled to allow access to
the reference file mapping at any point in time, and revert to any previous set
of reference files if desired. 


The CRDS context is usually set by default to always give the 'best' reference files
associated with a given pipeline version.
To use a specific CRDS context other than that automatically associated with a given pipeline version
(see https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/crds-migration-to-quarterly-calibration-updates),
the environment variable ``CRDS_CONTEXT`` can be used, e.g.

::

  $ export CRDS_CONTEXT='jwst_1293.pmap'

For all information about CRDS, including context lists, see the JWST CRDS
website:

    `https://jwst-crds.stsci.edu/ <https://jwst-crds.stsci.edu/>`_


CRDS Servers
------------
The CRDS server [1]_ can be found at

::

   https://jwst-crds.stsci.edu

To run the pipeline inside the STScI network, CRDS must be configured to find the CRDS server
by setting the environment variable

::

    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

This server will be used to determine the appropriate CRDS context for a given pipeline
version, and the pipeline will obtain individual reference files within this context from a local shared disk.

To run the pipeline outside the STScI network, CRDS must be configured by setting
two environment variables:

::

    export CRDS_PATH=$HOME/crds_cache/
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

This server will be used to determine the appropriate CRDS context for a given pipeline
version, and the pipeline will automatically download individual
reference files within this context to the local cache specified by ``CRDS_PATH``.

.. [1] Prior to November 10, 2022, there was a second CRDS server available to users,
   at ``https://jwst-crds-pub.stsci.edu``.  This PUB server was set up in anticipation of
   rapid reference file updates during commissioning and Cycle 1, but it was found to
   be unnecessary, due to the smooth transition to science operations.  It was
   decommissioned in March 2023.  The institute retains an internal archive of the
   reference files available on the PUB server at that time.  Observers who wish to
   use historical files from the PUB server will need to file a JWST Pipeline help
   desk ticket to access the information.


CRDS Cache Configuration for Developers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For most pipeline users, the above settings will suffice for establishing a consistent
local cache.  For pipeline developers or testers, however, it is important to be aware
that if you need to switch between CRDS servers (e.g. the `ops` and `test` servers), you
will need to establish a separate cache for each server.  Using the same cache for
more than one server will lead to a corrupted local cache.

For example, the recommended configuration for developers while using the `ops` server is :

::

    export CRDS_PATH=$HOME/crds_cache/jwst_ops
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

and while using the `test` server:

::

    export CRDS_PATH=$HOME/crds_cache/jwst_test
    export CRDS_SERVER_URL=https://jwst-test-crds.stsci.edu

If your cache does become corrupted, the best way to fix it is simply to remove
the local cache and allow subsequent pipeline runs to repopulate it as needed.
For example:

::

    rm -r $CRDS_PATH

For more information on CRDS configuration, see the
`CRDS user guide
<https://jwst-crds.stsci.edu/static/users_guide/environment.html>`__
posted to the JWST CRDS server.

.. _python_crds_variables:

Setting CRDS Environment Variables in Python
--------------------------------------------

The CRDS environment variables need to be defined *before* importing anything
from `jwst` or `crds`. The examples above show how to set an environment variable in
the shell, but this can also be done within a Python session by using `os.environ`.
In general, any scripts should assume the environment variables have been set before the scripts
have run. If one needs to define the CRDS environment variables within a script,
the following code snippet is the suggested method. These lines should be the first
executable lines:

::

   import os
   os.environ['CRDS_PATH'] = 'path_to_local_cache'
   os.environ['CRDS_SERVER_URL'] = 'url-of-server-to-use'

   # Now import anything else needed
   import jwst
