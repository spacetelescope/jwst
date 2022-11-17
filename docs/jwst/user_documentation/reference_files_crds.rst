.. _reference_files_crds:

=========================================
Reference Files, Parameter Files and CRDS
=========================================

The JWST pipeline uses version-controlled :ref:`reference files <crds_reference_files>` and
:ref:`parameter files <crds_parameter_files>` to supply pipeline steps with necessary data
and set pipeline/step parameters, respectivley. These files both use the ASDF format,
and are managed by the Calibration References Data System (:ref:`CRDS <crds>`) system.

.. warning::

   As of November 10, 2022, the process of deprecating the CRDS PUB Server will start.

   For details, refer to the :ref:`pub-deprecation` page.

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
steps as determined by the JWST instrument teams, based on insturment, observing
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


The CRDS context is usually set by default to always give access
to the most recent reference file deliveries and selection rules - i.e the
'best', most up-to-date set of reference files. On occasion it might be
necessary or desirable to use one of the non-default mappings in order to, for
example, run different versions of the pipeline software or use older versions
of the reference files. This can be accomplished by setting the environment
variable ``CRDS_CONTEXT`` to the desired project mapping version, e.g.

::

  $ export CRDS_CONTEXT='jwst_0421.pmap'

For all information about CRDS, including context lists, see the JWST CRDS
website:

    `https://jwst-crds.stsci.edu/ <https://jwst-crds.stsci.edu/>`_


CRDS Servers
------------
There are two servers available:

  - JWST OPS: https://jwst-crds.stsci.edu
  - JWST PUB: https://jwst-crds-pub.stsci.edu

JWST OPS supports the automatic processing pipeline at STScI. JWST PUB supports
the latest public release of the `jwst` package. Most often, the reference
contexts are one and the same. Regardless, if one wishes to calibrate using the
same exact information as the automatic processing, use JWST OPS. Otherwise, use
of JWST PUB is recommended.

Inside the STScI network, the pipeline defaults the CRDS setup to use JWST OPS with no modifications.
To run the pipeline outside the STScI network or to use a different server, CRDS must be configured by setting
two environment variables:

  - CRDS_PATH: Local folder where CRDS content will be cached.
  - CRDS_SERVER_URL: The server from which to pull reference information

To setup to use JWST OPS, use the following settings:

::

    export CRDS_PATH=$HOME/crds_cache/jwst_ops
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

To setup to use JWST PUB, use the following settings:

::

    export CRDS_PATH=$HOME/crds_cache/jwst_pub
    export CRDS_SERVER_URL=https://jwst-crds-pub.stsci.edu

.. _python_crds_variables:

Setting CRDS Environment Variables in Python
--------------------------------------------

The CRDS environment variables need to be defined *before* importing anything
from `jwst` or `crds`. The examples above show how to set an environment variable in
the shell, but this can also be done within a Python session by using `os.envron`.
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
