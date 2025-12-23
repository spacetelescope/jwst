More Information for JWST Users
===============================

The documentation hosted here focuses primarily on calibration software
implementation and use.

For more user-focused information on JWST instruments and calibration methods, please
see the
`JWST science calibration pipeline <https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline>`__
pages hosted on `JDox <https://jwst-docs.stsci.edu>`__.

More information on the
`latest build <https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/jwst-operations-pipeline-build-information>`__
in use for data processing operations is also available on JDox.

.. _jwst-version-scheme:

Versioning Scheme
-----------------

.. note::

   This section addresses the release version to PyPI.
   It is not to be confused with the DMS quarterly build version,
   which has a different numbering scheme (see
   `software vs. DMS builds <https://github.com/spacetelescope/jwst?tab=readme-ov-file#software-vs-dms-build-version-map>`_).

The ``jwst`` package does not follow strict `semantic versioning <https://semver.org/>`_,
in that a "minor" version bump may also contain non-backward compatible changes.
However, a "patch" release will strive to contain only bug fixes.

.. _jwst-public-vs-private-api:

API: Public vs Private
----------------------

As per Python convention, any API name that starts with underscore
(e.g., ``_my_private_function``) is considered private.
Any API not officially documented (i.e., you only found it after some extensive
code-diving) is also considered private.

However, contrary to common practice, the ``jwst`` package also considers what would
normally be public API (i.e., no leading underscore and fully documented)
to be private unless stated as follows. **Any API not listed below, even if
fully documented, could be removed or changed without notice.**

These are the only software interfaces that ``jwst`` considers public API and
thus follow :ref:`jwst-deprecation-policy`:

* Calibration step classes (e.g., Spec2Pipeline, JumpStep, AmiNormalizeStep)
* Command-line interface for calibration steps (e.g., ``strun calwebb_image3``)
* Data models (e.g., removal of a WCS transform converter causing I/O incompatibility)

.. _jwst-deprecation-policy:

API Deprecation Policy
----------------------

.. note::

   This section only covers API changes that are not backward-compatible.
   Adding a new API is usually backward compatible and thus not covered here.

Occasionally, there is a need to modify or remove a given API
(e.g., introduction of new algorithm or instrumentation changes).
For public API (see :ref:`jwst-public-vs-private-api`), when such a
breaking change happens, if possible, ``jwst`` would first deprecate
the item being removed or renamed; this could be a input keyword,
step configuration, function, class, etc.

A deprecated item would emit a ``DeprecationWarning`` via the Python warnings
system (and also to the log, where appropriate). Ideally, the deprecation
warning would provide at least the following information:

* which ``jwst`` version first deprecated it, and
* alternative or replacement, if available.

Deprecation usually cannot happen in a "patch" release (see :ref:`jwst-version-scheme`)
unless it somehow addresses critical security vulnerability or mission needs.
All deprecations would be stated clearly in the API documentation and
the release change log.

The deprecation period starts the moment it is first released.
Ideally, it would last at least two "minor" releases, or longer if the item
is high-impact. After the period elapses, then this deprecated item would
be removed completely in the following "minor" release.

Impacted users (e.g., scientists or notebook maintainers) should switch
away from the deprecated API early in its deprecation period, to give
yourselves ample time to address downstream changes. If you think you
are unable to perform this switch within the deprecation period, you could
either pin to a maximum version of ``jwst`` before the API removal, or
contact ``jwst`` developers to request an extension of the deprecation period.
Any extension request must have sufficient justification on why it is worth
the extra maintenance burden to do so.
