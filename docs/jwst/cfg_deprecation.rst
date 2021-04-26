.. _cfg_usage_deprecation_notice:

CFG Usage Deprecation Notice
============================

As of March 18, 2021, a significant change to how JWST pipelines operate was
completed and pushed to the JWST master branch on github. Theoretically the
change should be transparent. However, we are all familiar with the difference
between theory and practice and hence we want to alert all users.

Originally, how the pipelines operated was determined by a set of configuration
(CFG) files that were delivered as part of the JWST package. These configuration
files were retrieved using the ``collect_pipeline_cfgs`` command. The
configuration files were used to run each of the different pipelines using the
``strun`` command. For example::

$ collect_pipeline_cfgs ./
$ strun calwebb_spec2.cfg an_exposure_file.fits

The issue with the above process is that any changes, as determined by INS and
the Calibration Working Group, to the default operation of the pipeline requires
a code release. A better solution would be if the pipeline configurations could
come from reference files retrieved from CRDS.

As of the version of master introduced on March 18th, 2021, in conjunction with
CRDS context jwst_0712, the default pipeline configurations no longer depend on
the package-delivered configuration files. Instead, all default configuration
relies on settings in the pipeline code itself, using CRDS-retrieved parameter
reference files to modify any parameters that are data-dependent. There is no
longer any need to run ``collect_pipeline_cfgs`` and specify a configuration
file for the ``strun`` command. One only needs to specify a simplified pipeline
name. In most cases, this simple name, or alias, is the same as the name of the
old configuration file, but without the suffix ``.cfg``.

Taking the example above, to get the same operation, the single command would become::

$ strun calwebb_spec2 an_exposure_file.fits

The JWST documentation has been updated to account for this change in usage. To
get familiarized, it is best to start with the :ref:`Introduction<introduction>`

A list of the available pipeline aliases can be found in the :ref:`Pipeline
Stages<pipelines>` section.

An added benefit to removing the dependency on package-delivered configuration
files is that users, under normal circumstances, no longer need to be concerned
with configuration files and whether they are up-to-date. One only needs to
install the JWST package and start using the pipelines out-of-the-box.

Does this mean that everyone has to immediately change their behavior and code
if using the default configuration files? Short answer is “No”. If one wishes to
continue using the package-delivered configuration files from
``collect_pipeline_cfgs``, one may do so. However, these configuration files no
longer contain any parameter settings; only the class name of the pipeline to be
run. This allows the code-plus-CRDS-retrieved parameter reference files to
determine operation.

Since the configuration settings have simply been moved to CRDS, the results one
obtains should not change. If a change in behavior is noted, please report the
issue to the Help Desk, file a Github issue on the JWST Github repository, or
file a Jira issue against the JP project.

In the meantime, please consider deprecating the use of ``collect_pipeline_cfgs``
and the .cfg files in favor of simply specifying pipeline aliases, as the
documentation now describes.

For users that use their own, custom configuration files, there is no change to
functionality. However, there are changes to both how these files are documented
and their format.

Concerning documentation, there is a change of terminology. No longer are these
files referred to as “configuration files”, but are called “parameter files” or
“parameter reference files” when retrieved from CRDS.

In order to simplify integration with CRDS, the format of parameter files have
changed from the “cfg”, init-like format, to the ASDF format. All parameter
files in CRDS are in this format. Similarly, the tools provided by the JWST
package to create parameter files will create them in ASDF. “cfg”-formatted
files are still supported, but it is strongly suggested that users change to
using the ASDF form. For more information, please to refer to :ref:`config_asdf_files`

As always, if anyone finds any discrepancies or other issues with the
documentation, or actual operation of the pipelines, please contact the Help
Desk, or file issues directly against the Github repository or the JIRA “JP”
project.
