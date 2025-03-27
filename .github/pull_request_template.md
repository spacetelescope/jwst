<!-- If this PR closes a JIRA ticket, make sure the title starts with the JIRA issue number,
for example JP-1234: <Fix a bug> -->
Resolves [JP-nnnn](https://jira.stsci.edu/browse/JP-nnnn)

<!-- If this PR closes a GitHub issue, reference it here by its number -->
Closes #

<!-- describe the changes comprising this PR here -->
This PR addresses ...

<!-- if you can't perform these tasks due to permissions, please ask a maintainer to do them -->
## Tasks
- [ ] **request a review from someone specific**, to avoid making the maintainers review every PR
- [ ] add a build milestone, i.e. `Build 11.3` (use the [latest build](https://github.com/spacetelescope/jwst/milestones) if not sure)
- [ ] Does this PR change user-facing code / API? (if not, label with `no-changelog-entry-needed`)
  - [ ] write news fragment(s) in `changes/`: `echo "changed something" > changes/<PR#>.<changetype>.rst` (see below for change types) 
  - [ ] update or add relevant tests
  - [ ] update relevant docstrings and / or `docs/` page
  - [ ] [start a regression test](https://github.com/spacetelescope/RegressionTests/actions/workflows/jwst.yml) and include a link to the running job ([click here for instructions](https://github.com/spacetelescope/RegressionTests/blob/main/docs/running_regression_tests.md))
    - [ ] Do truth files need to be updated ("okified")?
      - [ ] **after the reviewer has approved these changes**, run `okify_regtests` to update the truth files
- [ ] if a JIRA ticket exists, [make sure it is resolved properly](https://github.com/spacetelescope/jwst/wiki/How-to-resolve-JIRA-issues)

<details><summary>news fragment change types...</summary>

- ``changes/<PR#>.general.rst``: infrastructure or miscellaneous change
- ``changes/<PR#>.docs.rst``
- ``changes/<PR#>.stpipe.rst``
- ``changes/<PR#>.datamodels.rst``
- ``changes/<PR#>.scripts.rst``
- ``changes/<PR#>.set_telescope_pointing.rst``
- ``changes/<PR#>.pipeline.rst``

## stage 1
- ``changes/<PR#>.group_scale.rst``
- ``changes/<PR#>.dq_init.rst``
- ``changes/<PR#>.emicorr.rst``
- ``changes/<PR#>.saturation.rst``
- ``changes/<PR#>.ipc.rst``
- ``changes/<PR#>.firstframe.rst``
- ``changes/<PR#>.lastframe.rst``
- ``changes/<PR#>.reset.rst``
- ``changes/<PR#>.superbias.rst``
- ``changes/<PR#>.refpix.rst``
- ``changes/<PR#>.linearity.rst``
- ``changes/<PR#>.rscd.rst``
- ``changes/<PR#>.persistence.rst``
- ``changes/<PR#>.dark_current.rst``
- ``changes/<PR#>.charge_migration.rst``
- ``changes/<PR#>.jump.rst``
- ``changes/<PR#>.clean_flicker_noise.rst``
- ``changes/<PR#>.ramp_fitting.rst``
- ``changes/<PR#>.gain_scale.rst``

## stage 2
- ``changes/<PR#>.assign_wcs.rst``
- ``changes/<PR#>.badpix_selfcal.rst``
- ``changes/<PR#>.msaflagopen.rst``
- ``changes/<PR#>.nsclean.rst``
- ``changes/<PR#>.imprint.rst``
- ``changes/<PR#>.background.rst``
- ``changes/<PR#>.extract_2d.rst``
- ``changes/<PR#>.master_background.rst``
- ``changes/<PR#>.wavecorr.rst``
- ``changes/<PR#>.srctype.rst``
- ``changes/<PR#>.straylight.rst``
- ``changes/<PR#>.wfss_contam.rst``
- ``changes/<PR#>.flatfield.rst``
- ``changes/<PR#>.fringe.rst``
- ``changes/<PR#>.pathloss.rst``
- ``changes/<PR#>.barshadow.rst``
- ``changes/<PR#>.photom.rst``
- ``changes/<PR#>.pixel_replace.rst``
- ``changes/<PR#>.resample_spec.rst``
- ``changes/<PR#>.residual_fringe.rst``
- ``changes/<PR#>.cube_build.rst``
- ``changes/<PR#>.extract_1d.rst``
- ``changes/<PR#>.resample.rst``

## stage 3
- ``changes/<PR#>.assign_mtwcs.rst``
- ``changes/<PR#>.mrs_imatch.rst``
- ``changes/<PR#>.tweakreg.rst``
- ``changes/<PR#>.skymatch.rst``
- ``changes/<PR#>.exp_to_source.rst``
- ``changes/<PR#>.outlier_detection.rst``
- ``changes/<PR#>.tso_photometry.rst``
- ``changes/<PR#>.stack_refs.rst``
- ``changes/<PR#>.align_refs.rst``
- ``changes/<PR#>.klip.rst``
- ``changes/<PR#>.spectral_leak.rst``
- ``changes/<PR#>.source_catalog.rst``
- ``changes/<PR#>.combine_1d.rst``
- ``changes/<PR#>.ami.rst``

## other
- ``changes/<PR#>.wfs_combine.rst``
- ``changes/<PR#>.white_light.rst``
- ``changes/<PR#>.cube_skymatch.rst``
- ``changes/<PR#>.engdb_tools.rst``
- ``changes/<PR#>.guider_cds.rst``
</details>
