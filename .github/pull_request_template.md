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
- [ ] Does this PR change user-facing code / API?
  - [ ] add an entry to `CHANGES.rst` within the relevant release section (otherwise add the `no-changelog-entry-needed` label to this PR)
  - [ ] update or add relevant tests
  - [ ] update relevant docstrings and / or `docs/` page
  - [ ] [start a regression test](https://github.com/spacetelescope/RegressionTests/actions/workflows/jwst.yml) and include a link to the running job ([click here for instructions](https://github.com/spacetelescope/RegressionTests/blob/main/docs/running_regression_tests.md))
    - [ ] Do truth files need to be updated ("okified")?
      - [ ] **after the reviewer has approved these changes**, run `okify_regtests` to update the truth files
- [ ] if a JIRA ticket exists, [make sure it is resolved properly](https://github.com/spacetelescope/jwst/wiki/How-to-resolve-JIRA-issues)
