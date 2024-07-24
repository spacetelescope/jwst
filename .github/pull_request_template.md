<!-- If this PR closes a JIRA ticket, make sure the title starts with the JIRA issue number,
for example JP-1234: <Fix a bug> -->
Resolves [JP-nnnn](https://jira.stsci.edu/browse/JP-nnnn)

<!-- If this PR closes a GitHub issue, reference it here by its number -->
Closes #

<!-- describe the changes comprising this PR here -->
This PR addresses ...

### Required
- [ ] **request a review from someone specific so it preferably doesn't fall back onto the maintainer team**
- [ ] add a build milestone, i.e. `Build 11.3` (use the [latest build](https://github.com/spacetelescope/jwst/milestones) if not sure)

### Optional
- [ ] this PR changes user-facing code / API
  - [ ] add an entry to `CHANGES.rst` within the relevant release section (otherwise add the `no-changelog-entry-needed` label to this PR)
  - [ ] update or added relevant tests
  - [ ] update relevant docstrings and / or `docs/` page
  - [ ] [start a regression test](https://github.com/spacetelescope/RegressionTests/actions/workflows/jwst.yml) and include a link to the running job
    - [ ] do truth files need to be updated ("okified")?
      - [ ] **after a maintainer has approved these changes**, run `okify_regtests` to update the truth files
- [ ] if a JIRA ticket exists, [make sure it is resolved properly](https://github.com/spacetelescope/jwst/wiki/How-to-resolve-JIRA-issues)
