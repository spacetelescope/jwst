<!-- If this PR closes a JIRA ticket, make sure the title starts with the JIRA issue number,
for example JP-1234: <Fix a bug> -->
Resolves [JP-nnnn](https://jira.stsci.edu/browse/JP-nnnn)

<!-- If this PR closes a GitHub issue, reference it here by its number -->
Closes #

<!-- describe the changes comprising this PR here -->
This PR addresses ...

**Checklist for PR authors (skip items if you don't have permissions or they are not applicable)**
- [ ] this PR changes user-facing code / API
  - [ ] added an entry to `CHANGES.rst` within the relevant release section
  - [ ] updated or added relevant tests
  - [ ] updated relevant documentation
- [ ] started a [regression test](https://github.com/spacetelescope/RegressionTests/actions/workflows/jwst.yml) and included a link to the running job
  - [ ] do truth files need to be updated ("okified")?
    - [ ] **after a maintainer has approved these changes**, run `okify_regtests` to update the truth files
- [ ] added a milestone, if relevant
- [ ] JIRA ticket is [resolved properly](https://github.com/spacetelescope/jwst/wiki/How-to-resolve-JIRA-issues)
