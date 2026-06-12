<!-- If this PR addresses a JIRA ticket: -->
<!-- Resolves [JP-nnnn](https://jira.stsci.edu/browse/JP-nnnn) -->

<!-- If this PR will close an existing GitHub issue (that is not already attached to a JIRA ticket): -->
<!-- Closes # -->

<!-- Describe your changes here: -->

## Description

This change ...

<!-- If you can't perform these tasks due to permissions, reach out to a maintainer. -->

## Tasks

- [ ] If you have a specific reviewer in mind, tag them.
- [ ] add a build milestone, i.e. `Build 12.0` (use the [latest build](https://github.com/spacetelescope/jwst/milestones) if not sure)
- [ ] update or add relevant tests
- [ ] update relevant docstrings and / or `docs/` page
- [ ] If this change affects user-facing code or public API, add news fragment file(s) to `changes/` (see [the changelog instructions](https://github.com/spacetelescope/jwst/blob/main/changes/README.md)).
      Otherwise, add the `no-changelog-entry-needed` label. - if your change breaks **step-level or public API** ([as defined in the docs](https://jwst.readthedocs.io/en/latest/jwst/user_documentation/more_information.html#api-public-vs-private)), also add a `changes/<PR#>.breaking.rst` news fragment
- [ ] [Start a regression test](https://github.com/spacetelescope/RegressionTests/actions/workflows/jwst.yml) and include a link to the running job ([click here for instructions](https://github.com/spacetelescope/RegressionTests/blob/main/docs/running_regression_tests.md))
- [ ] If truth files need to be updated ("okified"), then, **after the reviewer has approved these changes**, run `okify_regtests` to update the truth files.
- [ ] if a JIRA ticket exists, [make sure it is resolved properly](https://github.com/spacetelescope/jwst/wiki/How-to-resolve-JIRA-issues)

## Generative AI Disclosure

Were any generative AI or agentic LLMs used in the process of making this change?

- [ ] yes
- [ ] no
