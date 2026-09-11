# Writing News Fragments for the Changelog

This `changes/` directory contains "news fragments": small reStructuredText (`.rst`) files describing a change in a few sentences.
When making a release, run `towncrier build --version <VERSION>` to consume existing fragments in `changes/`
and insert them as a full changelog entry at the top of [`CHANGES.rst`](../CHANGES.rst) for the released version.

News fragment filenames consist of the pull request number and the changelog category.
Make a news fragment for every relevant category affected by your change.
A single change can have more than one news fragment, if it spans multiple categories.
Categories are listed in [`towncrier.toml`](../towncrier.toml).
