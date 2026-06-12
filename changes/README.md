# Writing News Fragments for the Changelog

This `changes/` directory contains "news fragments": small reStructuredText (`.rst`) files describing a change in a few sentences.
When making a release, run `towncrier build --version <VERSION>` to consume existing fragments in `changes/`
and insert them as a full changelog entry at the top of `CHANGES.rst` for the released version.

News fragment filenames consist of the pull request number and the changelog category (see below).
For example, https://github.com/spacetelescope/jwst/pull/10042 has two change notes: `tso_photometry` and `breaking`:

```
10028.other.rst
10042.breaking.rst
10042.tso_photometry.rst
10043.tweakreg.rst
10139.docs.rst
```

## Change Log Categories

Make a news fragment for every relevant category affected by your change.

- `<PR#>.breaking.rst`: changes that break [step-level or public API](https://jwst.readthedocs.io/en/latest/jwst/user_documentation/more_information.html#api-public-vs-private)
- `<PR#>.<stepmodulename>.rst`: changes to individual pipeline steps 

See [`towncrier.toml`](./towncrier.toml) for a full list of change log categories.

> [!NOTE]
> This README was adapted from a similar one in Astropy (under the terms of the BSD license),
> which was in turn adapted from Numpy (MIT license).
