# Writing and maintaining documentation

Documentation for `jwst` is written in [Sphinx reStructuredText (`.rst`)](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
in this `docs/` directory, and is hosted online at https://jwst.readthedocs.io### Building documentation locally

ReadTheDocs hosts several versions of the documentation.
To switch between versions, use the version selector:

- [`stable`](https://jwst-pipeline.readthedocs.io/en/stable/)
  is built from the last released version.
  If you successfully merge a PR with documentation changes,
  your changes will not be reflected in [`stable`](https://jwst-pipeline.readthedocs.io/en/stable/)
  until the next version of JWST is released.
- [`latest`](https://jwst-pipeline.readthedocs.io/en/latest/) is built from the latest commit on `main`.
- A static version is also built for each release.

ReadTheDocs will automatically build documentation for your branch when you push a commit to a pull request on this GitHub repository, and host a temporary build with a visual diff.
However, it is also good practice to build the docs locally if you are editing them, to reduce frustration from small errors.

To build the docs locally (assuming you have [set up your environment as described in `CONTRIBUTING.md`](../CONTRIBUTING.md#creating-a-development-environment)):

```shell
cd docs/
pip install ..[docs]
make clean
make html
```

The docs will build to `docs/_build/html/`.
Open `docs/_build/html/index.html` to view the pages in your browser.
