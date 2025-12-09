# Contributing to the JWST Calibration Pipeline

`jwst` is an open source package written in Python.
The source code is available in the [JWST Github repository](https://github.com/spacetelescope/jwst/).
New contributions and contributors are very welcome!
Do not hesitate to reach out to the package maintainers if you are new to open-source development or if you have any questions/concerns.
We only ask that all contributors adhere to the Space Telescope [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting bugs / requesting a new feature
To submit a bug report or feature request, [open a new issue in this repository](https://github.com/spacetelescope/jwst/issues).

## Suggesting code changes / contributions

> [!TIP]
> If you are new to GitHub, to `git`, or to version-control systems in general, refer to the [GitHub tutorial](https://docs.github.com/en/get-started/git-basics/set-up-git) and / or to the [`git` reference manual](https://git-scm.com/docs).

To suggest a specific code change, or to contribute new code:
1. [Fork this repository](https://github.com/spacetelescope/jwst/fork).
2. Install `pre-commit` to automatically check your changes for formatting issues:
    ```shell
    pip install pre-commit
    pre-commit install
    ```
    > [!TIP]
    > To run `pre-commit` checks manually, do `pre-commit run --all`.
3. Clone your fork to your local machine:
    ```shell
    git clone https://github.com/YOURUSERNAME/jwst
    ```
    > [!TIP]
    > When making changes, it is standard practice to create a new "branch" for each new feature or bug fix.
    > We recommend naming your new branch something like `feature/cool_new_feature`, `fix/thing_that_was_fixed`, `docs/updated_description`, etc:
    > ```shell
    > git switch -c docs/update_contributing_instructions
    > ```
4. Make your changes.
5. Test your changes by running unit tests (see `TESTING.md` for details).
6. Commit and push your changes to your fork as a new branch:
    ```shell
    git commit -a "description of changes"
    git push
    ```
    The [`git` reference manual](https://git-scm.com/docs) has details on what these commands do.
7. [Open a new Pull Request](https://github.com/spacetelescope/jwst/pulls) requesting that your changes be merged into the `main` branch of this repository.
8. Complete the items in the **Tasks** checklist (created when you open the pull request) to the best of your ability.

Once your pull request is created, it will need to be reviewed and approved by the code maintainer team.
They may require changes from you before your code can be merged,
in which case you will need to go back and make those changes, run tests again, and push the changes to the branch you made earlier.

# Advanced Contribution Instructions

## Keeping your development branch current - rebasing

As `jwst` is constantly evolving, you will often encounter the situation where you've
made changes to your branch off `main`, but in the time its taken you to make those
changes, `upstream/main` has evolved with new commits from other developers. In this
situation, you will want to make sure you incorporate these changes into your branch.
Rebasing allows you to do two things - 1. apply others changes on top of yours, and 2.
squash your commits, even if there aren't new changes to apply.

First, make sure the `upstream` remote on your local repository is set to this repository:
```shell
git remote add upstream https://github.com/spacetelescope/jwst
```

Periodically, while writing code, to keep your branch up to date you will want to
do an interactive rebase against `upstream/main` to apply any new changes on top of yours.
```shell
git rebase -i upstream/main
```

This will then prompt you to select commits and commit messages - if you select
just the top commit, this will `squash` the others and combine them into one. You
will be prompted to resolve any conflicts if you've made modifications to the same
file someone else has. Once you've completed your rebase, you must `force push` your
branch to origin/my_feature to make your local and remote match.

```shell
git push -f feature/cool_new_feature
```

Before merging a PR, we typically would like you to rebase and squash all of your
commits into a single one, which is also done with `git rebase`.

## Writing and building documentation

`jwst` uses [sphinx](https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html)
to generate documentation, which is then hosted online on [readthedocs](https://readthedocs.org/).

You can access two versions of the documentation on the
[JWST readthedocs website](https://readthedocs.org/projects/jwst-pipeline/)
-- the `latest` version is whatever is currently on the main branch, and the `stable`
version is the last released version. If you successfully merge a PR with documentation
changes, they will only appear on `latest` until the next JWST release.

All of the documentation resides in the `jwst/docs` subdirectories (mainly within
directories in `jwst/docs/jwst`, organized by module). The documentation is written
in `.rst` (reStructured text) files - if you wish to make changes or add to the documentation,
do so in these files in the appropriate location. reStructured text is the markup
language used by sphinx, for information on the syntax refer to the
[`sphinx` documentation](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

While writing documentation, you will want to make sure the documentation builds
successfully  (i.e., produces a working HTML file). This happens on the CI when you
open a pull request, but it is a good idea to build the documentation yourself locally
as well. Before you do this, you will need to make sure you have the correct dependencies
installed in your current environment. All of these optional dependencies are specified
in `pyproject.toml` and include things like the correct version of sphinx, as well as the
necessary sphinx themes that the project uses. These do not install automatically
when you install `jwst` unless directly specified. To do this, while in the top level
directory of `jwst` on your `feature/cool_new_feature` branch:
```shell
pip install -e ".[docs]"
```

(Note the doc dependencies are also included when installing `jwst` with the `[contrib]` tag).
Now, with the correct documentation dependencies installed, you can attempt to build
the documentation locally. To do this, enter into the `jwst/docs` subdirectory and do:
```shell
make html
```

If the documentation successfully builds, the output HTML files will be output in
the `_build/html` subdirectory. You can open the main `index.html` file in your browser
and explore the full documentation pages just like the readthedocs site. If there were
any errors or warnings, a traceback will be printed on your terminal. Small errors, like
a typo in a link to another section, can cause frustrating errors so it is good practice
to build the docs frequently when editing them.

## Simultaneously developing `jwst` and one of its dependencies

If you encounter the scenario where you wish to simultaneously make changes in
`jwst` and also in one of its dependencies like `stcal` or `stpipe`, we recommend
that you create a new environment with development versions of both of these packages.
To do this, you can follow the same workflow outlined in the 'Contributing code' section
of this guide. To summarize, you will want to create a new Python environment
(called, for example, `jwst_stcal_dev`), fork and clone a local copy of `jwst` if you
haven't already and install this (doing `pip install -e .` in the top level directory
of `jwst`), and then fork and clone `stcal` and install it in this environment in
the same way. To double check that you have the correct dev versions of these packages
in your environment, you can check their versions by doing:
```shell
conda list
```

When opening up two dependent pull requests in `jwst` and one of its dependency packages,
unit tests will not pass on the CI because the `pyproject.toml` file in `jwst` points to the
last released version of `stcal`, and stcal points to the last version of `jwst`, so the
issue becomes circular. What you will need to do is modify the `pyproject.toml` files in both
packages to point to the other to demonstrate that CI tests pass (and make a comment
noting this in your PR), and then change it back before the PR is merge so that changes
to `pyproject.toml` are not merged into main. In your `jwst` branch, to point to your
branch in the dependent package (in this example `stcal`), change the required `stcal`
version in `pyproject.toml` to:
```
    stcal @  git+https://github.com/<your_username>/stcal.git@<your_branch>
```

And similarly, in `stcal`, change the required `jwst` version to:
```
    jwst @  git+https://github.com/<your_username>/jwst.git@<your_branch>
```

Let the CI run and ensure it passes, comment this in your PR and make sure the reviewers
confirm, and then change the versions back before your PR is merged
(which will again cause the CI to fail, but that’s OK).

## Code style

We use a pre-commit CI workflow to ensure that the code and docstring style of the `jwst` repository
remains uniform and conforms to certain standards. We recommend checking these using `pre-commit`
as described in the [Installing JWST for Development](Step-3-Installing-jwst-for-development)
section. For additional information about any of these individual style checkers,
see their documentation linked below.
Our pre-commit Git hook, also described in the
[Installing JWST for Development](Step-3-Installing-jwst-for-development) section,
is designed to help contributors run all the checks on their contributions every time they commit.

The following style checks are performed:

* **PEP8-compliant code**

  The code style for the `jwst` repository generally conforms to
  [PEP8](https://peps.python.org/pep-0008/), and the code style rules are enforced
  using [Ruff](https://docs.astral.sh/ruff/). To run these checks standalone,
  use the command:
    ```shell
    pre-commit run ruff
    ```
    
  from within the `jwst` repository. Ruff will automatically pick up the appropriate configuration
  from the `.ruff.toml` and `pre-commit-config.yaml` files,
  and perform only the checks that are turned on for our repository. To run ruff's
  auto-formatter, which automatically fixes simple things like single vs double quotes,
  whitespace, etc., use the command:
    ```shell
    pre-commit run ruff-format
    ```

* **Numpy docstring style**

  The docstring style for the `jwst` repository generally conforms to the
  [Numpy style guide](https://numpydoc.readthedocs.io/en/latest/format.html), and the docstring
  style rules are enforced using [numpydoc-validation](https://numpydoc.readthedocs.io/en/latest/validation.html).

  To run these checks standalone, use the command:
    ```shell
    pre-commit run numpydoc-validation
    ```

* **Spell checking**

  We use [Codespell](https://github.com/codespell-project/codespell) to check for common
  misspellings in both our codebase and documentation.
  To run the spell checker standalone, use the command:
    ```shell
    pre-commit run codespell
    ```

* **PEP-compliant type hints**

  The majority of the `jwst` repository does *not* have any type hints, and type hints are *not*
  required for contributions. If type hints are used, though, their compliance with
  [PEP-484](https://peps.python.org/pep-0484/) standards
  is enforced using [mypy](https://mypy.readthedocs.io/en/stable/index.html).
  To run these checks locally, use the command:
    ```shell
    pre-commit run mypy
    ```
