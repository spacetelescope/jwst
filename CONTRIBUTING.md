# Contributing to the JWST Calibration Pipeline

`jwst` is an open source package written in Python.
The source code is available in the [JWST Github repository](https://github.com/spacetelescope/jwst/).
New contributions and contributors are very welcome!
Do not hesitate to reach out to the package maintainers if you are new to open-source development or if you have any questions/concerns.
We only ask that all contributors adhere to the Space Telescope [Code of Conduct](./CODE_OF_CONDUCT.md).

## Reporting bugs / requesting a new feature
To submit a bug report or feature request, [open a new issue in this repository](https://github.com/spacetelescope/jwst/issues).

## Suggesting code changes / contributions

> [!TIP]
> If you are new to GitHub, to `git`, or to version-control systems in general, refer to the [GitHub tutorial](https://docs.github.com/en/get-started/git-basics/set-up-git) and / or to the [`git` reference manual](https://git-scm.com/docs).

To suggest a specific code change, or to contribute new code:
1. [Fork this repository](https://github.com/spacetelescope/jwst/fork).
2. Clone your fork to your local machine:
    ```shell
    git clone https://github.com/YOUR_USERNAME/jwst
    cd jwst/
    ```

3. Add the `upstream` repository, as a remote, to your local clone:
    ```shell
    git remote add upstream https://github.com/spacetelescope/jwst
    ```

> [!TIP]
> When making changes, create a new "branch" for each new feature or bug fix.
> We recommend naming your new branch something like `feature/cool_new_feature`, `fix/thing_that_was_fixed`, `docs/updated_description_of_feature`, etc:
> ```shell
> git fetch upstream --tags
> git checkout upstream/main -b docs/update_contributing_instructions
> ```

4. Install `pre-commit` to automatically check your changes for formatting issues:
    ```shell
    pip install pre-commit
    pre-commit install
    ```

> [!TIP]
> To run `pre-commit` checks manually, do `pre-commit run --all`.

5. [Install `jwst` to your development environment.](#creating-a-development-environment)
6. Make your changes using your editor of choice.
7. Commit and push your changes to your fork as a new branch:
    ```shell
    git add changed_file.py
    git commit -m "description of changes"
    git push
    ```
    The [`git` reference manual](https://git-scm.com/docs) has details on what these commands do.
8. [Open a new Pull Request](https://github.com/spacetelescope/jwst/pulls) requesting that your changes be merged into the `main` branch of this repository.
9. Ensure that your change passes automated testing (see [`TESTING.md`](./TESTING.md) for details).
10. Complete the items in the **Tasks** checklist (created when you open the pull request) to the best of your ability.

Once your pull request is created, it will need to be reviewed and approved by the code maintainer team.
They may require changes from you before your code can be merged,
in which case go back and make those changes, run tests again, and push the changes to the branch you made earlier.

## Keeping your development branch current with `main`

As `jwst` is constantly evolving, you will often encounter the situation where you've made changes to your branch, but in that time there are new commits on `upstream/main` from other developers.
Incorporate those changes into your branch, either automatically with the button on the GitHub pull request webpage, or manually with `git rebase`.

### Incorporate upstream changes automatically with button on GitHub pull request page

Usually, GitHub can rebase a branch automatically.
If you see "**This branch is out-of-date with the base branch**", you will have the option to "**Update with merge commit**" or "**Update with rebase**". Updating with a merge commit is usually safer.

However, if the changes to `main` touch the same lines as your changes, you will see "**This branch has conflicts that must be resolved**". You will need to manually resolve these conflicts yourself; follow the steps described on the page.

### Incorporate upstream changes manually with `git rebase`

Rebase your current branch onto `upstream/main` to apply any new changes on top of yours:
```shell
git fetch --all
git rebase -i upstream/main
```

For more information on how to use `git rebase`, see [the `git rebase` documentation](https://git-scm.com/docs/git-rebase) or [Atlassian's tutorial on rebasing](https://www.atlassian.com/git/tutorials/rewriting-history/git-rebase).

Once you've completed your rebase, you will need to "force push" your branch to **overwrite your branch on GitHub**:

```shell
git push -u origin -f feature/cool_new_feature
```

## Creating a development environment

When developing `jwst` (or any other Python package), you should install the package locally to a development environment.

> [!TIP]
> Python "environments" are isolated Python installations, confined to a single directory, where you can install packages, dependencies, and tools without cluttering your system Python libraries.

The easiest way to create a development environment is with `virtualenv`:
```shell
virtualenv .venv/
source .venv/bin/activate
pip install -e .
hx .
```

Breaking down what these lines do:
1. Create a new empty Python environment in the `.venv/` directory:
   ```shell
   virtualenv .venv/
   ```
2. "Activate" the environment (change shell variables in the current session to point to the isolated Python installation):
   ```shell
   source .venv/bin/activate
   ```
3. Install the local package (`jwst`) to your environment in "editable mode", so that any code changes will be instantly reflected in the installed package (useful for testing):
   ```shell
   pip install -e . 
   ```
4. Run your editor of choice (in this example I use Helix `hx`):
   ```shell
   hx .
   ```

> [!TIP]
> There are other ways of managing environments.
> For instance, if you have `uv` installed, you can accomplish the same as above with a single command:
> ```shell
> uv run hx .
> ```

## Making simultaneous changes to `jwst` and one of its dependencies

If you need to make a change in `jwst` that requires a simultaneous change to one of its dependencies (`stcal`, `stdatamodels`, `stpipe`, etc.), also install that dependency from your local machine [to your development environment](#creating-a-development-environment).
For instance, assuming you've cloned the source code for both `jwst` and `stcal`, you can do the following from the `jwst/` directory:
```shell
cd jwst/
pip install -e .
pip install -e ../stcal
```

> [!TIP]
> It might be easier to use a separate Python environment (`virtualenv`, `uv`, `conda`, etc.) for this work:
> ```shell
> virtualenv .venv_jwst_stcal_dev
> source .venv_jwst_stcal_dev/bin/activate
> ```

Since we do not use a single repository (sometimes called a "monorepo") for these coupled packages, when making a change like this, you will need to make two pull requests: one in the `jwst` repository, and one in the dependency's repository.
However, unit tests will not automatically pass because the `pyproject.toml` file in `jwst` points to the last released version of `stcal`, which does not incorporate your changes.
To resolve this (temporarily, for testing) modify the `pyproject.toml` files in `jwst` to point to your branch, to demonstrate that unit tests pass. 
```toml
# jwst/pyproject.toml
[project]
...
dependencies = [
    ...
    "stcal @ git+https://github.com/YOUR_USERNAME/stcal.git@YOUR_BRANCH",
    ...
]
```
> [!WARNING]
> **REMEMBER TO REVERT THE CHANGES TO `pyproject.toml` BEFORE YOUR BRANCH IS MERGED INTO `main`.**

## Code style

We use `pre-commit` to enforce uniform code style and standards.

```shell
pip install pre-commit
pre-commit run
```

You can also install `pre-commit` locally, to run checks before every `git commit` action:
```shell
pre-commit install
```

The full configuration for `pre-commit` checks can be found in `.pre-commit-config.yaml`.


### PEP8 compliance

The code style for the `jwst` repository generally conforms to [PEP8](https://peps.python.org/pep-0008/), enforced using [`ruff`](https://docs.astral.sh/ruff/).
`ruff` will automatically pick up the appropriate configuration from `.ruff.toml` and perform only the checks that are turned on for our repository.

### Numpy docstring style

Docstrings in `jwst` conform to the [Numpy style guide](https://numpydoc.readthedocs.io/en/latest/format.html), enforced with [`numpydoc-validation`](https://numpydoc.readthedocs.io/en/latest/validation.html).

### Spell checking

We use [Codespell](https://github.com/codespell-project/codespell) to check for common misspellings in both our codebase and documentation.

### PEP-compliant type hints

The majority of the `jwst` repository does *not* have type hints, and type hints are *not* required for contributions. If type hints are used, though, their compliance with [PEP-484](https://peps.python.org/pep-0484/) is enforced with [`mypy`](https://mypy.readthedocs.io/en/stable/index.html).

## Writing and maintaining documentation

See [`docs/README.md`](./docs/README.md) for instructions.
