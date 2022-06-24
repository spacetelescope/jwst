# Basic Contribution Instructions

`jwst` is an open source python package, and the source code is available in the 
[JWST Github repository](https://github.com/spacetelescope/jwst/). New contributions
and contributors are very welcome! Do not hesitate to reach out to the package
maintainers if you are new to open-source development
or if you have any questions/concerns. We just ask that all contributors adhere
to the Space Telescope [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting bugs / requesting a new feature

If you would like to report a bug or request a new feature, this can be done by
opening a new [issue](https://github.com/spacetelescope/jwst/issues).

## Contributing code

If you would like to contribute code, this is done by submitting a [pull request](https://github.com/spacetelescope/jwst/pulls)
to the "master" branch of `spacetelescope/jwst`. To do this, we recommend the
following workflow (which assumes you already have a Github account / command line tools).
If you are also new to git, please refer to the [git reference manual](https://git-scm.com/docs)
for an overview of git basics.

### Step 1: Forking and cloning the `jwst` repository

First, to clarify some terms that will be commonly used here:

* `upstream` refers to the main `jwst` repository. this is where code from all
	contributors is ultimately merged into and where releases of the package will be made from.
* `origin` refers to the online fork you made of the upstream `jwst` repository.
* `local` refers to the clone you made of the origin on your computer.
* The term `remote` in this context can refer to origin, or upstream. in general, it means anything hosted online on Github.

The first step is to create your own 'remote' (online) and 'local' (on your machine)
clones of the central `spacetelescope/jwst` repository. You will make code changes
on your machine to your 'local' clone, push these to 'origin' (your online fork),
and finally, open a pull request to the ''master'' branch of `spacetelescope/jwst`.

1. On the 'spacetelescope/jwst' Github repository page, 'fork' the JWST repository
to your own account space by clicking the appropriate button on the upper right-hand
side. This will create an online clone of the main JWST repository but instead of
being under the 'spacetelescope' organization, it will be under your personal
account space. It is necessary to fork the repository so that many contributors
aren't making branches directly on 'spacetelescope/jwst'. 

2. Now that you have remotely forked `jwst`, it needs to be downloaded
to your machine. To create this 'local' clone, choose an area on your file system
and use the `git clone` command to dowload your remote fork on to your machine.

		>> cd directory
		>> git clone git@github.com:<your_username>/jwst.git

3. Make sure that your references to 'origin' and 'upstream' are set correctly - you will
need this to keep everything in sync and push your changes online. While your inital
local clone will be an exact copy of your remote, which is an exact copy of the 'upstream'
`spacetelescope/jwst`, these all must be kept in sync manually (via git fetch/pull/push).

	To check the current status of these references:

		>> git remote -v

After your inital clone, you will likely be missing the reference to 'upstream'
(which is just the most commonly used name in git to refer to the main project repository - you
can call this whatever you want but the origin/upstream conventions are most commonly used) - to 
set this, use the `add` git command:

If you are using an SSH key to authenticate.

	>> git remote add upstream git@github.com:spacetelescope/jwst.git
Otherwise, you can simply set it to the repository URL but you will have to
enter your password every time you fetch/push

	>> git remote add upstream https://github.com/spacetelescope/jwst

If you ever want to reset these URLs, add references to other remote forks of
`jwst` for collaboration, or change from URL to SSH, you can use the related
`git remote set-url` command.

### Step 2: Creating a branch for your changes

It is a standard practice in git to create a new 'branch' (off `upstream/master`)
for each new feature or bug fix. You can call this branch whatever you like - in
this example, we'll call it 'my_feature'. First, make sure you
have all recent changes to upstream by 'fetching' them:

		>> git fetch upstream

The following will create a new branch off local/master called 'my_feature', and automatically switch you over to your new branch.

		>> git checkout -b my_feature upstream/master

### Step 3: Installing `jwst` for development

Next, you will create a new Python environment where you will install your new
branch of `jwst`. We will show here how to use conda to manage environments, but
there are other options.

1.  Create a conda environment.

It is good practice to maintain different environments for different versions
of JWST and its dependencies. You will likely want to maintain one, for example,
for the latest released version of JWST (i.e. what you get by doing `pip install jwst`),
as well as one for development. Assuming the user has conda [installed](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html),
here we will create a new conda environment called 'jwst_dev' where we can install the new branch of our cloned repository.

	>> conda create -n jwst_dev python


Doing this will create a new environment with just some basic packages
(i.e setuptools, pip) installed.

2. Installing `jwst` in this environment

Make sure you are in your new environment:

	>> conda activate jwst_dev

And now in the top level of your local `jwst` repository, ensuring you're
on the 'my_feature' branch:

	>> pip install -e .
This will install `jwst` from this cloned source code in 'editable' mode,
meaning that you can import the code from this directory when within a Python
session. This makes it easier for development because you can have the code
you're editing somewhere convenient in your file system vs. with other packages
in 'site-packages'. If you cloned the repository on your Desktop, for example,
you can modify it there and Python will know that is where the source code is
when you're importing it within a Python session.

*Note* : If you use it, make sure to install iPython in your new environment
as well. Otherwise, it will pick up packages from the base environment instead.

### Step 4: Making code changes

Now that you've forked, cloned, made a new branch for your feature, and installed
it in a new environment for development of `jwst`, you are ready to make changes
to the code. As you make changes, make sure to `git commit -m <"some message">` frequently
(in case you need to undo something by reverting back to a previous commit - you
cant do this if you commit everything at once!). After you've made your desired
changes, and committed these changes, you will need to push them online to your
'remote' fork of `jwst`:

	>> git push origin my_feature

If the changes are significant, please make an entry in `CHANGES.rst` in the top
level `jwst` directory with a short description of the changes you've made and, once
you open a pull request, add the corresponding PR number.



### Step 4: Opening a pull request

Now, you can open a pull request on the master branch of the upstream `jwst` repository.

1. On the `spacetelescope/jwst` web page, after you push your changes you should
see a large green banner appear at the top prompting you to open a pull request
with your recently pushed changes. You can also open a pull request from the
[pull request tab](https://github.com/spacetelescope/jwst/pulls) on that page.
Select your fork and your 'my_feature' branch, and open a pull request against
the 'master' branch.

2. There is now a checklist of items that need to be done before your PR can be merged.
	* The continuous integration (CI) tests must complete and pass. The CI
	runs several different checks including running the unit tests, ensuring
	the documentation builds, checking for code style issues (see the [PEP8](https://peps.python.org/pep-0008/) style guide),
	and ensuring any changes are covered by unit tests. The CI runs upon opening
	a PR, and will re-run any time you push commits to that branch. 
	* You will need to add a change log entry in CHANGES.rst if your contribution
	is a new feature or bug fix. An entry is not required for small fixes like typos.
	* Your PR will need to be reviewed and approved by at least two maintainers.
	They may require changes from you before your code can be merged, in which
	case you will need to go back and make these changes and push them (they will
		automatically appear in the PR when they're pushed to origin/my_feature).


# Advanced Contribution Instructions

## Keeping your development branch current - rebasing

As `jwst` is constantly evolving, you will often encounter the situation where you've
made changes to your branch off 'master', but in the time its taken you to make those
changes, 'upstream/master' has evolved with new commits from other developers. In this
situation, you will want to make sure you incorporate these changes into your branch.
Rebasing allows you to do two things - 1. apply others changes on top of yours, and 2.
squash your commits, even if there aren't new changes to apply. 

Periodically, while writing code, to keep your branch up to date you will want to
do an interactive rebase against upstream/master to apply any new changes on top of yours:

	>> git rebase -i upstream/master

This will then prompt you to select commits and commit messages - if you select
just the top commit, this will 'squash' the others and combine them into one. You
will be prompted to resolve any conflicts if you've made modifications to the same
file someone else has. Once you've completed your rebase, you must `force push` your
branch to origin/my_feature to make your local and remote match.

	>> git push -f origin/my_feature

Before merging a PR, we typically would like you to rebase and squash all of your
commits into a single one, which is also done with `git rebase`

## Writing and building documentation

`jwst` uses [sphinx](https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html) to generate documentation, which is then hosted online on [readthedocs](https://readthedocs.org/).

You can access two versions of the documentation on the [JWST readthedocs website](https://readthedocs.org/projects/jwst-pipeline/)
- the 'latest' version is whatever is currently on the master branch, and the 'stable'
version is the last released version. If you successfully merge a PR with documentation
changes, they will only appear on 'latest' until the next JWST release.

All of the documentation resides in the `jwst/docs` subdirectories (mainly within
directories in `jwst/docs/jwst`, organized by module). The documentation is written
in `.rst` (reStructured text) files - if you wish to make changes or add to the documentation,
do so in these files in the appropriate location. reStructured text is the markup
language used by sphinx, for information on the syntax refer to the [sphinx documentation](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

While writing documentation, you will want to make sure the documentation builds
successfully  (i.e, produces a working .html file). This happens on the CI when you
open a pull request, but it is a good idea to build the documentation yourself locally
as well. Before you do this, you will need to make sure you have the correct dependencies
installed in your current environment. All of these optional dependencies are specified
in `setup.cfg` and include things like the correct version of sphinx, as well as the
necessary sphinx themes that the project uses. These do not install automatically
when you install `jwst` unless directly specified. To do this, while in the top level
directory of `jwst` on your my_feature branch:

	>> pip install -e ".[docs]"

Now, with the correct documentation dependencies installed, you can attempt to build
the documentation locally. To do this, enter into the `jwst/docs` subdirectory and do:

	>> make html

If the documentation successfully builds, the output HTML files will be output in
the `_build/html` subdirectory. You can open the main `index.html` file in your browser
and explore the full documentation pages just like the readthedocs site. If there were
any errors or warnings, a traceback will be printed on your terminal. Small errors, like
a typo in a link to another section, can cause frustrating errors so it is good practice
to build the docs frequently when editing them.

## Writing and running unit tests

Unit tests are located in each module in `jwst` in a `tests` subdirectory (for example,
`jwst/jwst/ramp_fitting/tests`). These tests run the code on simplified datasets to make
sure there are no breaking changes introduced. Most lines of `jwst` should be covered
by a unit test, so when adding code you will often need to write a new test or add
to an existing test to ensure adequate coverage.

Take a look around at the existing tests for a template - a majority of these tests
use a @pytest.fixture to set up a common dataset (usually a very simple `datamodel`
with some observational parameters and simplified data/dq arrays) for tests in that
file, and the test functions themselves set up a scenario to run a pipeline/step
under certain conditions, and culminate in a set of assertions that need to pass
(or fail if the test is marked as `xfail`).

The CI will run the unit tests on your branch when you open a pull request. They
will re-run every time you push commits to this remote branch as well. Unit tests
must pass on any changes you make, and if you're introducing new code that isn't
covered by an existing test, you will need to write one. The `codecoverage` check
that runs on the CI when you open a pull request will tell you if you've introduced
any new lines that aren't covered by a test.

You can also run unit tests locally, and you should do this periodically to test
your code. To do this, you will need to install the optional dependencies needed
for running tests.

	>> pip install -e ".[test]"

This will install the optional 'test' dependencies specified in `setup.cfg` that
don't install by default. The package `pytest` is one of these and is what's used
to run the tests. `pytest` searches through all the directories in your repository
(underneath the directory from which it was invoked command line) and looks for any
directories called 'test' or .py files with the word 'test' in the name. Functions
in these files will be executed.

To run all of the `jwst` unit tests, while in the `jwst/jwst` level directory, simply
run the command:

	>> pytest

If you want to run all the tests for a single module, for example `ramp_fitting`,
you can run this from either the 'jwst/jwst/ramp_fitting' OR the 'jwst/jwst/ramp_fitting/tests' directory/.

To run all tests within a single test file (for example, all tests in `jwst/jwst/jump/test_detect_jumps.py`).
	>> pytest test_detect_jumps.py

### Configuring pytest for unit tests

`pytest` can be configured with many different flags - see the [pytest documentation](https://docs.pytest.org/en/7.1.x/contents.html)
to see all of these. Here we will summarize a few of the most useful options.

If you are writing and debugging a test, you may find it helpful to have print statements
printed to the terminal while running tests (which doesn't happen by default). To do this,

	>> pytest -s

If you want to only run a specific test function within a file (or only tests with names
that contain a certain string):

	>> pytest -k testname.
This will search all files under the directory you're in for files or functions with
the string 'test' in them, and within those files run only the test functions that
contain the string `testname`.

Within the test files themselves, decorators can be used to control the behavior of the test. Some of the more useful decorators include:
	1. @pytest.mark.parametrize can be used to run a single test on multiples sets of input parameters
	2. @pytest.skip can be used to skip tests altogether, or under specific conditions (for example, only when being run by the CI)
	3. @pytest.fixture to declare a [fixture](https://docs.pytest.org/en/7.1.x/explanation/fixtures.html)
	4. @pytest.mark.xfail will make a test pass only if it fails.


## Simultaneously developing `jwst` and one of its dependencies

If you encounter the scenario where you wish to simultaneously make changes in
`jwst` and also in one of its dependencies like `stcal` or `stpipe`, we recommend
that you create a new environment with development versions of both of these packages.
To do this, you can follow the same workflow outlined in the 'Contributing code' section
of this guide. To summarize, you will want to create a new Python environment
(called, for example, jwst_stcal_dev), fork and clone a local copy of `jwst` if you
haven't already and install this (doing `pip install -e .` in the top level directory
of `jwst`), and then fork and clone `stcal` and install it in this environment in
the same way. To double check that you have the correct dev versions of these packages
in your environment, you can check their versions by doing:

	>> conda list
When opening up two dependent pull requests in `jwst` and one of its dependency packages,
unit tests will not pass on the CI because the setup.cfg file in `jwst` points to the
last released version of `stcal`, and stcal points to the last version of `jwst`, so the
issue becomes circular. What you will need to do is modify the setup.cfg files in both
packages to point to the other to demonstrate that CI tests pass (and make a comment
noting this in your PR), and then change it back before the PR is merge so that changes
to setup.cfg are not merged into master/main. In your `jwst` branch, to point to your
branch in the dependent package (in this example `stcal`), change the required `stcal`
version in setup.cfg to:

	>> stcal @  git+https://github.com/<your_username>/stcal.git@<your_branch>

And similarly, in `stcal`, change the required `jwst` version to:

	>> jwst @  git+https://github.com/<your_username>/jwst.git@<your_branch>
Let the CI run and ensure it passes, comment this in your PR and make sure the reviewers
confirm, and then change the versions back before your PR is merged (which will again cause the CI to fail, but that’s OK).
