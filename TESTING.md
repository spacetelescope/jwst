# Running tests for the JWST Calibration Pipeline

When editing the `jwst` code locally, you can install the package as "editable" to immediately propagate any code changes into the installed package:
```shell
pip install -e .
```

This is useful for making changes and then testing them, either by running the pipeline normally or by running the tests from the test suite. 

## Running tests

`jwst` uses `pytest` to collect and run tests.
To install `pytest`, along with other dependencies required to run the unit tests,
install the package with the `test` extra included:
```shell
pip install -e ".[test]"
```

Run simple unit tests by invoking `pytest`:
```shell
pytest
```

By default, simply running `pytest` will skip some tests that require access to large datasets.

To run these tests, you must
1. have access to the STScI internal network
2. set the environment variable `TEST_BIGDATA` to the STScI Artifactory server:
  ```shell
  export TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory
  ```
3. add the `--bigdata` flag to your `pytest` command:
  ```
  pytest --bigdata
  ```

Additionally, some more tests are skipped by default because they take a long time.
These can be run by adding the `--slow` flag to your `pytest` command:
```shell
pytest --slow
```

> [!TIP]
> You can control where test results are written by adding `--basetemp=<PATH>` to your `pytest` command.
> `pytest` will wipe this directory clean for each test session, so make sure it is a scratch area.


## Regression Tests

Regression tests mainly require access to large datasets:
```shell
pytest --bigdata --slow jwst/regtest
```

Latest regression test results can be found [here (currently viewable by STScI staff only)](https://github.com/spacetelescope/RegressionTests/actions/workflows/jwst.yml?query=event%3Aschedule).

See [Maintaining Regression Tests](https://github.com/spacetelescope/jwst/wiki/Maintaining-Regression-Tests) for instructions on updating datasets.

## Writing and running unit tests

Unit tests are located in each module in `jwst` in a `tests` subdirectory (for example,
`jwst/ramp_fitting/tests`). These tests run the code on simplified datasets to make
sure there are no breaking changes introduced. Most lines of `jwst` should be covered
by a unit test, so when adding code you will often need to write a new test or add
to an existing test to ensure adequate coverage.

See [Writing and Running Unit Tests](#writing-and-running-unit-tests) for more details on testing.

Take a look around at the existing tests for a template - a majority of these tests
use a `@pytest.fixture` to set up a common dataset (usually a very simple `datamodel`
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
for running tests:
```shell
pip install -e ".[test]"
```

This will install the optional `test` dependencies specified in `pyproject.toml` that
don't install by default. (Note these test dependencies are also included when
installing `jwst` with the `[contrib]` tag).
The package `pytest` is one of these and is what's used
to run the tests. `pytest` searches through all the directories in your repository
(underneath the directory from which it was invoked command line) and looks for any
directories called `test` or `.py` files with the word `test` in the name. Functions
in these files will be executed.

To run all of the `jwst` unit tests, while in the `jwst` level directory, simply
run the command:
```shell
pytest
```

If you want to run all the tests for a single module, for example `ramp_fitting`,
you can run this from either the `jwst/ramp_fitting` OR the
`jwst/ramp_fitting/tests` directory.

To run all tests within a single test file (for example, all tests in
`jwst/jump/tests/test_detect_jumps.py`):
```shell
pytest jwst/jump/tests/test_detect_jumps.py
```

### Configuring pytest for unit tests

`pytest` can be configured with many different flags - see the
[`pytest` documentation](https://docs.pytest.org/en/stable/)
to see all of these. Here we will summarize a few of the most useful options.

If you are writing and debugging a test, you may find it helpful to have print statements
printed to the terminal while running tests (which doesn't happen by default). To do this:
```shell
pytest -s
```

If you want to only run a specific test function within a file (or only tests with names
that contain a certain string):
```shell
pytest -k testname
```

This will search all files under the directory you're in for files or functions with
the string `test` in them, and within those files run only the test functions that
contain the string `testname`.

Within the test files themselves, decorators can be used to control the behavior of the test.
Some of the more useful decorators include:

1. `@pytest.fixture` to declare a
   [fixture](https://docs.pytest.org/en/stable/explanation/fixtures.html)
2. `@pytest.mark.bigdata` (provided by `ci-watson`) to declare a test accessing
   large remote data from Artifactory and only run when `--bigdata` flag is provided
3. `@pytest.mark.parametrize` to run a single test on multiples sets
   of input parameters
4. `@pytest.mark.skip` to skip tests altogether, or under specific conditions
   (for example, only when being run by the CI)
5. `@pytest.mark.slow` (provided by `ci-watson`) to declare a test taking significant
   amount of time to run and only run when `--slow` flag is provided
6. `@pytest.mark.xfail` will make a test pass only if it fails

