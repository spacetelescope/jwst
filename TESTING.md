# Running tests for the JWST Calibration Pipeline

`jwst` has several test suites to ensure that functionality remains consistent and does not break when code changes.
In order for a change you make to the code to be accepted and merged, that change must pass existing tests, as well as any new tests you write that cover new functionality.

`jwst` uses `pytest` to define and run tests. To install `pytest` and other required testing tools to your [development environment](./CONTRIBUTING.md#creating-a-development-environment), install `jwst` with the `test` extra:
```shell
pip install -e .[test]
```

To run tests, simply run `pytest`:
```shell
pytest
```

`pytest` recursively searches the given directory (by default `.`) for any files with a name like `test_*.py`, and runs all functions it finds that have a name like `test_*`.

For example, running `pytest` without any arguments will find and run this test from `jwst/associations/tests/test_load_from_asn.py` (along with many others):
```python
from stdatamodels.jwst.datamodels import ImageModel

from jwst.associations.load_as_asn import LoadAsLevel2Asn


def test_lv2_datamodel():
    model = ImageModel()
    model.meta.filename = "modelfile.fits"
    asn = LoadAsLevel2Asn.load(model)
    assert asn.filename == DEFAULT_NAME
    assert asn["program"] == DEFAULT_NAME
    assert asn["target"] == DEFAULT_NAME
    assert asn["asn_pool"] == DEFAULT_NAME
    assert len(asn["products"]) == 1
    assert asn["products"][0]["name"] == "modelfile"
```

To run tests from a specific directory or file, pass the path to `pytest`:
```shell
# only run tests in `jwst/associations/`
pytest jwst/associations/

# only run tests in `jwst/associations/tests/test_load_as_asn.py`
pytest jwst/associations/tests/test_load_as_asn.py
```

To run only specific tests that contain a given string in their name, use `-k <string>`:
```shell
# only run tests with `datamodel` in their name
pytest -k datamodel

# only run the `test_lv2_datamodel` test
pytest -k test_lv2_datamodel
```

See the [`pytest` documentation](https://docs.pytest.org) for more instructions on using `pytest`.

> [!NOTE]
> By default, simply running `pytest` will skip some tests that require access to large datasets.
> 
> To run these tests, you must have access to the STScI internal network, set the environment variable `TEST_BIGDATA` to the STScI Artifactory server, and use the `--bigdata` flag with `pytest`:
> ```shell
> TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory pytest --bigdata
> ```

> [!NOTE]
> Additionally, some more tests are skipped by default because they take a long time.
> These can be run by using the `--slow` flag with `pytest`:
> ```shell
> pytest --slow
> ```

> [!TIP]
> You can control where test results are written by adding `--basetemp=<PATH>` to your `pytest` command.
> `pytest` will wipe this directory clean for each test session, so make sure it is a scratch area.

# Unit Tests

Unit tests are located in each module in `jwst/*/tests/` (for example `jwst/ramp_fitting/tests/`).
These tests run the code on simplified datasets to make sure there are no breaking changes introduced. Most lines of `jwst` should be covered
by a unit test, so when adding code you will often need to write a new test, or add to an existing test, to ensure adequate coverage.

Most test modules use a shared `@pytest.fixture` to set up a common dataset (usually a very simple `datamodel` with some observational parameters and simplified data and / or DQ arrays) for tests in that file.
Most individual tests set up a scenario to run a pipeline / step under certain conditions and test a specific functionality, culminating in a set of assertions (`assert`) that must be true in order for the test to pass.

When you open a pull request (assuming you have [set up your environment as described in `CONTRIBUTING.md`](./CONTRIBUTING.md#creating-a-development-environment)), GitHub will automatically start a workflow to run tests on your branch.
These tests re-run every time you `git push`.

## Writing Unit Tests

If your change introduces new code that isn't covered by an existing test, you should write new test(s) to cover this new functionality.

> [!TIP]
> If you are writing and / or debugging a test, you may find it helpful to have print statements printed to the terminal while running tests (which doesn't happen by default):
> ```shell
> pytest -s
> ```

Within the test files themselves, decorators can be used to control the behavior of the test.
Some of the more useful decorators include:

1. `@pytest.fixture` to declare a
   [fixture](https://docs.pytest.org/en/stable/explanation/fixtures.html)
2. `@pytest.mark.bigdata` (provided by `ci-watson`) to only run a test when the `--bigdata` flag is provided
3. `@pytest.mark.parametrize` to run a test multiple times over a set (or sets) of parameters
4. `@pytest.mark.skipif` to skip a test under certain conditions
5. `@pytest.mark.slow` (provided by `ci-watson`) to only run a test when the `--slow` flag is provided
6. `@pytest.mark.xfail` will make a test pass only if it fails (useful in certain circumstances)

# Regression Tests

Whereas unit tests focus on functionality on individual features with small datasets, regression tests in `jwst/regtest/` are focused mainly on large-scale functionality and thus require access to large datasets:
```shell
pytest --bigdata --slow jwst/regtest/
```

Regression tests are run on STScI-provisioned runners from the [RegressionTests repository (currently viewable by STScI staff only)](https://github.com/spacetelescope/RegressionTests/actions/workflows/jwst.yml?query=event%3Aschedule). To run regression tests, read [the instructions](https://github.com/spacetelescope/RegressionTests/blob/main/docs/running_regression_tests.md).

See [Maintaining Regression Tests](https://github.com/spacetelescope/jwst/wiki/Maintaining-Regression-Tests) for instructions on updating datasets.

