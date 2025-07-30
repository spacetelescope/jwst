import json
from pathlib import Path

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.coron.tests.helpers import psf_nircam, target_nircam
from jwst.pipeline.calwebb_coron3 import Coron3Pipeline, to_container
from jwst.stpipe import Step


# Generate data
def make_container():
    """Create the container to test"""
    size = 5
    cube = datamodels.CubeModel((size, size, size))
    cube.meta.target.proposer_name = "JWST regression test"
    container = to_container(cube)
    return cube, container


class MockResample(Step):
    def process(self, input_data):
        cube = target_nircam()
        image = datamodels.ImageModel()
        image.data = cube.data[0]
        image.err = cube.err[0]
        image.dq = cube.dq[0]
        image.update(cube)
        cube.close()

        image.meta.cal_step.resample = "COMPLETE"
        return image


@pytest.fixture(scope="module")
def cube_and_container():
    cube, container = make_container()
    return cube, container


@pytest.fixture(scope="module")
def cube_models(cube_and_container):
    cube, container = cube_and_container
    models = [(cube, model) for model in container]
    return models


@pytest.fixture()
def coron3_input(tmp_path):
    n_model = 3
    asn = {
        "asn_id": "c1000",
        "target": "t001",
        "asn_pool": "test_coron3_pool.csv",
        "products": [{"name": "test_coron3", "members": []}],
    }
    for i in range(n_model):
        target_name = tmp_path / f"test_coron3_target_{i + 1}.fits"
        target = target_nircam()
        target.save(target_name)
        asn["products"][0]["members"].append({"expname": target_name.name, "exptype": "science"})

        psf_name = tmp_path / f"test_coron3_psf_{i + 1}.fits"
        psf = psf_nircam()
        psf.save(psf_name)
        asn["products"][0]["members"].append({"expname": psf_name.name, "exptype": "psf"})

    asn_name = tmp_path / "test_coron3_asn.json"
    with asn_name.open("w") as fh:
        json.dump(asn, fh)
    return asn_name


@pytest.fixture()
def coron3_science_only(tmp_path):
    asn = {
        "asn_id": "c1000",
        "target": "t001",
        "asn_pool": "test_coron3_pool.csv",
        "products": [{"name": "test_coron3", "members": []}],
    }
    target_name = tmp_path / "test_coron3_target_1.fits"
    target = target_nircam()
    target.save(target_name)
    asn["products"][0]["members"].append({"expname": target_name.name, "exptype": "science"})

    asn_name = tmp_path / "test_coron3_asn.json"
    with asn_name.open("w") as fh:
        json.dump(asn, fh)
    return asn_name


@pytest.fixture()
def coron3_psf_only(tmp_path):
    asn = {
        "asn_id": "c1000",
        "target": "t001",
        "asn_pool": "test_coron3_pool.csv",
        "products": [{"name": "test_coron3", "members": []}],
    }
    psf_name = tmp_path / "test_coron3_psf_1.fits"
    psf = psf_nircam()
    psf.save(psf_name)
    asn["products"][0]["members"].append({"expname": psf_name.name, "exptype": "psf"})

    asn_name = tmp_path / "test_coron3_asn.json"
    with asn_name.open("w") as fh:
        json.dump(asn, fh)
    return asn_name


def test_to_container():
    """Cover bug where IndexError would be raised when area extension of CubeModel
    has shape (x,y) instead of (nints,x,y). In this case area extension should
    be copied to each ImageModel in the ModelContainer.
    """
    shp = (10, 5, 5)
    cube = datamodels.CubeModel(shp)
    extensions_3d = [
        "data",
        "dq",
        "err",
    ]
    for extension in extensions_3d:
        setattr(cube, extension, np.random.rand(*shp))
    cube.area = np.random.rand(5, 5)

    container = to_container(cube)
    for i, model in enumerate(container):
        for extension in extensions_3d:
            assert np.all(getattr(model, extension) == getattr(cube, extension)[i])
        assert hasattr(model, "area")
        assert np.all(model.area == cube.area)


def test_to_container_wrong_dimensions():
    model = datamodels.CubeModel((10, 10, 10))
    model.zeroframe = np.zeros((10, 10, 10, 10))
    with pytest.raises(ValueError, match="Unexpected array shape"):
        to_container(model)


def test_container_shape(cube_and_container):
    """Test container shape"""
    cube, container = cube_and_container
    assert len(container) == cube.shape[0]


def test_meta(cube_models):
    """Test meta equivalency"""
    for model, cube in cube_models:
        assert model.meta.target.proposer_name == cube.meta.target.proposer_name


@pytest.mark.parametrize("array", ["data", "dq", "err"])
def test_shape(cube_models, array):
    """Test array shapes"""
    for cube, model in cube_models:
        assert model[array].shape == cube[array][0].shape


@pytest.mark.parametrize("array", ["zeroframe", "area", "con", "wht"])
def test_nonexistent_arrays(cube_models, array):
    """Test for non-existent arrays"""
    for cube, model in cube_models:
        with pytest.raises(AttributeError):
            model.getarray_noinit(array)


def test_run_coron3(tmp_cwd, coron3_input):
    steps = {
        # make sure outlier_detection runs
        "outlier_detection": {"skip": False},
        # skip resampling because it requires a wcs
        "resample": {"skip": True},
    }

    Coron3Pipeline.call(coron3_input, steps=steps, save_results=True)

    # Check for expected intermediate files:
    # psf_stack, psf_align, and psf_sub files are always saved;
    # outlier detection output is turned on by parameters
    basename = "test_coron3"
    expected = [f"{basename}_psfstack.fits"]
    for i in range(3):
        expected.append(f"{basename}_target_{i + 1}_c1000_psfalign.fits")
        expected.append(f"{basename}_target_{i + 1}_c1000_psfsub.fits")
        expected.append(f"{basename}_target_{i + 1}_c1000_crfints.fits")

    # Check for expected final output files: 9 "i2d" image models that
    # were input to resample, from 3 x 3 target cube models
    expected += [f"{basename}_{i}_i2d.fits" for i in range(9)]

    # Check all files
    for filename in expected:
        assert (tmp_cwd / filename).exists()
        with datamodels.open(filename) as model:
            if "i2d" in filename:
                assert isinstance(model, datamodels.ImageModel)
            else:
                assert isinstance(model, datamodels.CubeModel)

            # Make sure outlier detection ran
            assert model.meta.cal_step.outlier_detection == "COMPLETE"


def test_run_coron3_resample(monkeypatch, tmp_cwd, coron3_input):
    pipeline = Coron3Pipeline()
    # Skip irrelevant steps
    pipeline.outlier_detection.skip = True
    pipeline.stack_refs.skip = True
    pipeline.align_refs.skip = True
    pipeline.klip.skip = True

    # Mock resample
    monkeypatch.setattr(pipeline, "resample", MockResample())
    pipeline.run(coron3_input)

    # Check for expected final files: should be 1 i2d file
    filename = "test_coron3_i2d.fits"
    assert (tmp_cwd / filename).exists()
    with datamodels.open(filename) as model:
        assert isinstance(model, datamodels.ImageModel)

        # make sure association information is appended
        assert model.meta.asn.pool_name == "test_coron3_pool.csv"
        assert model.meta.asn.table_name == Path(coron3_input).name


def test_run_coron3_no_psf(tmp_cwd, caplog, coron3_science_only):
    Coron3Pipeline.call(coron3_science_only)
    assert "No reference PSF members" in caplog.text


def test_run_coron3_no_science(tmp_cwd, caplog, coron3_psf_only):
    # The expected warning is not reachable with the "call" function, because
    # CRDS checks will fail with no science members.
    # Test with the "run" function
    pipeline = Coron3Pipeline()
    pipeline.run(coron3_psf_only)
    assert "No science target members" in caplog.text
