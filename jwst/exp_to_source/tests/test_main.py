from glob import glob

import pytest

from jwst.exp_to_source.main import Main


@pytest.fixture
def save_input_files(mock_input, tmp_cwd):
    for model in mock_input:
        model.save(tmp_cwd / model.meta.filename)
    yield [model.meta.filename for model in mock_input]
    for model in mock_input:
        model.close()


def test_help(capsys):
    Main(["-h"])
    out, _err = capsys.readouterr()
    assert out.startswith("usage:")


def test_default_run(save_input_files, tmp_cwd, capsys):
    output_dir = tmp_cwd / "output"
    output_dir.mkdir(exist_ok=True)
    args = ["-o", str(output_dir)]
    args.extend(save_input_files)
    result = Main(args)
    assert len(result.sources) == 5
    files = glob(str(output_dir / "*.fits"))
    assert len(files) == 5


def test_dry_run(save_input_files, tmp_cwd, capsys):
    output_dir = tmp_cwd / "output"
    output_dir.mkdir(exist_ok=True)
    no_files = glob(str(output_dir / "*.fits"))
    assert len(no_files) == 0

    args = ["--dry-run"]
    args.extend(save_input_files)
    result = Main(args)
    assert len(result.sources) == 5
    no_files = glob(str(output_dir / "*.fits"))
    assert len(no_files) == 0
