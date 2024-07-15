from glob import glob

from jwst.exp_to_source.tests import helpers
from jwst.exp_to_source.main import Main


def test_help(capsys):
    Main(["-h"])
    out, err = capsys.readouterr()
    assert out.startswith("usage:")


def test_default_run(tmp_path, capsys):
    args = ["-o", str(tmp_path)]
    args.extend(helpers.INPUT_FILES)
    result = Main(args)
    assert len(result.sources) == 5
    files = glob(str(tmp_path / "*.fits"))
    assert len(files) == 5


def test_dry_run(tmp_cwd, capsys):
    no_files = glob("*.fits")
    assert len(no_files) == 0

    args = ["--dry-run"]
    args.extend(helpers.INPUT_FILES)
    result = Main(args)
    assert len(result.sources) == 5
    no_files = glob("*.fits")
    assert len(no_files) == 0
