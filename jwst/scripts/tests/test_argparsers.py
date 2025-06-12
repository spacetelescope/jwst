import pytest
from jwst.scripts.set_velocity_aberration import parse_args


def test_argparse_set_velocity_aberration():
    args = parse_args(["--force_level1bmodel", "0", "filename0", "filename1", "filename2"])
    assert args.filename == ["filename0", "filename1", "filename2"]
    assert isinstance(args.force_level1bmodel, int)
    assert not bool(args.force_level1bmodel)

    args2 = parse_args(["filename0", "filename1", "filename2", "--force_level1bmodel", "1"])
    assert args.filename == args2.filename

    args3 = parse_args(["filename0", "filename1", "filename2"])
    assert args3.force_level1bmodel is True


def test_argparse_set_velocity_aberration_bad_input(capsys):
    with pytest.raises(SystemExit):
        parse_args(["--force_level1bmodel", "1"])

    _, err = capsys.readouterr()
    assert err[:6] == "usage:"
