from __future__ import absolute_import, division, print_function

from os.path import dirname, join, abspath, isfile
import shutil
import tempfile

#from nose.tools import raises
import pytest
import numpy as np

from ..config_parser import ValidationError
from ..step import Step


def test_step():
    from .. import Step

    step_fn = join(dirname(__file__), 'steps', 'some_other_step.cfg')
    step = Step.from_config_file(step_fn)

    from .steps import AnotherDummyStep

    assert isinstance(step, AnotherDummyStep)
    assert step.name == 'SomeOtherStepOriginal'
    assert step.par2 == 'abc def'

    step.run(1, 2)


def test_step_from_python():
    from .steps import AnotherDummyStep

    step = AnotherDummyStep("SomeOtherStepOriginal", par1=42.0, par2="abc def")

    assert step.par1 == 42.0
    assert step.par2 == 'abc def'
    assert step.par3 is False

    result = step.run(1, 2)

    assert result == 3


def test_step_from_python_simple():
    from .steps import AnotherDummyStep

    result = AnotherDummyStep.call(1, 2, par1=42.0, par2="abc def")

    assert result == 3


def test_step_from_python_simple2():
    from .steps import AnotherDummyStep

    step_fn = join(dirname(__file__), 'steps', 'some_other_step.cfg')

    result = AnotherDummyStep.call(1, 2, config_file=step_fn)

    assert result == 3



def test_step_from_commandline():
    from .. import Step

    args = [
        abspath(join(dirname(__file__), 'steps', 'some_other_step.cfg')),
        '--par1=58', '--par2=hij klm'
        ]

    step = Step.from_cmdline(args)

    assert step.par1 == 58.
    assert step.par2 == 'hij klm'
    assert step.par3 is True

    step.run(1, 2)


def test_step_from_commandline_class():
    from .. import Step

    args = [
        'jwst.stpipe.tests.steps.AnotherDummyStep',
        '--par1=58', '--par2=hij klm'
        ]

    step = Step.from_cmdline(args)

    assert step.par1 == 58.
    assert step.par2 == 'hij klm'
    assert step.par3 is False

    step.run(1, 2)


#@raises(ValueError)
def test_step_from_commandline_invalid():
    from .. import Step
    args = [
            '__foo__'
        ]

    with pytest.raises(ValueError):
        step = Step.from_cmdline(args)


#@raises(ValueError)
def test_step_from_commandline_invalid2():
    from .. import Step

    args = [
        '__foo__.__bar__'
        ]
    with pytest.raises(ValueError):
        step = Step.from_cmdline(args)


#@raises(ValueError)
def test_step_from_commandline_invalid3():
    from .. import Step

    args = [
        'sys.foo'
        ]
    with pytest.raises(ValueError):
        step = Step.from_cmdline(args)


#@raises(ValueError)
def test_step_from_commandline_invalid4():
    from .. import Step

    args = [
        'sys.argv'
        ]
    with pytest.raises(ValueError):
        step = Step.from_cmdline(args)


def test_step_print_spec():
    import io
    buf = io.BytesIO()

    from .. import subproc

    subproc.SystemCall.print_configspec(buf)

    content = buf.getvalue()

    # TODO: Assert some things


def test_step_with_local_class():
    from .. import Step

    step_fn = join(dirname(__file__), 'steps', 'local_class.cfg')
    step = Step.from_config_file(step_fn)

    step.run(np.array([[0,0]]))


#@raises(ValidationError)
def test_extra_parameter():
    from .steps import AnotherDummyStep
    with pytest.raises(ValidationError):
        step = AnotherDummyStep("SomeOtherStepOriginal", par5='foo')


def test_crds_override():
    from .steps import AnotherDummyStep
    from ... import datamodels

    step = AnotherDummyStep(
        "SomeOtherStepOriginal",
        par1=42.0, par2="abc def",
        override_flat_field=join(dirname(__file__), 'data', 'flat.fits'))

    fd = step.get_reference_file(datamodels.open(), 'flat_field')
    assert fd == join(dirname(__file__), 'data', 'flat.fits')


def test_omit_ref_file():
    from .steps import OptionalRefTypeStep

    step = OptionalRefTypeStep(override_to_be_ignored_ref_type="")
    step.process()

'''
def test_save_model():
    tempdir = tempfile.mkdtemp()
    orig_filename = join(dirname(__file__), 'data', 'flat.fits')
    temp_filename = join(tempdir, 'flat_FOO.fits')
    shutil.copyfile(orig_filename, temp_filename)

    args = [
        'jwst.stpipe.tests.steps.SaveStep',
        temp_filename
    ]

    Step.from_cmdline(args)
    fname = join(tempdir, 'flat_processed.fits')
    assert isfile(fname)
'''