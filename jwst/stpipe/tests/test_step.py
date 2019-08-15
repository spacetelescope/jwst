from os.path import (
    abspath,
    dirname,
    join,
)

import pytest

from jwst import datamodels
from jwst.refpix import RefPixStep
from jwst.stpipe import Step
from jwst.stpipe.config_parser import ValidationError
from jwst.stpipe import crds_client

from .steps import MakeListStep
from .util import t_path

ParsModelWithPar3 = datamodels.StepParsModel(t_path(join('steps','jwst_generic_pars-makeliststep_0002.asdf')))
ParsModelWithPar3.parameters.instance.update({'par3': False})

REFPIXSTEP_CRDS_MIRI_PARS = {
    'class': 'jwst.refpix.refpix_step.RefPixStep',
    'name': 'refpix',
    'odd_even_columns': False,
    'odd_even_rows': False,
    'side_gain': 10.0,
    'side_smoothing_length': 21,
    'use_side_ref_pixels': False
}

def test_parameters_from_crds():
    """Test retrieval of parameters from CRDS"""
    data = datamodels.open(t_path(join('data', 'miri_data.fits')))
    pars = crds_client.get_parameters_from_reference(RefPixStep, data)
    assert pars == REFPIXSTEP_CRDS_MIRI_PARS


def test_parameters_from_crds_fail():
    """Test retrieval of parameters from CRDS"""
    data = datamodels.open(t_path(join('data', 'miri_data.fits')))
    data.meta.instrument.name = 'NIRSPEC'
    pars = crds_client.get_parameters_from_reference(RefPixStep, data)
    assert not len(pars)


@pytest.mark.parametrize(
    'cfg_file, expected',
    [
        ('local_class.cfg', 'TestLocallyInstalledStepClass'),
        ('jwst_generic_pars-makeliststep_0002.asdf', 'jwst_generic_pars-makeliststep_0002'),
    ]
)
def test_reftype(cfg_file, expected):
    """Test that reftype is produce as expected"""
    step = Step.from_config_file(t_path(join('steps', cfg_file)))
    assert step.pars_model.meta.reftype == 'pars-' + expected.lower()


def test_noproperty_pars():
    """Ensure that property parameters are excluded"""
    pars = Step.pars
    assert pars['input_dir'] is None


def test_saving_pars(tmpdir):
    """Save the step parameters from the commandline"""
    cfg_path = t_path(join('steps', 'jwst_generic_pars-makeliststep_0002.asdf'))
    saved_path = tmpdir.join('savepars.asdf')
    Step.from_cmdline([
        cfg_path,
        '--save-parameters',
        str(saved_path)
    ])
    assert saved_path.check()
    saved = datamodels.StepParsModel(str(saved_path))
    assert saved.parameters == ParsModelWithPar3.parameters


@pytest.mark.parametrize(
    'step_obj, expected',
    [
        (MakeListStep,
         datamodels.StepParsModel(
             {'parameters': {
                 'par1': 'float() # Control the frobulization',
                 'par2': 'string() # Reticulate the splines',
                 'par3': False,
                 'class': 'jwst.stpipe.tests.steps.MakeListStep',
                 'name': 'MakeListStep',
             }}
         )
        ),
        (MakeListStep(par1=0., par2='from args'),
         datamodels.StepParsModel({'parameters': {
             'par1': 0., 'par2': 'from args', 'par3': False,
             'class': 'jwst.stpipe.tests.steps.MakeListStep',
             'name': 'MakeListStep'
         }})
        ),
        (Step.from_config_file(t_path(join('steps', 'jwst_generic_pars-makeliststep_0002.asdf'))),
         ParsModelWithPar3
        )
    ]
)
def test_getpars_model(step_obj, expected):
    """Test retreiving of configuration parameters"""
    assert step_obj.pars_model.parameters == expected.parameters


@pytest.mark.parametrize(
    'step_obj, expected',
    [
        (MakeListStep, {
            'par1': 'float() # Control the frobulization',
            'par2': 'string() # Reticulate the splines',
            'par3': False
        }),
        (MakeListStep(par1=0., par2='from args'), {'par1': 0., 'par2': 'from args', 'par3': False}),
    ]
)
def test_getpars(step_obj, expected):
    """Test retreiving of configuration parameters"""
    assert step_obj.pars == expected


def test_hook():
    """Test the running of hooks"""
    step_fn = join(dirname(__file__), 'steps', 'stepwithmodel_hook.cfg')
    step = Step.from_config_file(step_fn)

    model = datamodels.ImageModel()
    result = step.run(model)

    assert result.pre_hook_run
    assert step.pre_hook_run
    assert result.post_hook_run
    assert step.post_hook_run


def test_hook_with_return():
    """Test the running of hooks"""
    step_fn = join(dirname(__file__), 'steps', 'stepwithmodel_hookreturn.cfg')
    step = Step.from_config_file(step_fn)

    model = datamodels.ImageModel()
    result = step.run(model)

    assert result == 'PostHookWithReturnStep executed'
    assert step.pre_hook_run
    assert step.post_hook_run


def test_step():
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
    args = [
        'jwst.stpipe.tests.steps.AnotherDummyStep',
        '--par1=58', '--par2=hij klm'
        ]

    step = Step.from_cmdline(args)

    assert step.par1 == 58.
    assert step.par2 == 'hij klm'
    assert step.par3 is False

    step.run(1, 2)


def test_step_from_commandline_invalid():
    args = [
            '__foo__'
        ]

    with pytest.raises(ValueError):
        Step.from_cmdline(args)


def test_step_from_commandline_invalid2():
    args = [
        '__foo__.__bar__'
        ]
    with pytest.raises(ValueError):
        Step.from_cmdline(args)


def test_step_from_commandline_invalid3():
    args = [
        'sys.foo'
        ]
    with pytest.raises(ValueError):
        Step.from_cmdline(args)


def test_step_from_commandline_invalid4():
    args = [
        'sys.argv'
        ]
    with pytest.raises(ValueError):
        Step.from_cmdline(args)


def test_step_with_local_class():
    step_fn = join(dirname(__file__), 'steps', 'local_class.cfg')
    step = Step.from_config_file(step_fn)

    step.run(datamodels.ImageModel((2, 2)))


def test_extra_parameter():
    from .steps import AnotherDummyStep
    with pytest.raises(ValidationError):
        AnotherDummyStep("SomeOtherStepOriginal", par5='foo')


def test_crds_override():
    from .steps import AnotherDummyStep

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


def test_search_attr():
    from .steps import SavePipeline

    value = '/tmp'
    pipeline = SavePipeline('afile.fits', output_dir=value)

    assert pipeline.search_attr('output_dir') == value
    assert pipeline.stepwithmodel.search_attr('output_dir') == value
    assert pipeline.search_attr('junk') is None
    assert pipeline.stepwithmodel.search_attr('junk') is None


def test_print_configspec():
    step = Step()
    step.print_configspec()
