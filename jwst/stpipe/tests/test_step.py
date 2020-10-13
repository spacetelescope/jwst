import os
from os.path import (
    abspath,
    dirname,
    join,
)

import pytest
import jwst
from jwst import datamodels
from jwst.extern.configobj.configobj import ConfigObj
from jwst.refpix import RefPixStep
from jwst.stpipe import Step, crds_client
from jwst.stpipe import cmdline
from jwst.stpipe.config_parser import ValidationError

from .steps import EmptyPipeline, MakeListPipeline, MakeListStep, ProperPipeline
from .util import t_path

from crds.core.exceptions import CrdsLookupError


ParsModelWithPar3 = datamodels.StepParsModel(t_path(join('steps','jwst_generic_pars-makeliststep_0002.asdf')))
ParsModelWithPar3.parameters.instance.update({'par3': False})

REFPIXSTEP_CRDS_MIRI_PARS = {
    'class': jwst.refpix.refpix_step.RefPixStep,
    'name': 'refpix',
    'odd_even_columns': False,
    'odd_even_rows': False,
    'side_gain': 10.0,
    'side_smoothing_length': 21,
    'use_side_ref_pixels': False
}

CRDS_ERROR_STRING = 'PARS-WITHDEFAULTSSTEP: No parameters found'

@pytest.fixture(scope='module')
def data_path():
    """Provide a test data model"""
    data_path = t_path(join('data', 'miri_data.fits'))
    return data_path


@pytest.mark.parametrize(
    'arg, env_set, expected_fn', [
        ('--disable-crds-steppars', None,   lambda stream: not CRDS_ERROR_STRING in stream),
        ('--verbose',               None,   lambda stream: CRDS_ERROR_STRING in stream),
        ('--verbose',               'true', lambda stream: not CRDS_ERROR_STRING in stream),
        ('--verbose',               'True', lambda stream: not CRDS_ERROR_STRING in stream),
        ('--verbose',               't',    lambda stream: not CRDS_ERROR_STRING in stream),
    ]
)
def test_disable_crds_steppars_cmdline(capsys, data_path, arg, env_set, expected_fn):
    """Test setting of disable_crds_steppars"""
    if env_set:
        os.environ['STPIPE_DISABLE_CRDS_STEPPARS'] = env_set

    try:
        step, step_class, positional, debug_on_exception = cmdline.just_the_step_from_cmdline(
            ['jwst.stpipe.tests.steps.WithDefaultsStep', data_path, arg]
        )
    finally:
        os.environ.pop('STPIPE_DISABLE_CRDS_STEPPARS', None)

    captured = capsys.readouterr()
    assert expected_fn(captured.err)


@pytest.mark.xfail(
    reason='Need to make regression with actual CRDS data',
    run=False,
)
def test_parameters_from_crds():
    """Test retrieval of parameters from CRDS"""
    step_class = REFPIXSTEP_CRDS_MIRI_PARS['class']
    data = datamodels.open(t_path(join('data', 'miri_data.fits')))
    pars = step_class.get_config_from_reference(data)
    assert pars == REFPIXSTEP_CRDS_MIRI_PARS


def test_parameters_from_crds_fail():
    """Test retrieval of parameters from CRDS"""
    with datamodels.open(t_path(join('data', 'miri_data.fits'))) as data:
        data.meta.instrument.name = 'NIRSPEC'
        pars = RefPixStep.get_config_from_reference(data)
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
    assert step.get_pars_model().meta.reftype == 'pars-' + expected.lower()


def test_noproperty_pars():
    """Ensure that property parameters are excluded"""
    pars = Step.get_pars()
    assert pars['input_dir'] is None


def test_saving_pars(tmpdir):
    """Save the step parameters from the commandline"""
    cfg_path = t_path(join('steps', 'jwst_generic_pars-makeliststep_0002.asdf'))
    saved_path = tmpdir.join('savepars.asdf')
    step = Step.from_cmdline([
        cfg_path,
        '--save-parameters',
        str(saved_path)
    ])
    assert saved_path.check()

    with datamodels.StepParsModel(str(saved_path)) as saved:
        assert saved.parameters == ParsModelWithPar3.parameters

    step.closeout()


@pytest.mark.parametrize(
    'step_obj, full_spec, expected',
    [
        # #############################################
        # Test `get_pars_model` with `full_spec = True`
        # #############################################
        #
        # Step Class with mix of required and optional parameters
        (MakeListStep, True,
         {'class': 'jwst.stpipe.tests.steps.MakeListStep',
          'name': 'MakeListStep',
          'pre_hooks': [],
          'post_hooks': [],
          'output_file': None,
          'output_dir': None,
          'output_ext': '.fits',
          'output_use_model': False,
          'output_use_index': True,
          'save_results': False,
          'skip': False,
          'suffix': None,
          'search_output_file': True,
          'input_dir': None,
          'par3': False,
          'par1': 'float() # Control the frobulization',
          'par2': 'string() # Reticulate the splines'}
        ),
        # Instance with parameters set
        (MakeListStep(par1=0., par2='from args'), True,
         {'class': 'jwst.stpipe.tests.steps.MakeListStep',
          'name': 'MakeListStep',
          'pre_hooks': [],
          'post_hooks': [],
          'output_file': None,
          'output_dir': None,
          'output_ext': '.fits',
          'output_use_model': False,
          'output_use_index': True,
          'save_results': False,
          'skip': False,
          'suffix': None,
          'search_output_file': True,
          'input_dir': '',
          'par1': 0.0,
          'par2': 'from args',
          'par3': False}
        ),
        # ##############################################
        # Test `get_pars_model` with `full_spec = False`
        # ##############################################
        #
        # Step Class with mix of required and optional parameters
        (MakeListStep, False,
         {'class': 'jwst.stpipe.tests.steps.MakeListStep',
          'name': 'MakeListStep',
          'par3': False,
          'par1': 'float() # Control the frobulization',
          'par2': 'string() # Reticulate the splines'}
        ),
        # Instance with parameters set
        (MakeListStep(par1=0., par2='from args'), False,
         {'class': 'jwst.stpipe.tests.steps.MakeListStep',
          'name': 'MakeListStep',
          'par1': 0.0,
          'par2': 'from args',
          'par3': False}
        ),
    ]
)
def test_getpars_model(step_obj, full_spec, expected):
    """Test retreiving of configuration parameters"""
    assert step_obj.get_pars_model(full_spec=full_spec).parameters.instance == expected


@pytest.mark.parametrize(
    'step_obj, full_spec, expected',
    [
        # #######################################
        # Test `get_pars` with `full_spec = True`
        # #######################################
        #
        # Step Class with mix of required and optional parameters
        (MakeListStep, True, {
            'pre_hooks': [],
            'post_hooks': [],
            'output_file': None,
            'output_dir': None,
            'output_ext': '.fits',
            'output_use_model': False,
            'output_use_index': True,
            'save_results': False,
            'skip': False,
            'suffix': None,
            'search_output_file': True,
            'input_dir': None,
            'par3': False,
            'par1': 'float() # Control the frobulization',
            'par2': 'string() # Reticulate the splines'
        }),
        # Instance with parameters set
        (MakeListStep(par1=0., par2='from args'), True, {
            'pre_hooks': [],
            'post_hooks': [],
            'output_file': None,
            'output_dir': None,
            'output_ext': '.fits',
            'output_use_model': False,
            'output_use_index': True,
            'save_results': False,
            'skip': False,
            'suffix': None,
            'search_output_file': True,
            'input_dir': '',
            'par1': 0.0,
            'par2': 'from args',
            'par3': False
        }),
        #
        # Pipeline Class with sub-step with mix of required and optional parameters
        #
        (MakeListPipeline, True, {
            'pre_hooks': [],
            'post_hooks': [],
            'output_file': None,
            'output_dir': None,
            'output_ext': '.fits',
            'output_use_model': False,
            'output_use_index': True,
            'save_results': False,
            'skip': False,
            'suffix': None,
            'search_output_file': True,
            'input_dir': None,
            'par1': 'Name the atomizer',
            'steps': {
                'make_list': {'pre_hooks': [],
                              'post_hooks': [],
                              'output_file': None,
                              'output_dir': None,
                              'output_ext': '.fits',
                              'output_use_model': False,
                              'output_use_index': True,
                              'save_results': False,
                              'skip': False,
                              'suffix': None,
                              'search_output_file': True,
                              'input_dir': None,
                              'par3': False,
                              'par1': 'float() # Control the frobulization',
                              'par2': 'string() # Reticulate the splines'
                }
            }
        }),
        #
        # Pipeline Instance with all parameters set
        #
        (MakeListPipeline(
            par1='Instantiated', steps={'make_list': {'par1': 0., 'par2': 'sub-instantiated'}}
        ), True, {
            'pre_hooks': [],
            'post_hooks': [],
            'output_file': None,
            'output_dir': None,
            'output_ext': '.fits',
            'output_use_model': False,
            'output_use_index': True,
            'save_results': False,
            'skip': False,
            'suffix': None,
            'search_output_file': True,
            'input_dir': '',
            'par1': 'Instantiated',
            'steps': {
                'make_list': {
                    'pre_hooks': [],
                    'post_hooks': [],
                    'output_file': None,
                    'output_dir': None,
                    'output_ext': '.fits',
                    'output_use_model': False,
                    'output_use_index': True,
                    'save_results': False,
                    'skip': False,
                    'suffix': None,
                    'search_output_file': True,
                    'input_dir': '',
                    'par1': 0.0,
                    'par2': 'sub-instantiated',
                    'par3': False
                }
            }
        }),
        #
        # Pipeline class without any sub-steps
        #
        (EmptyPipeline, True, {
            'pre_hooks': [],
            'post_hooks': [],
            'output_file': None,
            'output_dir': None,
            'output_ext': '.fits',
            'output_use_model': False,
            'output_use_index': True,
            'save_results': False,
            'skip': False,
            'suffix': None,
            'search_output_file': True,
            'input_dir': None,
            'par1': 'Name the atomizer',
            'steps': {}
        }),
        #
        # Pipeline instance without any sub-steps
        #
        (EmptyPipeline(par1='Instantiated'), True, {
            'pre_hooks': [],
            'post_hooks': [],
            'output_file': None,
            'output_dir': None,
            'output_ext': '.fits',
            'output_use_model': False,
            'output_use_index': True,
            'save_results': False,
            'skip': False,
            'suffix': None,
            'search_output_file': True,
            'input_dir': '',
            'par1': 'Instantiated',
            'steps': {}
        }),
        # ######################################
        # Test `get_pars` with `full_spec=False`
        # ######################################
        #
        # Step Class with mix of required and optional parameters
        #
        (MakeListStep, False, {
            'par3': False,
            'par1': 'float() # Control the frobulization',
            'par2': 'string() # Reticulate the splines'
        }),
        #
        # Instance with parameters set
        #
        (MakeListStep(par1=0., par2='from args'), False, {
            'par1': 0.0,
            'par2': 'from args',
            'par3': False
        }),
        #
        # Pipeline Class with sub-step with mix of required and optional parameters
        #
        (MakeListPipeline, False, {
            'par1': 'Name the atomizer',
            'steps': {
                'make_list': {
                    'par3': False,
                    'par1': 'float() # Control the frobulization',
                    'par2': 'string() # Reticulate the splines'
                }
            }
        }),
        #
        # Pipeline Instance with all parameters set
        #
        (MakeListPipeline(
            par1='Instantiated', steps={'make_list': {'par1': 0., 'par2': 'sub-instantiated'}}
        ), False, {
            'par1': 'Instantiated',
            'steps': {
                'make_list': {
                    'par1': 0.0,
                    'par2': 'sub-instantiated',
                    'par3': False
                }
            }
        }),
        #
        # Pipeline class without any sub-steps
        #
        (EmptyPipeline, False, {
            'par1': 'Name the atomizer',
            'steps': {}
        }),
        #
        # Pipeline instance without any sub-steps
        #
        (EmptyPipeline(par1='Instantiated'), False, {
            'par1': 'Instantiated',
            'steps': {}
        }),
    ]
)
def test_getpars(step_obj, full_spec, expected):
    """Test retreiving of configuration parameters"""
    assert step_obj.get_pars(full_spec=full_spec) == expected


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


@pytest.mark.parametrize(
    "command_line_pars, command_line_config_pars, reference_pars, expected_pars",
    [
        # If nothing else is present, we should use the spec defaults
        (
            None,
            None,
            None,
            {"par1": "default par1 value", "par2": "default par2 value",
             "par3": "default par3 value", "par4": "default par4 value"}
        ),
        # Reference file pars > spec defaults
        (
            None,
            None,
            {"par1": "reference par1 value", "par2": "reference par2 value",
             "par3": "reference par3 value", "par4": "reference par4 value"},
            {"par1": "reference par1 value", "par2": "reference par2 value",
             "par3": "reference par3 value", "par4": "reference par4 value"}
        ),
        # Command line config pars > reference pars
        (
            None,
            {"par1": "config par1 value", "par2": "config par2 value",
             "par3": "config par3 value", "par4": "config par4 value"},
            {"par1": "reference par1 value", "par2": "reference par2 value",
             "par3": "reference par3 value", "par4": "reference par4 value"},
            {"par1": "config par1 value", "par2": "config par2 value",
             "par3": "config par3 value", "par4": "config par4 value"}
        ),
        # Command line override pars > all other pars
        (
            {"par1": "override par1 value", "par2": "override par2 value",
             "par3": "override par3 value", "par4": "override par4 value"},
            {"par1": "config par1 value", "par2": "config par2 value",
             "par3": "config par3 value", "par4": "config par4 value"},
            {"par1": "reference par1 value", "par2": "reference par2 value",
             "par3": "reference par3 value", "par4": "reference par4 value"},
            {"par1": "override par1 value", "par2": "override par2 value",
             "par3": "override par3 value", "par4": "override par4 value"}
        ),
        # Test complex merging of parameters (one parameter per source)
        (
            {"par1": "override par1 value"},
            {"par2": "config par2 value"},
            {"par3": "reference par3 value"},
            {"par1": "override par1 value", "par2": "config par2 value",
             "par3": "reference par3 value", "par4": "default par4 value"}
        ),
        # Test complex merging of parameters (more parameters specified on the lower-precedence sources)
        (
            {"par1": "override par1 value"},
            {"par1": "config par1 value", "par2": "config par2 value"},
            {"par1": "reference par1 value", "par2": "reference par2 value", "par3": "reference par3 value"},
            {"par1": "override par1 value", "par2": "config par2 value",
             "par3": "reference par3 value", "par4": "default par4 value"}
        ),
    ]
)
def test_step_from_commandline_par_precedence(command_line_pars, command_line_config_pars,
                                              reference_pars, expected_pars, tmp_path, monkeypatch):
    args = []

    class_name = "jwst.stpipe.tests.steps.WithDefaultsStep"
    config_name = "WithDefaultsStep"
    reference_type = f"pars-{config_name.lower()}"
    input_path = join(dirname(__file__), "data", "science.fits")

    if command_line_config_pars:
        command_line_config_path = tmp_path/"with_defaults_step.cfg"
        config = ConfigObj(str(command_line_config_path))
        config["class"] = class_name
        config["name"] = config_name
        for key, value in command_line_config_pars.items():
            config[key] = value
        config.write()
        args.append(str(command_line_config_path.absolute()))
    else:
        args.append(class_name)

    args.append(input_path)

    if command_line_pars:
        for key, value in command_line_pars.items():
            args.append(f"--{key}={value}")

    reference_file_map = {}
    if reference_pars:
        reference_path = tmp_path/f"{reference_type}.asdf"
        parameters = {
            "class": class_name,
            "name": config_name,
        }
        parameters.update(reference_pars)
        model = datamodels.StepParsModel({"parameters": parameters})
        model.save(reference_path)

        reference_file_map[reference_type] = str(reference_path)

    def mock_get_reference_file(dataset, reference_file_type, observatory=None, asn_exptypes=None):
        if reference_file_type in reference_file_map:
            return reference_file_map[reference_file_type]
        else:
            raise CrdsLookupError(f"Error determining best reference for '{reference_file_type}'  = \
  Unknown reference type '{reference_file_type}'")
    monkeypatch.setattr(crds_client, "get_reference_file", mock_get_reference_file)

    step = Step.from_cmdline(args)

    for key, value in expected_pars.items():
        assert getattr(step, key) == value


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


def test_call_with_config(caplog, _jail):
    """Test call using a config file with substeps

    In particular, from JP-1482, there was a case where a substep parameter
    was not being overridden. Test for that case.
    """
    cfg = t_path(join('data', 'proper_pipeline.asdf'))
    model = t_path(join('data', 'flat.fits'))

    ProperPipeline.call(model, config_file=cfg)

    assert "'par1': 'newpar1'" in caplog.text
