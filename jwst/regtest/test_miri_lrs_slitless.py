import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step
from jwst.associations.asn_from_list import asn_from_list

DATASET_ID = "jw00623026001_03106_00005_mirimage"
PRODUCT_NAME = "jw00623-a3001_t001_miri_p750l-slitlessprism"
PROGRAM = "00623"


@pytest.fixture(scope="module")
def run_tso1_pipeline(jail, rtdata_module):
    """Run the calwebb_tso1 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module
    collect_pipeline_cfgs("config")
    rtdata.get_data(f"miri/lrs/{DATASET_ID}_uncal.fits")

    args = [
        "config/calwebb_tso1.cfg",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.rscd.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.jump.rejection_threshold=150",
        "--steps.jump.save_results=True"
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso_spec2_pipeline(run_tso1_pipeline, jail, rtdata_module):
    """Run the calwebb_tso-spec2 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module

    rtdata.input = f"{DATASET_ID}_rateints.fits"

    collect_pipeline_cfgs("config")

    args = [
        "config/calwebb_tso-spec2.cfg",
        rtdata.input,
        "--steps.flat_field.save_results=true",
        "--steps.srctype.save_results=true",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def generate_tso3_asn():
    """Generate an association file that references the output of run_spec2_pipeline."""
    asn = asn_from_list([f"{DATASET_ID}_calints.fits"], product_name=PRODUCT_NAME)
    asn.data["program"] = PROGRAM
    asn.data["asn_type"] = "tso3"
    asn.sequence = 1

    name, serialized = asn.dump(format="json")
    with open(name, "w") as f:
        f.write(serialized)

    return name, asn["asn_id"]


@pytest.fixture(scope="module")
def run_tso3_pipeline(run_tso_spec2_pipeline, generate_tso3_asn):
    """Run the calwebb_tso3 pipeline on the output of run_spec2_pipeline."""
    asn_filename, _ = generate_tso3_asn

    args = [
        "config/calwebb_tso3.cfg",
        asn_filename,
        "--steps.outlier_detection.save_results=true",
    ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("step_suffix", ['dq_init', 'saturation', 'refpix',
    'rscd', 'linearity', 'dark_current', 'jump', 'rate', 'rateints'])
def test_miri_lrs_slitless_tso1(run_tso1_pipeline, rtdata_module, fitsdiff_default_kwargs, step_suffix):
    """
    Regression test of tso1 pipeline performed on MIRI LRS slitless TSO data.
    """
    rtdata = rtdata_module
    output_filename = f"{DATASET_ID}_{step_suffix}.fits"
    rtdata.output = output_filename

    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso1/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("step_suffix", ["flat_field", "srctype", "calints", "x1dints"])
def test_miri_lrs_slitless_tso_spec2(run_tso_spec2_pipeline, rtdata_module, fitsdiff_default_kwargs,
    step_suffix):
    """Compare the output of a MIRI LRS slitless calwebb_tso-spec2 pipeline."""
    rtdata = rtdata_module

    output_filename = f"{DATASET_ID}_{step_suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso_spec2/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("step_suffix, filename_extension, is_product",
    [
        ("outlier_detection", "fits", False),
        ("crfints", "fits", False),
        ("x1dints", "fits", True),
        ("whtlt", "ecsv", True),
    ],
)
def test_miri_lrs_slitless_tso3(run_tso3_pipeline, generate_tso3_asn, rtdata_module,
    fitsdiff_default_kwargs, diff_astropy_tables, step_suffix, filename_extension,
    is_product):
    """Compare the output of a MIRI LRS slitless calwebb_tso3 step."""
    rtdata = rtdata_module
    _, asn_id = generate_tso3_asn

    if is_product:
        output_filename = f"{PRODUCT_NAME}_{step_suffix}.{filename_extension}"
    else:
        output_filename = f"{DATASET_ID}_{asn_id}_{step_suffix}.{filename_extension}"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso3/{output_filename}")

    if filename_extension == "fits":
        diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
        assert diff.identical, diff.report()
    else:
        diff = diff_astropy_tables(rtdata.output, rtdata.truth)
        assert len(diff) == 0, "\n".join(diff)
