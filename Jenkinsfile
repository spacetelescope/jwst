if (utils.scm_checkout()) return

python_ver = '3.6'
pip_numpy_ver = '1.16'
conda_numpy_conda_ver = '1.15'

def test_env = [
    "CRDS_SERVER_URL=https://jwst-crds.stsci.edu",
    "CRDS_PATH=./crds_cache",
]

// Define conda packages needed for build bc2.  Some come from astroconda-dev
def conda_packages = [
    "asdf",
    "astropy",
    "crds",
    "drizzle",
    "flake8",
    "gwcs",
    "jsonschema",
    "jplephem",
    "matplotlib",
    "photutils",
    "scipy",
    "spherical-geometry",
    "stsci.image",
    "stsci.imagestats",
    "stsci.stimage",
    "verhawk",
    "pytest",
]

// Pip related setup
def pip_index = "https://bytesalad.stsci.edu/artifactory/api/pypi/datb-pypi-virtual/simple"
def pip_install_args = "--index-url ${pip_index} --progress-bar=off"

// Generate distributions
bc0 = new BuildConfig()
bc0.nodetype = 'linux'
bc0.name = 'dist'
bc0.conda_packages = ["python=${python_ver}"]
bc0.build_cmds = [
    "pip install ${pip_install_args} numpy==${matrix_numpy[0]}",
    "pip wheel ${pip_install_args} .",
    "python setup.py sdist",
]

// Generate pip build/test with released upstream dependencies
bc1 = utils.copy(bc0)
bc1.name = "released dependencies"
bc1.env_vars = test_env
bc1.build_cmds = [
    "pip install ${pip_install_args} numpy~=${pip_numpy_ver}",
    "pip install ${pip_install_args} -e .[test]",
    "python setup.py develop",
]
bc1.test_cmds = ["pytest -r sx --basetemp=test_results --junitxml=results.xml"]

// Generate conda build/test with astroconda-dev dependencies
bc2 = utils.copy(bc0)
bc2.name = "astroconda-dev"
bc2.env_vars = test_env
bc2.conda_channels = ['http://ssb.stsci.edu/astroconda-dev']
bc2.conda_packages += conda_packages + ["numpy=${conda_numpy_ver}"]
bc2.build_cmds = [
    "pip install ${pip_install_args} -e .[test]",
    "python setup.py develop",
]
bc2.test_cmds = ["pytest -r sx --basetemp=test_results --junitxml=results.xml"]

utils.run([bc0, bc1, bc2])
