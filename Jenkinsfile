if (utils.scm_checkout()) return

python_version = '3.6'

env_vars = [
    "CRDS_SERVER_URL=https://jwst-crds.stsci.edu",
    "CRDS_PATH=./crds_cache",
    "CRDS_CONTEXT=jwst_0502.pmap",
]

// Pip related setup
def pip_index = "https://bytesalad.stsci.edu/artifactory/api/pypi/datb-pypi-virtual/simple"
def pip_install_args = "--index-url ${pip_index} --progress-bar=off"


// Generate distributions build
bc0 = new BuildConfig()
bc0.nodetype = 'linux'
bc0.name = 'wheel-sdist'
bc0.conda_ver = '4.6.8'
bc0.conda_packages = ["python=${python_version}"]
bc0.build_cmds = [
    "pip install ${pip_install_args} numpy",
    "pip wheel ${pip_install_args} .",
    "python setup.py sdist",
]

// Generate pip build/test with released upstream dependencies
bc1 = utils.copy(bc0)
bc1.name = "stable-deps"
bc1.env_vars = env_vars
bc1.build_cmds = [
    "pip install ${pip_install_args} numpy",
    "pip install ${pip_install_args} -e .[test]",
    "python setup.py develop",
]
bc1.test_cmds = ["pytest -r sx --basetemp=test_results --junitxml=results.xml"]

// Generate conda build/test with astroconda-dev dependencies
bc2 = utils.copy(bc0)
bc2.name = "astroconda-dev"
bc2.env_vars = env_vars
bc2.conda_channels = [
    "http://ssb.stsci.edu/astroconda-dev"
]
bc2.conda_packages = [
    "python=${python_version}",
    "numpy",
    "nomkl",
    "jwst",
]
bc2.build_cmds = [
    "python setup.py develop",
    "pip install -e .[test]",
]
bc2.test_cmds = ["pytest -r sx --basetemp=test_results --junitxml=results.xml"]

utils.run([bc0, bc1, bc2])
