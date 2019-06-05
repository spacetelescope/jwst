if (utils.scm_checkout()) return

python_version = '3.6'

env_vars = [
    "CRDS_SERVER_URL=https://jwst-crds.stsci.edu",
    "CRDS_PATH=./crds_cache",
]

// pip related setup
def pip_index = "https://bytesalad.stsci.edu/artifactory/api/pypi/datb-pypi-virtual/simple"
def pip_install_args = "--index-url ${pip_index} --progress-bar=off"


// Generate distributions build
bc0 = new BuildConfig()
bc0.nodetype = 'linux'
bc0.name = 'wheel-sdist'
bc0.conda_ver = '4.6.14'
bc0.conda_packages = [
    "python=${python_version}",
]
bc0.build_cmds = [
    "pip install numpy",
    "pip wheel .",
    "python setup.py sdist",
]

// Generate pip build/test with released upstream dependencies
bc1 = utils.copy(bc0)
bc1.name = "stable-deps"
bc1.env_vars = env_vars
bc1.build_cmds = [
    "pip install -e .[test]",
]
bc1.test_cmds = [
    "pytest -r sx --junitxml=results.xml"
]

// Generate pip build/test with dev upstream dependencies
bc2 = utils.copy(bc1)
bc2.name = "dev-deps"
bc2.build_cmds = [
    "pip install -r requirements-dev.txt -e .[test]",
]

// Generate conda-free build with python 3.7
bc3 = new BuildConfig()
bc3.nodetype = 'python3.7'
bc3.name = 'python3.7'
bc3.build_cmds = [
    "pip install -e .[test]",
]
bc3.test_cmds = [
    "pytest -r sx --junitxml=results.xml"
]

utils.run([bc0, bc1, bc2, bc3])
