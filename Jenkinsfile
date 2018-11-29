if (utils.scm_checkout()) return

matrix_python = ['3.6']
matrix_numpy = ['1.15']
matrix_astropy = ['4']
matrix = []

def test_env = [
    "CRDS_SERVER_URL=https://jwst-crds.stsci.edu",
    "CRDS_PATH=./crds_cache",
]

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
    "namedlist",
    "photutils",
    "scipy",
    "spherical-geometry",
    "stsci.image",
    "stsci.imagestats",
    "stsci.stimage",
    "stsci.tools",
    "verhawk",
    "pytest"
]
def conda_packages_docs = [
    "sphinx",
    "sphinx_rtd_theme",
    "stsci_rtd_theme"
]

// Pip related setup
def pip_packages_docs = "sphinx-automodapi"
def pip_packages_tests = "requests_mock ci_watson"
def pip_index = "https://bytesalad.stsci.edu/artifactory/api/pypi/datb-pypi-virtual/simple"
def pip_install_args = "--index-url ${pip_index} --src /tmp --progress-bar=off"

// Generate distributions
dist = new BuildConfig()
dist.nodetype = 'linux'
dist.name = 'dist'
dist.build_cmds = [
    "pip install numpy==${matrix_numpy[0]}",
    "python setup.py sdist",
    "python setup.py bdist_wheel"
]
matrix += dist

/*
// Compile documentation
docs = new BuildConfig()
docs.nodetype = 'linux'
docs.name = 'docs'
docs.conda_channels = ['http://ssb.stsci.edu/astroconda-dev']
docs.conda_packages = conda_packages + conda_packages_docs + ["numpy=${matrix_numpy[0]}"] + ["python=${matrix_python[0]}"]
docs.build_cmds = [
    "pip install -q ${pip_packages_docs}",
    "python setup.py build_sphinx"
]
matrix += docs
*/

// Generate pip build and test matrix
for (python_ver in matrix_python) {
    for (numpy_ver in matrix_numpy) {
        for (astropy_ver in matrix_astropy) {
            def name = "pip_py${python_ver}np${numpy_ver}ap${astropy_ver}"
            bc = new BuildConfig()
            bc.nodetype = 'linux'
            bc.env_vars = test_env
            bc.name = name
            bc.conda_packages = ["python=${python_ver}"]
            bc.build_cmds = [
                "pip install ${pip_install_args} numpy==${matrix_numpy[0]}",
                "pip install ${pip_install_args} ${pip_packages_tests}",
                "pip install ${pip_install_args} -r requirements-dev.txt .",
            ]
            bc.test_cmds = ["pytest -r s --basetemp=test_results --junitxml=results.xml"]
            matrix += bc
        }
    }
}

/*
// Generate conda build and test matrix
for (python_ver in matrix_python) {
    for (numpy_ver in matrix_numpy) {
        for (astropy_ver in matrix_astropy) {
            def name = "conda_py${python_ver}np${numpy_ver}ap${astropy_ver}"
            bc = new BuildConfig()
            bc.nodetype = 'linux'
            bc.env_vars = test_env
            bc.name = name
            bc.conda_channels = ['http://ssb.stsci.edu/astroconda-dev']
            bc.conda_packages = conda_packages + ["python=${python_ver}"] + ["numpy=${numpy_ver}"]
            bc.build_cmds = [
                "pip install -q ${pip_packages_tests}",
                "python setup.py install"
            ]
            bc.test_cmds = ["pytest -r s --basetemp=test_results --junitxml=results.xml"]
            matrix += bc
        }
    }
}
*/
utils.run(matrix)
