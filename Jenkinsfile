if (utils.scm_checkout()) return

matrix_python = ['3.6']
matrix_numpy = ['1.15']
matrix_astropy = ['4']
matrix_astropy_pip = ['3.1']
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
bc_dist = new BuildConfig()
bc_dist.nodetype = 'linux'
bc_dist.name = 'dist'
bc_dist.conda_packages = ["python=${matrix_python[0]}"]
bc_dist.build_cmds = [
    "pip install ${pip_install_args} numpy==${matrix_numpy[0]}",
    "pip wheel ${pip_install_args} .",
    "python setup.py sdist",
]
matrix += bc_dist

// Generate pip build and test matrix with released upstream dependencies
for (python_ver in matrix_python) {
    for (numpy_ver in matrix_numpy) {
        for (astropy_ver in matrix_astropy_pip) {
            def name = "pip_py${python_ver}np${numpy_ver}ap${astropy_ver}"
            bc = new BuildConfig()
            bc.nodetype = 'linux'
            bc.env_vars = test_env
            bc.name = name
            bc.conda_packages = ["python=${python_ver}", "git"]
            bc.build_cmds = [
                "pip install ${pip_install_args} numpy",
                "pip install ${pip_install_args} -e .[test]",
                "git clean -xdf",
                "python setup.py develop",
            ]
            bc.test_cmds = ["pytest -r s --basetemp=test_results --junitxml=results.xml"]
            matrix += bc
        }
    }
}

// Generate conda build and test matrix with astroconda-dev dependencies
for (python_ver in matrix_python) {
    for (numpy_ver in matrix_numpy) {
        for (astropy_ver in matrix_astropy) {
            def name = "conda_py${python_ver}np${numpy_ver}ap${astropy_ver}"
            bc = new BuildConfig()
            bc.nodetype = 'linux'
            bc.env_vars = test_env
            bc.name = name
            bc.conda_channels = ['http://ssb.stsci.edu/astroconda-dev']
            bc.conda_packages = conda_packages + ["python=${python_ver}", "numpy=${numpy_ver}", git]
            bc.build_cmds = [
                "pip install ${pip_install_args} -e .[test]",
                "git clean -xdf",
                "python setup.py develop",
            ]
            bc.test_cmds = ["pytest -r s --basetemp=test_results --junitxml=results.xml"]
            matrix += bc
        }
    }
}

utils.run(matrix)
