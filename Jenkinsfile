// Perform initial clone and enable [skip ci] feature
if (utils.scm_checkout()) return

// Build matrix setup
def PY_VERSIONS = ['3.6']
def NPY_VERSIONS = ['1.14']
def ASTROPY_VERSIONS = ['4']  // dev channel is (major + 1)
def matrix = []  // used by 'utils'

// Shell environment setup
// NOTE: './' or '.' are replaced at runtime with $WORKSPACE
def ENV_SETUP = [
    "CRDS_SERVER_URL=https://jwst-crds.stsci.edu",
    "CRDS_PATH=./crds_cache",
]

// Conda related setup
def CONDA_ARGS = "-q -y"
def CONDA_CREATE = "conda create ${CONDA_ARGS}"
def CONDA_INST = "conda install ${CONDA_ARGS}"
def CONDA_CHANNEL = "http://ssb.stsci.edu/astroconda-dev"
def CONDA_DEPS = "asdf \
                  astropy \
                  crds \
                  dask \
                  drizzle \
                  flake8 \
                  gwcs \
                  jsonschema \
                  jplephem \
                  matplotlib \
                  namedlist \
                  numpy \
                  photutils \
                  scipy \
                  six \
                  spherical-geometry \
                  stsci.image \
                  stsci.imagestats \
                  stsci.stimage \
                  stsci.tools \
                  verhawk"
def CONDA_DOC_DEPS = "sphinx \
                      sphinx_rtd_theme \
                      stsci_rtd_theme"
def CONDA_TEST_DEPS = "pytest"

// Pip related setup
def PIP_ARGS = "-q"
def PIP_INST = "pip install ${PIP_ARGS}"
def PIP_DEPS = ""
def PIP_DOC_DEPS = "sphinx-automodapi"
def PIP_TEST_DEPS = "requests_mock git+https://github.com/spacetelescope/ci_watson.git"

// Pytest wrapper
def PYTEST = "pytest \
              -r s \
              -v \
              --basetemp=./test_results \
              --junit-xml=results.xml"

// Python related setup
def PY_SETUP = "python setup.py"


//
// JOB SPECIFICATION
//

// Generate distributions
dist = new BuildConfig()
dist.nodetype = 'linux'
dist.name = 'dist'
dist.build_cmds = [
    "${CONDA_INST} numpy",
    "${PY_SETUP} sdist",
    "${PY_SETUP} bdist_egg",
    "${PY_SETUP} bdist_wheel"
]
matrix += dist


// Compile documentation
docs = new BuildConfig()
docs.nodetype = 'linux'
docs.name = 'docs'
docs.build_cmds = [
    "conda config --add channels ${CONDA_CHANNEL}",
    "${CONDA_CREATE} -n ${docs.name} ${CONDA_DOC_DEPS} ${CONDA_DEPS}",
    "with_env -n ${docs.name} ${PIP_INST} ${PIP_DOC_DEPS}",
    "with_env -n ${docs.name} ${PY_SETUP} build_sphinx"
]
matrix += docs


// Generate the build and test matrix
for (py in PY_VERSIONS) {
    for (npy in NPY_VERSIONS) {
        for (apy in ASTROPY_VERSIONS) {
            def NAME = "py${py}np${npy}ap${apy}"
            def WRAPPER = "with_env -n ${NAME}"

            bc = new BuildConfig()
            bc.nodetype = 'linux'
            bc.env_vars = ENV_SETUP
            bc.name = NAME
            bc.build_cmds = [
                "conda config --add channels ${CONDA_CHANNEL}",
                "${CONDA_CREATE} -n ${NAME} \
                    python=${py} numpy=${npy} astropy=${apy} ${CONDA_DEPS}",
                "${WRAPPER} ${PY_SETUP} install"
            ]
            bc.test_cmds = [
                "${WRAPPER} ${CONDA_INST} ${CONDA_TEST_DEPS}",
                "${WRAPPER} ${PIP_INST} ${PIP_TEST_DEPS}",
                "${WRAPPER} ${PYTEST}"
            ]
            matrix += bc
        }
    }
}

utils.run(matrix)
