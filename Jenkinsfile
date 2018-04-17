if (utils.scm_checkout()) return
def ENV_SETUP = [
    "CRDS_SERVER_URL=https://jwst-crds.stsci.edu",
    "CRDS_PATH=./crds_cache",
]
def PY_VERSIONS = ['3.6']
def NPY_VERSIONS = ['1.14']
def ASTROPY_VERSIONS = ['3']
def CONDA_CHANNEL = "http://ssb.stsci.edu/astroconda-dev"
def CONDA_DEPS = "asdf \
                  astropy \
                  crds \
                  dask \
                  drizzle \
                  fitsblender \
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
def CONDA_ARGS = "-q -y"
def CONDA_CREATE = "conda create ${CONDA_ARGS}"
def CONDA_INST = "conda install ${CONDA_ARGS}"
def PIP_INST = "pip install -q"
def PIP_DEPS = ""
def PIP_DOC_DEPS = "sphinx-automodapi"
def PIP_TEST_DEPS = "requests_mock"
def PY_SETUP = "python setup.py"
def PYTEST = "pytest -r s -v --basetemp=./test_results --junit-xml=results.xml"
def matrix = []


sdist = new BuildConfig()
sdist.nodetype = 'linux-stable'
sdist.build_mode = 'dist'
sdist.build_cmds = [
    "${CONDA_INST} numpy",
    "${PY_SETUP} sdist"
]
matrix += sdist


docs = new BuildConfig()
docs.nodetype = 'linux-stable'
docs.build_mode = 'docs'
docs.build_cmds = [
    "conda config --add channels ${CONDA_CHANNEL}",
    "${CONDA_CREATE} -n ${docs.build_mode} ${CONDA_DOC_DEPS} ${CONDA_DEPS}",
    "with_env -n ${docs.build_mode} ${PIP_INST} ${PIP_DOC_DEPS}",
    "with_env -n ${docs.build_mode} ${PY_SETUP} build_sphinx"
]
matrix += docs


for (py in PY_VERSIONS) {
    for (npy in NPY_VERSIONS) {
        for (apy in ASTROPY_VERSIONS) {
            def NAME = "py${py}np${npy}ap${apy}"
            def WRAPPER = "with_env -n ${NAME}"

            bc = new BuildConfig()
            bc.nodetype = 'linux-stable'
            bc.env_vars = ENV_SETUP
            bc.build_mode = NAME
            bc.build_cmds = [
                "conda config --add channels ${CONDA_CHANNEL}",
                "${CONDA_CREATE} -n ${NAME} \
                    python=${py} numpy=${npy} astropy=${apy} ${CONDA_DEPS}",
                "${WRAPPER} ${PY_SETUP} develop"
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
