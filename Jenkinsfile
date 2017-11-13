if (utils.scm_checkout()) return
def ENV_SETUP = [
    "CRDS_SERVER_URL=https://jwst-crds.stsci.edu",
    "CRDS_PATH=./crds_cache",
]
def PY_VERSIONS = ['3.5', '3.6']
def NPY_VERSIONS = ['1.13', '1.14']
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
def CONDA_CHANNEL = "http://ssb.stsci.edu/astroconda-dev"
def PY_SETUP = "python setup.py"
def PYTEST = "pytest -r s -v --basetemp=./test_results --junit-xml=results.xml"
def matrix = []


for (py in PY_VERSIONS) {
    for (npy in NPY_VERSIONS) {
        def NAME = "py${py}np${npy}"
        def WRAPPER = "with_env -n ${NAME}"

        bc = new BuildConfig()
        bc.nodetype = 'linux-stable'
        bc.env_vars = ENV_SETUP
        bc.build_mode = NAME
        bc.build_cmds = [
            "conda config --add channels ${CONDA_CHANNEL}",
            "conda create -y -q -n ${NAME} python=${py} numpy=${npy} ${CONDA_DEPS}",
            "${WRAPPER} ${PY_SETUP} develop"
        ]
        bc.test_cmds = [
            "${WRAPPER} ${PYTEST}"
        ]
        matrix += bc
    }
}

utils.run(matrix)
