import os
import sys
import pkgutil
from os.path import basename
from subprocess import check_call, CalledProcessError
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.test import test as TestCommand
from setuptools.command.build_ext import build_ext
from glob import glob

try:
    from numpy import get_include as np_include
except ImportError:
    print('Unable to import "numpy".\n'
          'Please install "numpy" and try again.', file=sys.stderr)
    exit(1)

if sys.version_info < (3, 5):
    error = """
    JWST 0.9+ does not support Python 2.x, 3.0, 3.1, 3.2, 3.3 or 3.4.
    Beginning with JWST 0.9, Python 3.5 and above is required.

    This may be due to an out of date pip

    Make sure you have pip >= 9.0.1.

    """
    sys.exit(error)


try:
    from sphinx.cmd.build import build_main
    from sphinx.setup_command import BuildDoc

    class BuildSphinx(BuildDoc):
        """Build Sphinx documentation after compiling C source files"""

        description = 'Build Sphinx documentation'

        user_options = BuildDoc.user_options[:]

        user_options.append(
            ('keep-going', 'k',
             'Parses the sphinx output and sets the return code to 1 if there '
             'are any warnings. Note that this will cause the sphinx log to '
             'only update when it completes, rather than continuously as is '
             'normally the case.'))


        def initialize_options(self):
            BuildDoc.initialize_options(self)

        def finalize_options(self):
            BuildDoc.finalize_options(self)

        def run(self):
            build_cmd = self.reinitialize_command('build_ext')
            build_cmd.inplace = 1
            self.run_command('build_ext')
            retcode = build_main(['-W', '--keep-going', '-b', 'html', './docs', './docs/_build/html'])
            if retcode != 0:
                sys.exit(retcode)

except ImportError:
    class BuildSphinx(Command):
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
            print('!\n! Sphinx is not installed!\n!', file=sys.stderr)
            exit(1)


NAME = 'jwst'
SCRIPTS = [s for s in glob('scripts/*') if basename(s) != '__pycache__']
PACKAGE_DATA = {
    '': [
        '*.fits',
        '*.txt',
        '*.inc',
        '*.cfg',
        '*.csv',
        '*.yaml',
        '*.json',
        '*.asdf'
    ]
}
DOCS_REQUIRE = [
    'matplotlib',
    'sphinx',
    'sphinx-automodapi',
    'sphinx-rtd-theme',
    'stsci-rtd-theme',
    'sphinx-astropy',
]
TESTS_REQUIRE = [
    'ci-watson>=0.3.0',
    'pytest',
    'pytest-doctestplus',
    'requests_mock',
    'pytest-astropy',
]

def get_transforms_data():
    # Installs the schema files in jwst/transforms
    # Because the path to the schemas includes "stsci.edu" they
    # can't be installed using setuptools.
    transforms_schemas = []
    root = os.path.join(NAME, 'transforms', 'schemas')
    for node, dirs, files in os.walk(root):
        for fname in files:
            if fname.endswith('.yaml'):
                transforms_schemas.append(
                    os.path.relpath(os.path.join(node, fname), root))
    # In the package directory, install to the subdirectory 'schemas'
    transforms_schemas = [os.path.join('schemas', s) for s in transforms_schemas]
    return transforms_schemas


transforms_schemas = get_transforms_data()
PACKAGE_DATA['jwst.transforms'] = transforms_schemas


class PyTest(TestCommand):

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = [NAME]

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        try:
            import pytest
        except ImportError:
            print('Unable to run tests...')
            print('To continue, please install "pytest":')
            print('    pip install pytest')
            exit(1)

        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


if not pkgutil.find_loader('relic'):
    relic_local = os.path.exists('relic')
    relic_submodule = (relic_local and
                       os.path.exists('.gitmodules') and
                       not os.listdir('relic'))
    try:
        if relic_submodule:
            check_call(['git', 'submodule', 'update', '--init', '--recursive'])
        elif not relic_local:
            check_call(['git', 'clone', 'https://github.com/spacetelescope/relic.git'])

        sys.path.insert(1, 'relic')
    except CalledProcessError as e:
        print(e)
        exit(1)

import relic.release

version = relic.release.get_info()
relic.release.write_template(version, NAME)

entry_points = dict(asdf_extensions=['jwst_pipeline = jwst.transforms.jwextension:JWSTExtension'])

setup(
    name=NAME,
    version=version.pep386,
    author='JWST Pipeline developers',
    description='Python library for science observations from the James Webb Space Telescope',
    long_description=('The JWST Data Reduction Pipeline is a Python '
                      'software suite that automatically processes the '
                      'data taken by the JWST instruments NIRCam, NIRSpec, '
                      'NIRISS, MIRI, and FGS to remove instrumental signatures '
                      'from the observations.'),
    url='https://github.com/spacetelescope/jwst',
    license='BSD',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: C',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    python_requires='>=3.5',
    scripts=SCRIPTS,
    packages=find_packages(),
    package_data=PACKAGE_DATA,
    ext_modules=[
        Extension('jwst.tweakreg.chelp',
                  glob('src/tweakreg/*.c'),
                  include_dirs=[np_include()],
                  define_macros=[('NUMPY', '1')]),
    ],
    install_requires=[
        'asdf>=2.3.2',
        'astropy>=3.1',
        'crds>=7.2.7',
        'drizzle>=1.12',
        'gwcs>=0.10',
        'jsonschema>=2.3,<=2.6',
        'numpy>=1.13',
        'photutils>=0.4',
        'scipy>=1.0',
        'spherical-geometry>=1.2',
        'stsci.image>=2.3',
        'stsci.imagestats>=1.4',
        'stsci.stimage>=0.2',
        'verhawk',
    ],
    extras_require={
        'docs': DOCS_REQUIRE,
        'ephem': ['pymssql>=2.1', 'jplephem>=2.8'],
        'test': TESTS_REQUIRE,
    },
    tests_require=TESTS_REQUIRE,
    cmdclass={
        'test': PyTest,
        'build_sphinx': BuildSphinx
    },
    entry_points=entry_points,
)
