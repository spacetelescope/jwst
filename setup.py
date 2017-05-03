import os
import subprocess
import sys
from setuptools import setup, find_packages, Extension
from setuptools.command.test import test as TestCommand
from numpy import get_include as np_include
from glob import glob

# hack building the sphinx docs with C source
from setuptools.command.build_ext import build_ext
import sphinx
from sphinx.setup_command import BuildDoc


class BuildSphinx(BuildDoc):
    """Build Sphinx documentation after compiling C source files"""

    description = 'Build Sphinx documentation'

    def initialize_options(self):
        BuildDoc.initialize_options(self)

    def finalize_options(self):
        BuildDoc.finalize_options(self)

    def run(self):
        build_cmd = self.reinitialize_command('build_ext')
        build_cmd.inplace = 1
        self.run_command('build_ext')
        sphinx.build_main(['setup.py', '-b', 'html', './docs', './docs/_build/html'])



NAME = 'jwst'
SCRIPTS = glob('scripts/*')
PACKAGE_DATA = {
    '': [
        '*.fits',
        '*.txt',
        '*.inc',
        '*.cfg',
        '*.csv',
        '*.yaml',
        '*.json'
    ]
}


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


if os.path.exists('relic'):
    sys.path.insert(1, 'relic')
    import relic.release
else:
    try:
        import relic.release
    except ImportError:
        try:
            subprocess.check_call(['git', 'clone',
                'https://github.com/jhunkeler/relic.git'])
            sys.path.insert(1, 'relic')
            import relic.release
        except subprocess.CalledProcessError as e:
            print(e)
            exit(1)


version = relic.release.get_info()
relic.release.write_template(version, NAME)


setup(
    name=NAME,
    version=version.pep386,
    author='OED/SSB, etc',
    author_email='help@stsci.edu',
    description='JWST',
    url='http://ssb.stsci.edu',
    license='BSD',
    classifiers=[
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    scripts=SCRIPTS,
    packages=find_packages(),
    package_data=PACKAGE_DATA,
    ext_modules=[
        Extension('jwst.tweakreg.chelp',
            glob('src/tweakreg/*.c'),
            include_dirs=[np_include()],
            define_macros=[('NUMPY', '1')]),
    ],
    tests_require=[
        'backports.tempfile',
        'pytest',
        'requests_mock',
        'pytest-catchlog'
    ],
    cmdclass={
        'test': PyTest,
        'build_sphinx': BuildSphinx
    },
)
