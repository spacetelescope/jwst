import os
from os.path import basename
import subprocess
import sys
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.test import test as TestCommand
from numpy import get_include as np_include
from glob import glob

if sys.version_info < (3, 5):
    error = """
    JWST 0.9+ does not support Python 2.x, 3.0, 3.1, 3.2, 3.3 or 3.4.
    Beginning with JWST 0.9, Python 3.5 and above is required.

    This may be due to an out of date pip

    Make sure you have pip >= 9.0.1.

    """
    sys.exit(error)

# hack building the sphinx docs with C source
from setuptools.command.build_ext import build_ext

try:
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

entry_points = dict(asdf_extensions='jwst_pipeline = jwst.transforms.jwextension:JWSTExtension')

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
        'namedlist'
    ],
    tests_require=[
        'pytest',
        'requests_mock'
    ],
    cmdclass={
        'test': PyTest,
        'build_sphinx': BuildSphinx
    },
    entry_points=entry_points,
)
