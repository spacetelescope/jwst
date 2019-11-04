import os
import sys
from os.path import basename
from setuptools import setup, find_packages
from glob import glob


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
    'sphinx-asdf',
]
TESTS_REQUIRE = [
    'ci-watson>=0.3.0',
    'pytest',
    'pytest-doctestplus',
    'requests_mock',
    'pytest-openfiles',
    'pytest-cov',
    'codecov',
]
AWS_REQUIRE = [
    'stsci-aws-utils @ git+https://github.com/spacetelescope/stsci-aws-utils@0.1.1'
]
ENTRY_POINTS = dict(asdf_extensions=['jwst_pipeline = jwst.transforms.jwextension:JWSTExtension',
                                     'jwst_datamodel = jwst.datamodels.extension:DataModelExtension'])

transforms_schemas = get_transforms_data()
PACKAGE_DATA['jwst.transforms'] = transforms_schemas

setup(
    name=NAME,
    use_scm_version=True,
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
    python_requires='>=3.6',
    scripts=SCRIPTS,
    packages=find_packages(),
    package_data=PACKAGE_DATA,
    setup_requires=[
        'setuptools_scm',
    ],
    install_requires=[
        'asdf>=2.4',
        'astropy @ git+https://github.com/nden/astropy@11386e2df7af8',
        'crds>=7.2.7',
        'drizzle>=1.13',
        'gwcs @ git+https://github.com/spacetelescope/gwcs@ace1c2c30a658',
        'jsonschema>=2.3,<4',
        'numpy>=1.16',
        'photutils>=0.7',
        'scipy>=1.0',
        'spherical-geometry>=1.2',
        'stsci.image>=2.3.3',
        'stsci.imagestats>=1.4',
        'tweakwcs>=0.5.1',
    ],
    extras_require={
        'docs': DOCS_REQUIRE,
        'ephem': ['pymssql==2.1.4', 'jplephem==2.9'], # for timeconversion
        'test': TESTS_REQUIRE,
        'aws': AWS_REQUIRE,
    },
    tests_require=TESTS_REQUIRE,
    entry_points=ENTRY_POINTS,
)
