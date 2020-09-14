import os
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
PACKAGE_DATA['jwst.transforms'] = get_transforms_data()

setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    scripts=SCRIPTS,
    packages=find_packages(),
    package_data=PACKAGE_DATA,

    setup_requires=[
        'setuptools_scm',
    ],
    install_requires=[
        'asdf>=2.5',
        'astropy>=4.0',
        'crds>=7.4.1.3',
        'drizzle>=1.13',
        'gwcs>=0.13.0',
        'jsonschema>=3.0.1',
        'numpy>=1.16',
        'photutils>=0.7',
        'poppy>=0.9.1',
        'pyparsing>=2.2',
        'requests>=2.22',
        'scipy>=1.1.0',
        'spherical-geometry>=1.2.2',
        'stsci.image>=2.3.3',
        'tweakwcs>=0.6.4',
        'uncertainties>=3.1.4',
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
