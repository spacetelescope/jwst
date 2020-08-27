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
)
