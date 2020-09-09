from os.path import basename
from setuptools import setup, find_packages
from glob import glob


scripts = [s for s in glob('scripts/*') if basename(s) != '__pycache__']
package_data = {
    '': [
        '*.fits', '**/*.fits','**/**/*.fits','**/**/**/*.fits',
        '*.txt',
        '*.inc',
        '*.cfg',
        '*.csv', '**/*.csv','**/**/*.csv','**/**/**/*.csv',
        '*.yaml',
        '*.json', '**/*.json','**/**/*.json','**/**/**/*.json',
        '*.asdf',
        '*.ecsv',
        '*.prop', '**/*.prop','**/**/*.prop','**/**/**/*.prop',
        '*.db', '**/*.db','**/**/*.db','**/**/**/*.db',
    ],
    "jwst.transforms": ['*.yaml', '**/*.yaml', '**/**/*.yaml', '**/**/**/*.yaml'],
}

setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    scripts=scripts,
    packages=find_packages(),
    package_data=package_data,
)
