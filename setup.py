import os
import subprocess
import sys
from setuptools import setup, find_packages, Extension
from numpy import get_include as np_include
from glob import glob


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
relic.release.write_template(version, 'jwst/')


SCRIPTS=glob('scripts/*')
PACKAGE_DATA={
    '': ['*.fits', 
        '*.txt',
        '*.inc',
        '*.json',
        '*.cfg']
}

setup(
    name='jwst',
    version=version.pep386,
    author='OED/SSB, etc',
    author_email='help@stsci.edu',
    description='JWST',
    url='http://ssb.stsci.edu',
    license='BSD',
    classifiers = [
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
            define_macros=[('NUMPY','1')]),
    ],
)
    
