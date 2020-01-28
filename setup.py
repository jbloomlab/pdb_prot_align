"""Setup script for ``pdb_prot_align``."""


import re
import sys

from setuptools import setup


if not (sys.version_info[0] == 3 and sys.version_info[1] >= 6):
    raise RuntimeError(
                'pdb_prot_align requires Python >=3.6.\n'
                f"You are using {sys.version_info[0]}.{sys.version_info[1]}.")

# get metadata from package `__init__.py` file as here:
# https://packaging.python.org/guides/single-sourcing-package-version/
metadata = {}
init_file = 'pdb_prot_align/__init__.py'
with open(init_file) as f:
    init_text = f.read()
for dataname in ['version', 'author', 'email', 'url']:
    matches = re.findall(
            '__' + dataname + r'__\s+=\s+[\'"]([^\'"]+)[\'"]',
            init_text)
    if len(matches) != 1:
        raise ValueError(f"found {len(matches)} matches for {dataname} "
                         f"in {init_file}")
    else:
        metadata[dataname] = matches[0]

with open('README.rst') as f:
    readme = f.read()

# main setup command
setup(
    name='pdb_prot_align',
    version=metadata['version'],
    author=metadata['author'],
    author_email=metadata['email'],
    url=metadata['url'],
    download_url='https://github.com/jbloomlab/pdb_prot_align/tarball/' +
                 metadata['version'],  # tagged version on GitHub
    description='Align proteins to PDB chains',
    # Install "pdb_prot_align" program, calls pdb_prot_align.__main__.main()
    # https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
    entry_points={
            'console_scripts': [
                    'pdb_prot_align = pdb_prot_align.main:_run_main',
                    ],
            },
    long_description=readme,
    license='GPLv3',
    install_requires=[
            'biopython>=1.74',
            'natsort>=6.0.0',
            'pandas>=0.25.3',
            ],
    platforms='Linux and Mac OS X.',
    packages=['pdb_prot_align'],
    package_dir={'pdb_prot_align': 'pdb_prot_align'},
    )
