#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='nf-wgs_amr',
    __version__="0.1",

    description='xxx',

    url='https://github.com/jhayer/nf-wgs_amr',
    download_url='https://github.com/jhayer/nf-wgs_amr/archive/v' + __version__ +'.tar.gz',
    author='Juliette Hayer, Jacques Dainat',

    license='GPL-3.0',
    packages=find_packages(),

    install_requires=['pyyaml', 'gitpython', 'requests', 'tarfile', 'gzip', 'biopython>=1.78', 'numpy>=1.22', 'python_version>="3.8.0"' ],
    include_package_data=True,

)
