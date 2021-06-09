#!/usr/bin/env python3
# encoding: utf-8

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

with open('README.md') as file:
    readme = file.read()

# Setup
setup(
    name='dnassembly',
    version='0.1.9',
    author='James Lucas',
    author_email='james.lucas@berkeley.edu',
    description='Digest and assemble DNA',
    long_description=readme,
    url='https://github.com/jaaamessszzz/DNAssembly',
    keywords=[
        'DNA',
        'assembly',
        'cloning',
        'plasmid',
        'digest',
        'restriction',
        'enzyme',
        'molecular',
        'biology'
    ],
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    packages=find_packages(),
    install_requires=[
        'networkx',
        'biopython'
    ],
    entry_points={
        'console_scripts': [
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
