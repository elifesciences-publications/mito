#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0',
    'numpy>=1.9',
    'scipy>=0.15',
    'Cython>=0.25.2',
    'scikit-image>=0.12.3',
    'pandas>=0.19.2',
    'seaborn>=0.7.1',
    'javabridge==1.0.15',
    'python-bioformats>=1.0.8',
    'pyome>=0.1.0'
]

test_requirements = [
    'nose>=1.3',
    'numpy>=1.9',
    'scipy>=0.15',
    'Cython>=0.25.2',
    'scikit-image>=0.12.3',
    'pandas>=0.19.2',
    'seaborn>=0.7.1',
    'javabridge==1.0.15',
    'python-bioformats>=1.0.8',
    'pyome>=0.1.0'
]

setup(
    name='mito',
    version='0.1.0',
    description="Tools for analysing mitochondrial distribution in fluoresence microscopy images",
    long_description=readme + '\n\n' + history,
    author="Keith Schulze",
    author_email='keith.schulze@monash.edu',
    url='https://gitlab.erc.monash.edu.au/skeith/mito.git',
    packages=[
        'mito',
    ],
    package_dir={'mito':
                 'mito'},
    entry_points={
        'console_scripts': [
            'mito=mito.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    dependency_links=[
        'https://github.com/LeeKamentsky/python-javabridge/archive/master.zip#egg=javabridge-1.0.15',
        'https://github.com/keithschulze/pyome/archive/master.zip#egg=pyome-0.1.0'
    ],
    license="MIT license",
    zip_safe=False,
    keywords='mito',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
