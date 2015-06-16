from setuptools import setup, find_packages
import os

VERSION = '0.0.x'
LONG_DESCRIPTION = os.linesep.join([open('README.md').read()])

setup(
    name='py-vcf-parser',
    version=VERSION,

    author='Kensuke Numakura',
    author_email='knmkr3gma+pip@gmail.com',

    description='Simple VCF (variant call format) parser for python.',
    long_description=LONG_DESCRIPTION,
    url='',
    classifiers=[
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 2 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],

    packages=find_packages(),
)
