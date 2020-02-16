"""Install script"""
from pyimzml import __version__
from setuptools import setup


DESCRIPTION = "Reader and writer of imzML 1.1.0 files"
VERSION = __version__
with open("README.md") as f:
    LONG_DESCRIPTION = f.read()

PACKAGE_NAME = "pyimzML"
MAINTAINER = "Lukasz G. Migas"
MAINTAINER_EMAIL = "lukas.migas@yahoo.com"
URL = "https://github.com/lukasz-migas/pyimzML"
LICENSE = "Apache license 2.0"
DOWNLOAD_URL = "https://github.com/lukasz-migas/pyimzML"
with open("requirements/requirements-std.txt") as f:
    INSTALL_REQUIRES = f.readlines()

CLASSIFIERS = [
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: Apache Software License",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
    "Operating System :: Microsoft :: Windows",
    "Natural Language :: English",
]

setup(
    name=PACKAGE_NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    download_url=DOWNLOAD_URL,
    url=URL,
    author="Dominik Fay",
    author_email="dominik.fay@embl.de",
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    keywords="bioinformatics imaging mass spectrometry parser imzML",
    packages=["pyimzml"],
    install_requires=INSTALL_REQUIRES,
)
