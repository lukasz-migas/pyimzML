# pyimzML

[![Build Status](https://readthedocs.org/projects/pyimzml/badge/?version=latest&style=flat)](https://readthedocs.org/projects/pyimzml/badge/?version=latest)
![Tests](https://github.com/lukasz-migas/pyimzML/workflows/Tests/badge.svg)
![Style](https://github.com/lukasz-migas/pyimzML/workflows/Style/badge.svg)

**FORK**: [https://github.com/alexandrovteam/pyimzML](https://github.com/alexandrovteam/pyimzML)

**ORIGINAL AUTHOR**: Dominik Fay / [Alexandrov Team](https://github.com/alexandrovteam) 

## Description

A parser for the imzML format used in imaging mass spectrometry. See [specification](http://imzml.org/download/imzml/specifications_imzML1.1.0_RC1.pdf)

Designed for use with imzML version 1.1.0. Outputs data as python lists and dicts.

## What's different from the original version

### Implemented

- changed the test framework
- added GitHub workflows for CI testing
- changed few method names
- tidied-up code and renamed few unclear variables

### Planned

- multi-thread/core support for extracting ion images
- speed-up XML parsing (if possible?)
- add support for parsing compressed imzML files as currently decompression is omitted by the parser
- changed the structure of byte offset arrays to single 2D array

## Installation

pyimzML is available on [PyPI](https://pypi.python.org/pypi/pyimzML) where you can install it using several different ways

Simply run the `pip install` command which will install the version developed by [Alexandrov Team](https://github.com/alexandrovteam/pyimzML)

```python
pip install pyimzml
```

The updated version is available from this GitHub repository and can be installed by invoking

```python
pip install git+git://github.com/lukasz-migas/pyimzML.git 
```

## Dependencies

pyimzML has an optional dependency to [lxml](http://lxml.de/index.html). If `lxml` is not installed, pyimzML will instead use the built-in cElementTree or ElementTree package.

## Testing

Tests have been rewritten and use `pytest` framework. You can find them inside the `tests/` directory. If you would 
like to run the tests, make sure to install the additional dependencies (see `requirements/requirements-dev.txt`).

```python
pip install -r requirements/requirements-dev.txt
pytest .
```

## Documentation

Documentation is available on [ReadTheDocs](http://pyimzml.readthedocs.org/en/latest)