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

### Main changes

- renamed few functions to follow the `snake_case` style. Where appropriate, the old method names remained and alias functions have been created
- the `PortableSpectrumReader` parser now has nearly identical interface to the `ImzMLParser`. This was done by replacing the
way binary (`.ibd`) data was read from handle to `with open(FILENAME.ibd) as ibd_handle` meaning that it can be easily pickled. This also enables 
multicore support to ion image extraction
- moved the `browse` and the two classes to separate file as it didn't fit with the parsers
- the `ImzMLParser` now has two additonal keyword parameters (`as_threads` and `pool_size` which control how multicore/thread)
data is read
- the `ImzMLParser` no longer uses handles to the `ibd_file` but simply extracts the filename (attribute `.name`) and handles
it like that 

### Implemented so far

- changed the test framework
- added GitHub workflows for CI testing
- changed few method names
- tidied-up code and renamed few unclear variables
- multi-thread/core support for extracting ion images (see `get_async_ion_image`)

### Planned features

- speed-up XML parsing (if possible?)
- add support for parsing compressed imzML files as currently decompression is omitted by the parser
- changed the structure of byte offset arrays to single 2D array - somewhat tricky since we do not know the number of pixels
in the the dataset so pre-allocation of large array is not convenient

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