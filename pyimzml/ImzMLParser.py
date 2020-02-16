# -*- coding: utf-8 -*-

# Copyright 2015 Dominik Fay
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Standard library imports
import math
from bisect import bisect_left, bisect_right
import re
from pathlib import Path
from warnings import warn
import logging
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import wait

# Third-party imports
import numpy as np

PRECISION_DICT = {
    "32-bit float": "f",
    "64-bit float": "d",
    "32-bit integer": "i",
    "64-bit integer": "l",
}
SIZE_DICT = {"f": 4, "d": 8, "i": 4, "l": 8}

param_group_elname = "referenceableParamGroup"
data_processing_elname = "dataProcessing"
instrument_confid_elname = "instrumentConfiguration"

LOGGER = logging.getLogger(__name__)


def choose_iterparse(parse_lib=None):
    if parse_lib == "ElementTree":
        from xml.etree.ElementTree import iterparse
    elif parse_lib == "lxml":
        from lxml.etree import iterparse
    else:
        try:
            from lxml.etree import iterparse
        except ImportError:
            from xml.etree.ElementTree import iterparse
    return iterparse


class ImzMLParser:
    """
    Parser for imzML 1.1.0 files (see specification here:
    http://imzml.org/download/imzml/specifications_imzML1.1.0_RC1.pdf).

    Iteratively reads the .imzML file into memory while pruning the per-spectrum metadata (everything in
    <spectrumList> elements) during initialization. Returns a spectrum upon calling getspectrum(i). The binary file
    is read in every call of getspectrum(i). Use enumerate(parser.coordinates) to get all coordinates with their
    respective index. Coordinates are always 3-dimensional. If the third spatial dimension is not present in
    the data, it will be set to zero.

    *pyimzML* has limited support for the metadata embedded in the imzML file. For some general metadata, you can use
    the parser's ``ìmzmldict`` attribute. You can find the exact list of supported metadata in the documentation of the
    ``__readimzmlmeta`` method.
    """

    def __init__(
        self, filename, parse_lib=None, ibd_file=None, as_threads=False, pool_size=4
    ):
        """
        Opens the two files corresponding to the file name, reads the entire .imzML
        file and extracts required attributes. Does not read any binary data, yet.

        :param filename:
            name of the XML file. Must end with .imzML. Binary data file must be named equally but ending with .ibd
            Alternatively an open file or Buffer Protocol object can be supplied, if ibd_file is also supplied
        :param parse_lib:
            XML-parsing library to use: 'ElementTree' or 'lxml', the later will be used if argument not provided
        :param ibd_file:
            File or Buffer Protocol object for the .ibd file. Leave blank to infer it from the imzml filename.
            Set to None if no data from the .ibd file is needed (getspectrum calls will not work)
        """
        # custom map sizes are currently not supported, therefore mapsize is hardcoded.
        mapsize = 0
        # ElementTree requires the schema location for finding tags (why?) but
        # fails to read it from the root element. As this should be identical
        # for all imzML files, it is hard-coded here and prepended before every tag
        self.sl = "{http://psi.hupo.org/ms/mzml}"
        # maps each imzML number format to its struct equivalent
        self.precisionDict = dict(PRECISION_DICT)

        # maps each number format character to its amount of bytes used
        self.sizeDict = dict(SIZE_DICT)
        self.filename = filename

        self.mzOffsets = []
        self.intensityOffsets = []
        self.mzLengths = []
        self.intensityLengths = []
        self.byte_array = None

        # list of all (x,y,z) coordinates as tuples.
        self.coordinates = []
        self.root = None
        self.mzGroupId = (
            self.intGroupId
        ) = self.mzPrecision = self.intensityPrecision = None
        self._parse_lib = parse_lib
        # self.iterparse = choose_iterparse(parse_lib)

        # setup few attributes
        self._n_pixels = None
        self._idx = 0

        # select iterator
        self.__iter_read_spectrum_meta()

        # get binary data handle
        if ibd_file is None:
            # name of the binary file
            ibd_filename = self._infer_bin_filename(self.filename)
            self.ibd_filename = ibd_filename
        elif isinstance(ibd_file, str):
            self.ibd_filename = ibd_file
        else:
            if hasattr(ibd_file, "name"):
                self.ibd_filename = ibd_file.name
            else:
                raise ValueError(
                    "The `ibd_file` signature was changed. Please provide filename or open object"
                )

        # Dict for basic imzML metadata other than those required for reading spectra. See method __readimzmlmeta()
        self.imzmldict = self.__read_imzml_metadata()
        self.imzmldict["max count of pixels z"] = np.asarray(self.coordinates)[
            :, 2
        ].max()

        # multi-core support for image generation
        self._as_threads = as_threads
        self._pool_size = 2 if self._as_threads else pool_size

        # setup executor
        self.executor = (
            ThreadPoolExecutor(self._pool_size)
            if self._as_threads
            else ProcessPoolExecutor(self._pool_size)
        )
        self.root = None

    def __repr__(self):
        return f"ImzMLParser<filename: {self.filename}; no. pixels: {self.n_pixels}>"

    def __enter__(self):
        """Enter - system method used `with ... as`"""
        return self

    def __exit__(self, exc_t, exc_v, trace):
        """Exit - system method used `with ... as `"""

    def __iter__(self):
        return self

    def __next__(self):
        """Get next spectrum"""
        idx = self._idx
        if idx < self.n_pixels:
            self._idx += 1
            return self[idx]
        else:
            self._idx = 0
            raise StopIteration

    def __getitem__(self, item):
        """Retrieve spectrum"""
        try:
            return self.get_spectrum(item)
        except IndexError:
            LOGGER.warning(f"Could not retrieve {item}")

    @property
    def n_pixels(self):
        """Return number of pixels in the dataset"""
        if self._n_pixels is None:
            self._n_pixels = len(self.coordinates)
        return self._n_pixels

    @staticmethod
    def _infer_bin_filename(imzml_path):
        imzml_path = Path(imzml_path)
        ibd_path = [
            f
            for f in imzml_path.parent.glob("*")
            if re.match(r".+\.ibd", str(f), re.IGNORECASE) and f.stem == imzml_path.stem
        ][0]
        return str(ibd_path)

    def __iter_read_spectrum_meta(self):
        """
        This method should only be called by __init__. Reads the data formats, coordinates and offsets from
        the .imzML file and initializes the respective attributes. While traversing the XML tree, the per-spectrum
        metadata is pruned, i.e. the <spectrumList> element(s) are left behind empty.

        Supported accession values for the number formats: "MS:1000521", "MS:1000523", "IMS:1000141" or
        "IMS:1000142". The string values are "32-bit float", "64-bit float", "32-bit integer", "64-bit integer".
        """
        mz_group = int_group = None

        # get iterator
        iterparse = choose_iterparse(self._parse_lib)
        elem_iterator = iterparse(self.filename, events=("start", "end"))

        slist = None
        _, self.root = next(elem_iterator)
        for event, elem in elem_iterator:
            if elem.tag == self.sl + "spectrumList" and event == "start":
                slist = elem
            elif elem.tag == self.sl + "spectrum" and event == "end":
                self.__process_spectrum(elem)
                slist.remove(elem)
            elif elem.tag == self.sl + "referenceableParamGroup" and event == "end":
                for param in elem:
                    if param.attrib["name"] == "m/z array":
                        self.mzGroupId = elem.attrib["id"]
                        mz_group = elem
                    elif param.attrib["name"] == "intensity array":
                        self.intGroupId = elem.attrib["id"]
                        int_group = elem

        # cleanup
        self.__assign_precision(int_group, mz_group)
        self.__fix_offsets()

    def __fix_offsets(self):
        """Fix errors introduced by incorrect signed 32bit integers when unsigned 64bit was appropriate"""

        def fix(array):
            fixed = []
            delta = 0
            prev_value = float("nan")
            for value in array:
                if value < 0 <= prev_value:
                    delta += 2 ** 32
                fixed.append(value + delta)
                prev_value = value
            return fixed

        self.mzOffsets = fix(self.mzOffsets)
        self.intensityOffsets = fix(self.intensityOffsets)

    def __assign_precision(self, int_group, mz_group):
        valid_accession_strings = (
            "MS:1000521",
            "MS:1000523",
            "IMS:1000141",
            "IMS:1000142",
            "MS:1000519",
            "MS:1000522",
        )
        mz_precision = int_precision = None
        for s in valid_accession_strings:
            param = mz_group.find('%scvParam[@accession="%s"]' % (self.sl, s))
            if param is not None:
                mz_precision = self.precisionDict[param.attrib["name"]]
                break
        for s in valid_accession_strings:
            param = int_group.find('%scvParam[@accession="%s"]' % (self.sl, s))
            if param is not None:
                int_precision = self.precisionDict[param.attrib["name"]]
                break
        if (mz_precision is None) or (int_precision is None):
            raise RuntimeError(
                "Unsupported number format: mz = %s, int = %s"
                % (mz_precision, int_precision)
            )
        self.mzPrecision, self.intensityPrecision = mz_precision, int_precision

    def __process_spectrum(self, elem):
        array_list_item = elem.find("%sbinaryDataArrayList" % self.sl)
        element_list = list(array_list_item)
        element_list_sorted = [None, None]
        for element in element_list:
            ref = element.find("%sreferenceableParamGroupRef" % self.sl).attrib["ref"]
            if ref == self.mzGroupId:
                element_list_sorted[0] = element
            elif ref == self.intGroupId:
                element_list_sorted[1] = element

        mz_offset_elem = element_list_sorted[0].find(
            '%scvParam[@accession="IMS:1000102"]' % self.sl
        )
        self.mzOffsets.append(int(mz_offset_elem.attrib["value"]))

        mz_length_elem = element_list_sorted[0].find(
            '%scvParam[@accession="IMS:1000103"]' % self.sl
        )
        self.mzLengths.append(int(mz_length_elem.attrib["value"]))

        intensity_offset_elem = element_list_sorted[1].find(
            '%scvParam[@accession="IMS:1000102"]' % self.sl
        )
        self.intensityOffsets.append(int(intensity_offset_elem.attrib["value"]))

        intensity_length_elem = element_list_sorted[1].find(
            '%scvParam[@accession="IMS:1000103"]' % self.sl
        )
        self.intensityLengths.append(int(intensity_length_elem.attrib["value"]))

        scan_elem = elem.find("%sscanList/%sscan" % (self.sl, self.sl))
        x = scan_elem.find('%scvParam[@accession="IMS:1000050"]' % self.sl).attrib[
            "value"
        ]
        y = scan_elem.find('%scvParam[@accession="IMS:1000051"]' % self.sl).attrib[
            "value"
        ]
        try:
            z = scan_elem.find('%scvParam[@accession="IMS:1000052"]' % self.sl).attrib[
                "value"
            ]
            self.coordinates.append((int(x), int(y), int(z)))
        except AttributeError:
            self.coordinates.append((int(x), int(y), 1))

    def __read_imzml_metadata(self):
        """
        This method should only be called by __init__. Initializes the imzmldict with frequently used metadata from
        the .imzML file.

        This method reads only a subset of the available meta information and may be extended in the future. The keys
        are named similarly to the imzML names. Currently supported keys: "max dimension x", "max dimension y",
        "pixel size x", "pixel size y", "matrix solution concentration", "wavelength", "focus diameter x",
        "focus diameter y", "pulse energy", "pulse duration", "attenuation".

        If a key is not found in the XML tree, it will not be in the dict either.

        :return d:
            dict containing above mentioned meta data
        :rtype:
            dict
        :raises Warning:
            if an xml attribute has a number format different from the imzML specification
        """

        def check_meta(param, accession, elem_list):
            for idx, _ in enumerate(param):
                acc, attr = accession[idx]
                elem = elem_list.find('.//%scvParam[@accession="%s"]' % (self.sl, acc))
                if elem is None:
                    break
                name, T = param[idx]
                try:
                    metadata_dict[name] = T(elem.attrib[attr])
                except ValueError:
                    warn(
                        Warning(
                            'Wrong data type in XML file. Skipped attribute "%s"' % name
                        )
                    )

        metadata_dict = {}
        scan_settings_list_elem = self.root.find("%sscanSettingsList" % self.sl)
        instrument_config_list_elem = self.root.find(
            "%sinstrumentConfigurationList" % self.sl
        )
        supported_params_1 = [
            ("max count of pixels x", int),
            ("max count of pixels y", int),
            ("max dimension x", int),
            ("max dimension y", int),
            ("pixel size x", float),
            ("pixel size y", float),
            ("matrix solution concentration", float),
        ]
        supported_params_2 = [
            ("wavelength", float),
            ("focus diameter x", float),
            ("focus diameter y", float),
            ("pulse energy", float),
            ("pulse duration", float),
            ("attenuation", float),
        ]
        supported_accession_1 = [
            ("IMS:1000042", "value"),
            ("IMS:1000043", "value"),
            ("IMS:1000044", "value"),
            ("IMS:1000045", "value"),
            ("IMS:1000046", "value"),
            ("IMS:1000047", "value"),
            ("MS:1000835", "value"),
        ]
        supported_accession_2 = [
            ("MS:1000843", "value"),
            ("MS:1000844", "value"),
            ("MS:1000845", "value"),
            ("MS:1000846", "value"),
            ("MS:1000847", "value"),
            ("MS:1000848", "value"),
        ]
        check_meta(supported_params_1, supported_accession_1, scan_settings_list_elem)
        check_meta(
            supported_params_2, supported_accession_2, instrument_config_list_elem
        )
        return metadata_dict

    def get_physical_coordinates(self, i):
        """
        For a pixel index i, return the real-world coordinates in nanometers.

        This is equivalent to multiplying the image coordinates of the given pixel with the pixel size.

        :param i: the pixel index
        :return: a tuple of x and y coordinates.
        :rtype: Tuple[float]
        :raises KeyError: if the .imzML file does not specify the attributes "pixel size x" and "pixel size y"
        """
        try:
            pixel_size_x = self.imzmldict["pixel size x"]
            pixel_size_y = self.imzmldict["pixel size y"]
        except KeyError:
            raise KeyError("Could not find all pixel size attributes in imzML file")
        image_x, image_y = self.coordinates[i][:2]
        return image_x * pixel_size_x, image_y * pixel_size_y

    def get_spectrum(self, index):
        """
        Reads the spectrum at specified index from the .ibd file.

        :param index:
            Index of the desired spectrum in the .imzML file

        Output:

        mz_array: numpy.ndarray
            Sequence of m/z values representing the horizontal axis of the desired mass
            spectrum
        intensity_array: numpy.ndarray
            Sequence of intensity values corresponding to mz_array
        """
        mz_bytes, intensity_bytes = self._read_spectrum(index)
        mz_array = np.frombuffer(mz_bytes, dtype=self.mzPrecision)
        intensity_array = np.frombuffer(intensity_bytes, dtype=self.intensityPrecision)
        return mz_array, intensity_array

    # add alias
    getspectrum = get_spectrum

    def _read_spectrum(self, index):
        """
        Reads m/z array and intensity array of the spectrum at specified location
        from the binary file as a byte string. The string can be unpacked by the struct
        module. To get the arrays as numbers, use getspectrum

        :param index:
            Index of the desired spectrum in the .imzML file
        :rtype: Tuple[str, str]

        Output:

        mz_string:
            string where each character represents a byte of the mz array of the
            spectrum
        intensity_string:
            string where each character represents a byte of the intensity array of
            the spectrum
        """
        offsets = [self.mzOffsets[index], self.intensityOffsets[index]]
        lengths = [self.mzLengths[index], self.intensityLengths[index]]
        lengths[0] *= self.sizeDict[self.mzPrecision]
        lengths[1] *= self.sizeDict[self.intensityPrecision]
        with open(self.ibd_filename, "rb") as ibd_handle:
            ibd_handle.seek(offsets[0])
            mz_string = ibd_handle.read(lengths[0])
            ibd_handle.seek(offsets[1])
            intensity_string = ibd_handle.read(lengths[1])

        return mz_string, intensity_string

    def get_ion_image(self, mz_value, tol=0.1, z=1, reduce_func=sum):
        """Get an image representation of the intensity distribution of the ion with specified m/z value.

        By default, the intensity values within the tolerance region are summed.

        :param p:
            the ImzMLParser (or anything else with similar attributes) for the desired dataset
        :param mz_value:
            m/z value for which the ion image shall be returned
        :param tol:
            Absolute tolerance for the m/z value, such that all ions with values
            mz_value-|tol| <= x <= mz_value+|tol| are included. Defaults to 0.1
        :param z:
            z Value if spectrogram is 3-dimensional.
        :param reduce_func:
            the bahaviour for reducing the intensities between mz_value-|tol| and mz_value+|tol| to a single value. Must
            be a function that takes a sequence as input and outputs a number. By default, the values are summed.

        :return:
            numpy matrix with each element representing the ion intensity in this
            pixel. Can be easily plotted with matplotlib
        """
        return get_ion_image(self, mz_value, tol, z, reduce_func)

    def get_async_ion_image(self, mz_value, tol=0.1, z=1, reduce_func=sum):
        """Same as `get_ion_image`, however the action is executed in parallel using Thread/ProcessPool"""
        from pyimzml.utilities import chunks

        # create pickleable object
        portable = self.portable_spectrum_reader()

        # determine chunksize based on the pool executor
        chunk_size = math.ceil(self.n_pixels / self._pool_size)
        futures, start_idx = [], 0
        for coordinate_chunk in chunks(self.coordinates, chunk_size):
            futures.append(
                self.executor.submit(
                    _get_ion_image,
                    portable,
                    mz_value,
                    tol,
                    z,
                    reduce_func,
                    coordinate_chunk,
                    start_idx,
                )
            )
            start_idx += len(coordinate_chunk)

        # extract features
        futures_done, __ = wait(futures)

        # recreate image based on the results
        im = np.zeros(
            (
                self.imzmldict["max count of pixels y"],
                self.imzmldict["max count of pixels x"],
            )
        )
        for future in futures_done:
            im += future.result()
        return im

    def portable_spectrum_reader(self):
        """
        Builds a PortableSpectrumReader that holds the coordinates list and spectrum offsets in the .ibd file
        so that the .ibd file can be read without opening the .imzML file again.

        The PortableSpectrumReader can be safely pickled and unpickled, making it useful for reading the spectra
        in a distributed environment such as PySpark or PyWren.
        """
        return PortableSpectrumReader(
            self.coordinates,
            self.mzPrecision,
            self.mzOffsets,
            self.mzLengths,
            self.intensityPrecision,
            self.intensityOffsets,
            self.intensityLengths,
            self.ibd_filename,
            self.imzmldict,
        )


def get_ion_image(p, mz_value, tol=0.1, z=1, reduce_func=sum):
    """Get an image representation of the intensity distribution of the ion with specified m/z value.

    By default, the intensity values within the tolerance region are summed.

    :param p:
        the ImzMLParser (or anything else with similar attributes) for the desired dataset
    :param mz_value:
        m/z value for which the ion image shall be returned
    :param tol:
        Absolute tolerance for the m/z value, such that all ions with values
        mz_value-|tol| <= x <= mz_value+|tol| are included. Defaults to 0.1
    :param z:
        z Value if spectrogram is 3-dimensional.
    :param reduce_func:
        the bahaviour for reducing the intensities between mz_value-|tol| and mz_value+|tol| to a single value. Must
        be a function that takes a sequence as input and outputs a number. By default, the values are summed.

    :return:
        numpy matrix with each element representing the ion intensity in this
        pixel. Can be easily plotted with matplotlib
    """
    tol = abs(tol)
    im = np.zeros(
        (p.imzmldict["max count of pixels y"], p.imzmldict["max count of pixels x"])
    )
    for i, (x, y, z_) in enumerate(p.coordinates):
        if z_ == 0:
            UserWarning(
                "z coordinate = 0 present, if you're getting blank images set getionimage(.., .., z=0)"
            )

        if z_ == z:
            mzs, ints = map(lambda x: np.asarray(x), p.get_spectrum(i))
            min_i, max_i = _bisect_spectrum(mzs, mz_value, tol)
            im[y - 1, x - 1] = reduce_func(ints[min_i : max_i + 1])
    return im


def _get_ion_image(p, mz_value, tol, z, reduce_func, coordinates, start_idx):
    """Utility method used by mutlicore/thread image reader"""
    tol = abs(tol)
    im = np.zeros(
        (p.imzmldict["max count of pixels y"], p.imzmldict["max count of pixels x"])
    )
    for i, (x, y, z_) in enumerate(coordinates, start=start_idx):
        if z_ == 0:
            UserWarning(
                "z coordinate = 0 present, if you're getting blank images set getionimage(.., .., z=0)"
            )
        if z_ == z:
            mzs, ints = map(lambda x: np.asarray(x), p.get_spectrum(i))
            min_i, max_i = _bisect_spectrum(mzs, mz_value, tol)
            im[y - 1, x - 1] = reduce_func(ints[min_i : max_i + 1])
    return im


# keep alias to previously named functions
getionimage = get_ion_image


def browse(p):
    """
    Create a per-spectrum metadata browser for the parser.
    Usage::

        # get a list of the instrument configurations used in the first pixel
        instrument_configurations = browse(p).for_spectrum(0).get_ids("instrumentConfiguration")

    Currently, ``instrumentConfiguration``, ``dataProcessing`` and ``referenceableParamGroup`` are supported.

    For browsing all spectra iteratively, you should by all means use **ascending** indices. Doing otherwise can result
    in quadratic runtime. The following example shows how to retrieve all unique instrumentConfigurations used::

        browser = browse(p)
        all_config_ids = set()
        for i, _ in enumerate(p.coordinates):
            all_config_ids.update(browser.for_spectrum(i).get_ids("instrumentConfiguration"))

    This is a list of ids with which you can find the corresponding ``<instrumentConfiguration>`` tag in the xml tree.

    :param p: the parser
    :return: the browser
    """
    return _ImzMLMetaDataBrowser(p.root, p.filename, p.sl)


def _bisect_spectrum(mzs, mz_value, tol):
    ix_l, ix_u = bisect_left(mzs, mz_value - tol), bisect_right(mzs, mz_value + tol) - 1
    if ix_l == len(mzs):
        return len(mzs), len(mzs)
    if ix_u < 1:
        return 0, 0
    if ix_u == len(mzs):
        ix_u -= 1
    if mzs[ix_l] < (mz_value - tol):
        ix_l += 1
    if mzs[ix_u] > (mz_value + tol):
        ix_u -= 1
    return ix_l, ix_u


class _ImzMLMetaDataBrowser:
    def __init__(self, root, fn, sl):
        self._root = root
        self._sl = sl
        self._fn = fn
        self._iter, self._previous, self._list_elem = None, None, None
        self.iterparse = choose_iterparse()

    def for_spectrum(self, i):
        if self._previous is None or i <= self._previous:
            self._iter = self.iterparse(self._fn, events=("start", "end"))
        for event, s in self._iter:
            if s.tag == self._sl + "spectrumList" and event == "start":
                self._list_elem = s
            elif s.tag == self._sl + "spectrum" and event == "end":
                self._list_elem.remove(s)
                if s.attrib["index"] == str(i):
                    self._previous = i
                    return _SpectrumMetaDataBrowser(self._root, self._sl, s)


class _SpectrumMetaDataBrowser:
    def __init__(self, root, sl, spectrum):
        self._root = root
        self._sl = sl
        self._spectrum = spectrum

    def get_ids(self, element):
        param_methods = {
            param_group_elname: self._find_referenceable_param_groups,
            data_processing_elname: self._find_data_processing,
            instrument_confid_elname: self._find_instrument_configurations,
        }
        try:
            return param_methods[element]()
        except KeyError as e:
            raise ValueError("Unsupported element: " + str(element))

    def _find_referenceable_param_groups(self):
        param_group_refs = self._spectrum.findall(
            "%sreferenceableParamGroupRef" % self._sl
        )
        ids = map(lambda g: g.attrib["ref"], param_group_refs)
        return ids

    def _find_instrument_configurations(self):
        ids = None
        scan_list = self._spectrum.find("%sscanList" % self._sl)
        if scan_list:
            scans = scan_list.findall("%sscan[@instrumentConfigurationRef]" % self._sl)
            ids = map(lambda s: s.attrib["instrumentConfigurationRef"], scans)
        if not ids:
            run = self._root.find("%srun")
            try:
                return [run.attrib["defaultInstrumentConfigurationRef"]]
            except KeyError as _:
                return list()
        else:
            return ids

    def _find_data_processing(self):
        try:
            return self._spectrum.attrib["dataProcessingRef"]
        except KeyError as _:
            spectrum_list = self._root.find(
                "%srun/%sspectrumList" % tuple(2 * [self._sl])
            )
            try:
                return [spectrum_list.attrib["defaultDataProcessingRef"]]
            except KeyError as _:
                return []


class PortableSpectrumReader:
    """A pickle-able class for holding the minimal set of data required for reading, without holding any
    references to open files that wouldn't survive pickling"""

    def __init__(
        self,
        coordinates,
        mzPrecision,
        mzOffsets,
        mzLengths,
        intensityPrecision,
        intensityOffsets,
        intensityLengths,
        ibd_filename=None,
        imzmldict=None,
    ):
        self.coordinates = coordinates
        self.mzPrecision = mzPrecision
        self.mzOffsets = mzOffsets
        self.mzLengths = mzLengths
        self.intensityPrecision = intensityPrecision
        self.intensityOffsets = intensityOffsets
        self.intensityLengths = intensityLengths

        self.ibd_filename = ibd_filename
        self.imzmldict = imzmldict
        self._n_pixels = None

    @property
    def n_pixels(self):
        """Return number of pixels in the dataset"""
        if self._n_pixels is None:
            self._n_pixels = len(self.coordinates)
        return self._n_pixels

    def read_spectrum_from_file(self, file, index):
        """
        Reads the spectrum at specified index from the .ibd file.

        :param file:
            File or file-like object for the .ibd file
        :param index:
            Index of the desired spectrum in the .imzML file

        Output:

        mz_array: numpy.ndarray
            Sequence of m/z values representing the horizontal axis of the desired mass
            spectrum
        intensity_array: numpy.ndarray
            Sequence of intensity values corresponding to mz_array
        """
        file.seek(self.mzOffsets[index])
        mz_bytes = file.read(self.mzLengths[index] * SIZE_DICT[self.mzPrecision])
        file.seek(self.intensityOffsets[index])
        intensity_bytes = file.read(
            self.intensityLengths[index] * SIZE_DICT[self.intensityPrecision]
        )

        mz_array = np.frombuffer(mz_bytes, dtype=self.mzPrecision)
        intensity_array = np.frombuffer(intensity_bytes, dtype=self.intensityPrecision)

        return mz_array, intensity_array

    def get_spectrum(self, index):
        """Reads the spectrum at specified index from the .ibd file.

        :param index:
            Index of the desired spectrum in the .imzML file

        Output:

        mz_array: numpy.ndarray
            Sequence of m/z values representing the horizontal axis of the desired mass
            spectrum
        intensity_array: numpy.ndarray
            Sequence of intensity values corresponding to mz_array
        """
        mz_bytes, intensity_bytes = self._read_spectrum(index)
        mz_array = np.frombuffer(mz_bytes, dtype=self.mzPrecision)
        intensity_array = np.frombuffer(intensity_bytes, dtype=self.intensityPrecision)
        return mz_array, intensity_array

    def _read_spectrum(self, index):
        """Reads m/z array and intensity array of the spectrum at specified location
        from the binary file as a byte string. The string can be unpacked by the struct
        module. To get the arrays as numbers, use getspectrum

        :param index:
            Index of the desired spectrum in the .imzML file
        :rtype: Tuple[str, str]

        Output:

        mz_string:
            string where each character represents a byte of the mz array of the
            spectrum
        intensity_string:
            string where each character represents a byte of the intensity array of
            the spectrum
        """
        offsets = [self.mzOffsets[index], self.intensityOffsets[index]]
        lengths = [self.mzLengths[index], self.intensityLengths[index]]
        lengths[0] *= SIZE_DICT[self.mzPrecision]
        lengths[1] *= SIZE_DICT[self.intensityPrecision]
        with open(self.ibd_filename, "rb") as ibd_handle:
            ibd_handle.seek(offsets[0])
            mz_string = ibd_handle.read(lengths[0])
            ibd_handle.seek(offsets[1])
            intensity_string = ibd_handle.read(lengths[1])

        return mz_string, intensity_string
