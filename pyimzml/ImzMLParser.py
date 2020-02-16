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

# Feb 2020 - modified by Lukasz G. Migas

# Standard library imports
import math
import time
import re
from pathlib import Path
from warnings import warn
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import wait
import logging

# Third-party imports
import numpy as np

# Local imports
from pyimzml.utilities import chunks
from pyimzml.utilities import choose_iterparse
from pyimzml.utilities import bisect_spectrum
from pyimzml.utilities import format_time


PRECISION_DICT = {
    "32-bit float": "f",
    "64-bit float": "d",
    "32-bit integer": "i",
    "64-bit integer": "l",
}
SIZE_DICT = {"f": 4, "d": 8, "i": 4, "l": 8}

LOGGER = logging.getLogger(__name__)


class ImzMLParserBase:
    def __init__(self):

        self.filename = None
        self.ibd_filename = None
        self.coordinates = []
        self.mzOffsets = []
        self.intensityOffsets = []
        self.mzLengths = []
        self.intensityLengths = []
        self.mzPrecision = None
        self.intensityPrecision = None
        self.imzmldict = dict()

        self._idx = 0
        self._n_pixels = None
        self._as_threads = None
        self._pool_size = None
        self.executor = None

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
        t_start = time.time()
        mz_bytes, intensity_bytes = self._read_spectrum(index)
        mz_array = np.frombuffer(mz_bytes, dtype=self.mzPrecision)
        intensity_array = np.frombuffer(intensity_bytes, dtype=self.intensityPrecision)
        LOGGER.debug(f"Retrieved spectrum in {format_time(time.time()-t_start)}")
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
        lengths[0] *= SIZE_DICT[self.mzPrecision]
        lengths[1] *= SIZE_DICT[self.intensityPrecision]
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
        t_start = time.time()
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
        LOGGER.debug(f"Retrieved ion image in {format_time(time.time() - t_start)}")
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


class ImzMLParser(ImzMLParserBase):
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
        ImzMLParserBase.__init__(self)
        t_start = time.time()
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

        self.root = None
        self.mzGroupId = self.intGroupId = None
        self._parse_lib = parse_lib

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
        LOGGER.debug(f"Initilized parser in {format_time(time.time()-t_start)}")

    def __repr__(self):
        return f"ImzMLParser<filename: {self.filename}; no. pixels: {self.n_pixels}>"

    def __enter__(self):
        """Enter - system method used `with ... as`"""
        return self

    def __exit__(self, exc_t, exc_v, trace):
        """Exit - system method used `with ... as `"""

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
        LOGGER.debug("Setup metadata")

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
        LOGGER.debug("Fixed offsets")

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
        LOGGER.debug("Setup precision")

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
    t_start = time.time()
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
            min_i, max_i = bisect_spectrum(mzs, mz_value, tol)
            im[y - 1, x - 1] = reduce_func(ints[min_i : max_i + 1])
    LOGGER.debug(f"Retrieved ion image in {format_time(time.time() - t_start)}")
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
            min_i, max_i = bisect_spectrum(mzs, mz_value, tol)
            im[y - 1, x - 1] = reduce_func(ints[min_i : max_i + 1])
    return im


# keep alias to previously named functions
getionimage = get_ion_image


class PortableSpectrumReader(ImzMLParserBase):
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
        ImzMLParserBase.__init__(self)
        self.coordinates = coordinates
        self.mzPrecision = mzPrecision
        self.mzOffsets = mzOffsets
        self.mzLengths = mzLengths
        self.intensityPrecision = intensityPrecision
        self.intensityOffsets = intensityOffsets
        self.intensityLengths = intensityLengths

        self.ibd_filename = ibd_filename
        self.imzmldict = imzmldict

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
