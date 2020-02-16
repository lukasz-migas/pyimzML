from __future__ import print_function

# Standard library imports
import os
import uuid
import hashlib
import sys
import getopt
from collections import namedtuple, OrderedDict, defaultdict

# Third-party imports
import numpy as np
from wheezy.template import Engine, CoreExtension, DictLoader

# Local imports
from pyimzml.compression import NoCompression, ZlibCompression
from pyimzml.template import IMZML_TEMPLATE


class _MaxlenDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        self.maxlen = kwargs.pop('maxlen', None)
        OrderedDict.__init__(self, *args, **kwargs)

    def __setitem__(self, key, value):
        if self.maxlen is not None and len(self) >= self.maxlen:
            self.popitem(0)  # pop oldest
        OrderedDict.__setitem__(self, key, value)


# todo: change named tuple to dict and parse xml template properly (i.e. remove hardcoding so parameters can be
#       optional)
_Spectrum = namedtuple(
    '_Spectrum',
    "coords mz_len mz_offset mz_enc_len int_len int_offset "
    "int_enc_len mz_min mz_max mz_base int_base int_tic userParams")


class ImzMLWriter(object):
    """
        Create an imzML+ibd file.

        :param output_filename:
            is used to make the base name by removing the extension (if any).
            two files will be made by adding ".ibd" and ".imzML" to the base name
        :param intensity_dtype:
            The numpy data type to use for saving intensity values
        :param mz_dtype:
            The numpy data type to use for saving mz array values
        :param mode:

            * "continuous" mode will save the first mz array only
            * "processed" mode save every mz array separately
            * "auto" mode writes only mz arrays that have not already been written
        :param intensity_compression:
            How to compress the intensity data before saving
            must be an instance of :class:`~pyimzml.compression.NoCompression` or :class:`~pyimzml.compression.ZlibCompression`
        :param mz_compression:
            How to compress the mz array data before saving
    """

    def __init__(self, output_filename,
                 mz_dtype=np.float64, intensity_dtype=np.float32, mode="auto", spec_type="centroid",
                 scan_direction="top_down", line_scan_direction="line_left_right", scan_pattern="one_way",
                 scan_type="horizontal_line",
                 mz_compression=NoCompression(), intensity_compression=NoCompression(),
                 polarity=None):

        self.mz_dtype = mz_dtype
        self.intensity_dtype = intensity_dtype
        self.mode = mode
        self.spec_type = spec_type
        self.mz_compression = mz_compression
        self.intensity_compression = intensity_compression
        self.run_id = os.path.splitext(output_filename)[0]
        self.filename = self.run_id + ".imzML"
        self.ibd_filename = self.run_id + ".ibd"
        self.xml = open(self.filename, "w")
        self.ibd = open(self.ibd_filename, "wb+")
        self.sha1 = hashlib.sha1()
        self.uuid = uuid.uuid4()

        self.scan_direction = scan_direction
        self.scan_pattern = scan_pattern
        self.scan_type = scan_type
        self.line_scan_direction = line_scan_direction

        self._write_ibd(self.uuid.bytes)

        self.wheezy_engine = Engine(loader=DictLoader({"imzml": IMZML_TEMPLATE}), extensions=[CoreExtension()])
        self.imzml_template = self.wheezy_engine.get_template("imzml")
        self.spectra = []
        self.first_mz = None
        self.hashes = defaultdict(list)  # mz_hash -> list of mz_data (disk location)
        self.lru_cache = _MaxlenDict(maxlen=10)  # mz_array (as tuple) -> mz_data (disk location)
        self._setPolarity(polarity)

    @staticmethod
    def _np_type_to_name(dtype):
        if dtype.__name__.startswith("float"):
            return "%s-bit float" % dtype.__name__[5:]
        elif dtype.__name__.startswith("int"):
            return "%s-bit integer" % dtype.__name__[3:]

    def _setPolarity(self, polarity):
        if polarity:
            if polarity.lower() in ["positive", "negative"]:
                self.polarity = polarity.lower()
            else:
                raise ValueError(
                    "value for polarity must be one of 'positive', 'negative'. Received: {}".format(polarity))
        else:
            self.polarity = ""

    def _write_xml(self):
        spectra = self.spectra
        mz_data_type = self._np_type_to_name(self.mz_dtype)
        int_data_type = self._np_type_to_name(self.intensity_dtype)
        obo_codes = {"32-bit integer": "1000519",
                     "16-bit float": "1000520",
                     "32-bit float": "1000521",
                     "64-bit integer": "1000522",
                     "64-bit float": "1000523",
                     "continuous": "1000030",
                     "processed": "1000031",
                     "zlib compression": "1000574",
                     "no compression": "1000576",
                     "line_bottom_up": "1000492",
                     "line_left_right": "1000491",
                     "line_right_left": "1000490",
                     "line_top_down": "1000493",
                     "bottom_up": "1000400",
                     "left_right": "1000402",
                     "right_left": "1000403",
                     "top_down": "1000401",
                     "meandering": "1000410",
                     "one_way": "1000411",
                     "random_access": "1000412",
                     "horizontal_line": "1000480",
                     "vertical_line": "1000481"}
        obo_names = {"line_bottom_up": "line scan bottom up",
                     "line_left_right": "line scan left right",
                     "line_right_left": "line scan right left",
                     "line_top_down": "line scan top down",
                     "bottom_up": "bottom up",
                     "left_right": "left right",
                     "right_left": "right left",
                     "top_down": "top down",
                     "meandering": "meandering",
                     "one_way": "one way",
                     "random_access": "random access",
                     "horizontal_line": "horizontal line scan",
                     "vertical_line": "vertical line scan"}

        uuid = ("{%s}" % self.uuid).upper()
        sha1sum = self.sha1.hexdigest().upper()
        run_id = self.run_id
        if self.mode == "auto":
            mode = "processed" if len(self.lru_cache) > 1 else "continuous"
        else:
            mode = self.mode
        spec_type = self.spec_type
        mz_compression = self.mz_compression.name
        int_compression = self.intensity_compression.name
        polarity = self.polarity
        scan_direction = self.scan_direction
        scan_pattern = self.scan_pattern
        scan_type = self.scan_type
        line_scan_direction = self.line_scan_direction

        self.xml.write(self.imzml_template.render(locals()))

    def _write_ibd(self, bytes):
        self.ibd.write(bytes)
        self.sha1.update(bytes)
        return len(bytes)

    def _encode_and_write(self, data, dtype=np.float32, compression=NoCompression()):
        data = np.asarray(data, dtype=dtype)
        offset = self.ibd.tell()
        bytes = data.tobytes()
        bytes = compression.compress(bytes)
        return offset, data.shape[0], self._write_ibd(bytes)

    def _read_mz(self, mz_offset, mz_len, mz_enc_len):
        """reads a mz array from the currently open ibd file"""
        self.ibd.seek(mz_offset)
        data = self.ibd.read(mz_enc_len)
        self.ibd.seek(0, 2)
        data = self.mz_compression.decompress(data)
        return tuple(np.fromstring(data, dtype=self.mz_dtype))

    def _get_previous_mz(self, mzs):
        """given an mz array, return the mz_data (disk location)
        if the mz array was not previously written, write to disk first"""
        mzs = tuple(mzs)  # must be hashable
        if mzs in self.lru_cache:
            return self.lru_cache[mzs]

        # mz not recognized ... check hash
        mz_hash = "%s-%s-%s" % (hash(mzs), sum(mzs), len(mzs))
        if mz_hash in self.hashes:
            for mz_data in self.hashes[mz_hash]:
                test_mz = self._read_mz(*mz_data)
                if mzs == test_mz:
                    self.lru_cache[test_mz] = mz_data
                    return mz_data
        # hash not recognized
        # must be a new mz array ... write it, add it to lru_cache and hashes
        mz_data = self._encode_and_write(mzs, self.mz_dtype, self.mz_compression)
        self.hashes[mz_hash].append(mz_data)
        self.lru_cache[mzs] = mz_data
        return mz_data

    def add_spectrum(self, mz_x, mz_y, coords, userParams=[]):
        """
        Add a mass spectrum to the file.

        :param mz_x:
            mz array
        :param mz_y:
            intensity array
        :param coords:

            * 2-tuple of x and y position OR
            * 3-tuple of x, y, and z position

            note some applications want coords to be 1-indexed
        """
        # must be rounded now to allow comparisons to later data
        # but don't waste CPU time in continuous mode since the data will not be used anyway
        if self.mode != "continuous" or self.first_mz is None:
            mz_x = self.mz_compression.rounding(mz_x)
        mz_y = self.intensity_compression.rounding(mz_y)

        if self.mode == "continuous":
            if self.first_mz is None:
                self.first_mz = self._encode_and_write(mz_x, self.mz_dtype, self.mz_compression)
            mz_data = self.first_mz
        elif self.mode == "processed":
            mz_data = self._encode_and_write(mz_x, self.mz_dtype, self.mz_compression)
        elif self.mode == "auto":
            mz_data = self._get_previous_mz(mz_x)
        else:
            raise TypeError("Unknown mode: %s" % self.mode)
        mz_offset, mz_len, mz_enc_len = mz_data

        int_offset, int_len, int_enc_len = self._encode_and_write(mz_y, self.intensity_dtype,
                                                                  self.intensity_compression)
        mz_min = np.min(mz_x)
        mz_max = np.max(mz_x)
        ix_max = np.argmax(mz_y)
        mz_base = mz_x[ix_max]
        int_base = mz_y[ix_max]
        int_tic = np.sum(mz_y)
        s = _Spectrum(coords, mz_len, mz_offset, mz_enc_len, int_len, int_offset, int_enc_len, mz_min, mz_max, mz_base,
                      int_base, int_tic, userParams)
        self.spectra.append(s)

    addSpectrum = add_spectrum

    def close(self):  # 'close' is a more common use for this
        """Writes the XML file and closes all files.

        Will be called automatically if ``with``-pattern is used.
        """
        self.finish()

    def finish(self):
        """alias of close()"""
        self.ibd.close()
        self._write_xml()
        self.xml.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_t, exc_v, trace):
        if exc_t is None:
            self.finish()
        else:
            self.ibd.close()
            self.xml.close()


def _main(argv):
    from pyimzml.ImzMLParser import ImzMLParser
    inputfile = ""
    outputfile = ""
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print("test.py -i <inputfile> -o <outputfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("test.py -i <inputfile> -o <outputfile>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    if inputfile == "":
        print("test.py -i <inputfile> -o <outputfile>")
        raise IOError("input file not specified")
    if outputfile == "":
        outputfile = inputfile + ".imzML"
    imzml = ImzMLParser(inputfile)
    spectra = []
    with ImzMLWriter(outputfile, mz_dtype=np.float32, intensity_dtype=np.float32) as writer:
        for i, coords in enumerate(imzml.coordinates):
            mzs, intensities = imzml.get_spectrum(i)
            writer.add_spectrum(mzs, intensities, coords)
            spectra.append((mzs, intensities, coords))

    imzml = ImzMLParser(outputfile)
    spectra2 = []
    for i, coords in enumerate(imzml.coordinates):
        mzs, intensities = imzml.get_spectrum(i)
        spectra2.append((mzs, intensities, coords))

    print(spectra[0] == spectra2[0])


if __name__ == "__main__":
    _main(sys.argv[1:])
