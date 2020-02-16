"""Test pyimzml.imzwriter"""
# Standard library imports
import os

import numpy as np

# Third-party imports
import pytest
from numpy.testing import assert_array_almost_equal

from pyimzml.ImzMLParser import ImzMLParser

# Local imports
from pyimzml.ImzMLWriter import ImzMLWriter
from pyimzml.compression import NoCompression, ZlibCompression


class TestImzMLWriter:
    @staticmethod
    def test_writer_single_pixel(get_temp_path):
        mz_x = np.linspace(100, 1000, 20)
        mz_y = np.random.rand(mz_x.shape[0])
        coordinates = [1, 1, 1]

        output_filename = os.path.join(get_temp_path, "test.imzML")
        with ImzMLWriter(output_filename, mode="processed") as imzml:
            imzml.add_spectrum(mz_x, mz_y, coords=coordinates)

        with ImzMLParser(output_filename) as parser:
            _mz_x, _mz_y = parser.get_spectrum(0)
            assert_array_almost_equal(_mz_x, mz_x, 4)
            assert_array_almost_equal(_mz_y, mz_y, 4)
            assert parser.n_pixels == 1

    @staticmethod
    @pytest.mark.parametrize(
        "compression", (NoCompression(), ZlibCompression(), None, "None", "zlib")
    )
    def test_writer_with_compression(get_temp_path, compression):
        mz_x = np.linspace(100, 1000, 20)
        mz_y = np.random.rand(mz_x.shape[0])
        coordinates = [1, 1, 1]

        output_filename = os.path.join(get_temp_path, "test.imzML")
        with ImzMLWriter(
            output_filename,
            mode="processed",
            mz_compression=compression,
            intensity_compression=compression,
        ) as imzml:
            imzml.add_spectrum(mz_x, mz_y, coords=coordinates)

    @staticmethod
    @pytest.mark.parametrize("round_digits", (None, 4))
    def test_writer_zlib_compression_round(get_temp_path, round_digits):
        mz_x = np.linspace(100, 1000, 20)
        mz_y = np.random.rand(mz_x.shape[0])
        coordinates = [1, 1, 1]

        output_filename = os.path.join(get_temp_path, "test.imzML")
        compression = ZlibCompression(round_digits)
        with ImzMLWriter(
            output_filename,
            mode="processed",
            mz_compression=compression,
            intensity_compression=compression,
        ) as imzml:
            imzml.add_spectrum(mz_x, mz_y, coords=coordinates)
