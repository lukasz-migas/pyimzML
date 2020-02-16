"""Test pyimzml.compression.py"""
# Third-party imports
import numpy as np
from numpy.testing import assert_array_equal

# Local imports
from pyimzml.compression import NoCompression
from pyimzml.compression import ZlibCompression


class TestNoCompression:
    @staticmethod
    def test_compress_float():
        count = 100
        array = np.random.rand(count)
        compressor = NoCompression()
        compressed = compressor.compress(array.tobytes())
        decompressed = np.frombuffer(compressor.decompress(compressed), count=count)
        assert_array_equal(array, decompressed)

    @staticmethod
    def test_compress_int():
        count = 100
        array = np.random.randint(0, 1000, count, dtype=np.int32)
        compressor = NoCompression()
        compressed = compressor.compress(array.tobytes())
        decompressed = np.frombuffer(
            compressor.decompress(compressed), count=count, dtype=np.int32
        )
        assert_array_equal(array, decompressed)


class TestZlibCompression:
    @staticmethod
    def test_compress_float():
        count = 100
        array = np.random.rand(count)
        compressor = ZlibCompression()
        compressed = compressor.compress(array.tobytes())
        decompressed = np.frombuffer(compressor.decompress(compressed), count=count)
        assert_array_equal(array, decompressed)

    @staticmethod
    def test_compress_int():
        count = 100
        array = np.random.randint(0, 1000, count, dtype=np.int32)
        compressor = ZlibCompression()
        compressed = compressor.compress(array.tobytes())
        decompressed = np.frombuffer(
            compressor.decompress(compressed), count=count, dtype=np.int32
        )
        assert_array_equal(array, decompressed)
