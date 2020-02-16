"""Test pyimzml.imzmlparser.py"""
# Standard library imports
from pathlib import Path
import pickle

# Third-party imports
import pytest
import numpy as np
from numpy.testing import assert_equal

# Local imports
from pyimzml.ImzMLParser import ImzMLParser
from pyimzml.ImzMLParser import _bisect_spectrum

# Example files from https://ms-imaging.org/wp/imzml/example-files-test/
CONTINUOUS_IMZML_PATH = str(Path(__file__).parent / 'data/Example_Continuous.imzML')
CONTINUOUS_IBD_PATH = str(Path(__file__).parent / 'data/Example_Continuous.ibd')
PROCESSED_IMZML_PATH = str(Path(__file__).parent / 'data/Example_Processed.imzML')
PROCESSED_IBD_PATH = str(Path(__file__).parent / 'data/Example_Processed.ibd')


class TestImzMLParser:
    @staticmethod
    @pytest.mark.parametrize("data_path", (CONTINUOUS_IMZML_PATH, PROCESSED_IMZML_PATH))
    @pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
    def test_parser_init_paths(data_path, parse_lib):
        parser = ImzMLParser(data_path, parse_lib=parse_lib)
        assert len(parser.coordinates) == 9
        assert parser.n_pixels == 9

        mz_x, mz_y = parser.get_spectrum(0)
        assert len(mz_x) == len(mz_y)
        assert len(mz_x) > 0
        assert len(mz_y) > 0

        mz_x, mz_y = parser.get_spectrum(4)
        assert len(mz_x) == len(mz_y)
        assert len(mz_x) == 8399
        assert len(mz_y) == 8399
        assert np.all(mz_x > 100.0)
        assert np.all(mz_x < 800.0)
        assert np.all(mz_y >= 0.0)
        assert np.all(mz_y < 3.0)

    @staticmethod
    @pytest.mark.parametrize("data_path", (CONTINUOUS_IMZML_PATH, PROCESSED_IMZML_PATH))
    @pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
    def test_parser_init_paths_as_with(data_path, parse_lib):
        with ImzMLParser(data_path, parse_lib=parse_lib) as parser:
            assert len(parser.coordinates) == 9
            assert parser.n_pixels == 9

            mz_x, mz_y = parser.get_spectrum(0)
            assert len(mz_x) == len(mz_y)
            assert len(mz_x) > 0
            assert len(mz_y) > 0

    @staticmethod
    @pytest.mark.parametrize("imzml_path, ibd_path", ([CONTINUOUS_IMZML_PATH, CONTINUOUS_IBD_PATH],
                                                      [PROCESSED_IMZML_PATH, PROCESSED_IBD_PATH]))
    @pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
    def test_parser_init_ibd_as_file(imzml_path, ibd_path, parse_lib):
        with open(ibd_path, "rb") as ibd_file:
            with ImzMLParser(imzml_path, parse_lib=parse_lib, ibd_file=ibd_file) as parser:
                assert len(parser.coordinates) == 9
                assert parser.n_pixels == 9

                mz_x, mz_y = parser.get_spectrum(0)
                assert len(mz_x) == len(mz_y)
                assert len(mz_x) > 0
                assert len(mz_y) > 0

    @staticmethod
    @pytest.mark.parametrize("data_path", (CONTINUOUS_IMZML_PATH, PROCESSED_IMZML_PATH))
    @pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
    def test_parser_get_spectrum(data_path, parse_lib):
        parser = ImzMLParser(data_path, parse_lib=parse_lib)
        for px in range(parser.n_pixels):
            mz_x, mz_y = parser.get_spectrum(px)
            assert len(mz_x) == len(mz_y)
            assert len(mz_x) > 0
            assert len(mz_y) > 0

    @staticmethod
    @pytest.mark.parametrize("data_path", (CONTINUOUS_IMZML_PATH, PROCESSED_IMZML_PATH))
    @pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
    def test_parser_iter(data_path, parse_lib):
        parser = ImzMLParser(data_path, parse_lib=parse_lib)

        count = 0
        for px, (mz_x, mz_y) in enumerate(parser):
            _mz_x, _mz_y = parser.get_spectrum(px)
            assert len(mz_x) == len(mz_y)
            assert len(mz_x) == len(_mz_x)
            assert len(mz_y) == len(_mz_y)
            assert_equal(_mz_x, mz_x)
            assert_equal(_mz_y, mz_y)
            count += 1

        assert count == parser.n_pixels

    @staticmethod
    @pytest.mark.parametrize("data_path", (CONTINUOUS_IMZML_PATH, PROCESSED_IMZML_PATH))
    @pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
    def test_parser_get_image(data_path, parse_lib):
        parser = ImzMLParser(data_path, parse_lib=parse_lib)
        im = parser.get_ion_image(500, 100)
        assert im.sum() > 0

    @staticmethod
    @pytest.mark.parametrize("data_path", (CONTINUOUS_IMZML_PATH, PROCESSED_IMZML_PATH))
    @pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
    def test_parser_get_physical_coordinates(data_path, parse_lib):
        parser = ImzMLParser(data_path, parse_lib=parse_lib)

        x, y = parser.get_physical_coordinates(0)
        assert x == 100.
        assert y == 100.


class TestPortableSpectrumReader:

    @staticmethod
    @pytest.mark.parametrize("imzml_path, ibd_path", ([CONTINUOUS_IMZML_PATH, CONTINUOUS_IBD_PATH],
                                                      [PROCESSED_IMZML_PATH, PROCESSED_IBD_PATH]))
    @pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
    def test_portable_init(imzml_path, ibd_path, parse_lib):
        # get normal parser
        parser = ImzMLParser(imzml_path, parse_lib=parse_lib)

        # get detached parser and get handle of the portable reader
        detached_parser = ImzMLParser(imzml_path, parse_lib=parse_lib)
        portable_reader = detached_parser.portable_spectrum_reader()

        # pickle and unpickle to ensure it survives for its intended use case
        portable_reader = pickle.loads(pickle.dumps(portable_reader))

        with open(ibd_path, "rb") as ibd_file:
            for idx in range(parser.n_pixels):
                mz_x, mz_y = parser.get_spectrum(idx)
                _mz_x, _mz_y = portable_reader.read_spectrum_from_file(ibd_file, idx)
                assert np.all(mz_x == _mz_x)
                assert np.all(mz_y == _mz_y)


class TestBisect:
    @staticmethod
    def test_bisect():
        mz_x = [100., 201.89, 201.99, 202.0, 202.01, 202.10000001, 400.]
        test_mz = 202.0
        test_tol = 0.1
        ix_l, ix_u = _bisect_spectrum(mz_x, test_mz, test_tol)
        assert ix_l == 2
        assert ix_u == 4
        assert ix_l <= ix_u
        assert mz_x[ix_l] >= test_mz - test_tol
        assert mz_x[ix_u] <= test_mz + test_tol
