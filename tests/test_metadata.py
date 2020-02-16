"""Test pyimzml.metadata"""
# Standard library imports
from pathlib import Path


# Third-party imports
import pytest

# Local imports
from pyimzml.ImzMLParser import ImzMLParser
from pyimzml.ImzMLMetadata import browse


# Example files from https://ms-imaging.org/wp/imzml/example-files-test/
CONTINUOUS_IMZML_PATH = str(Path(__file__).parent / "data/Example_Continuous.imzML")
CONTINUOUS_IBD_PATH = str(Path(__file__).parent / "data/Example_Continuous.ibd")
PROCESSED_IMZML_PATH = str(Path(__file__).parent / "data/Example_Processed.imzML")
PROCESSED_IBD_PATH = str(Path(__file__).parent / "data/Example_Processed.ibd")


@pytest.mark.parametrize("data_path", (CONTINUOUS_IMZML_PATH, PROCESSED_IMZML_PATH))
@pytest.mark.parametrize("parse_lib", ("lxml", "ElementTree"))
@pytest.mark.parametrize(
    "item_ids", ("referenceableParamGroup", "dataProcessing", "instrumentConfiguration")
)
def test_browse(data_path, parse_lib, item_ids):
    parser = ImzMLParser(data_path, parse_lib=parse_lib)
    browser = browse(parser)
    assert browser

    all_item_ids = set()
    for i in range(parser.n_pixels):
        all_item_ids.update(browser.for_spectrum(i).get_ids(item_ids))

    assert len(all_item_ids) != 0
