"""Config file for pytest"""
# Standard library imports

# Third-party imports
import pytest

print(">Sd")
@pytest.fixture(scope="session", autouse=True)
def get_temp_path(tmpdir_factory):
    output_dir = str(tmpdir_factory.mktemp("output"))
    return output_dir
