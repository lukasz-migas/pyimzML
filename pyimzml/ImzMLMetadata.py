"""Metadata"""
# Local imports
from pyimzml.utilities import choose_iterparse

PARAM_GROUP_ELEMENT_NAME = "referenceableParamGroup"
PARAM_DATA_PROCESSING_ELEMENT_NAME = "dataProcessing"
PARAM_INSTRUMENT_CONFIG_ID_ELEMENT_NAME = "instrumentConfiguration"


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
            PARAM_GROUP_ELEMENT_NAME: self._find_referenceable_param_groups,
            PARAM_DATA_PROCESSING_ELEMENT_NAME: self._find_data_processing,
            PARAM_INSTRUMENT_CONFIG_ID_ELEMENT_NAME: self._find_instrument_configurations,
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
        if scan_list is not None:
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
