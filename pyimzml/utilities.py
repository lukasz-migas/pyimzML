"""Utility functions"""
from bisect import bisect_left
from bisect import bisect_right


def chunks(item_list, n_items):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(item_list), n_items):
        yield item_list[i : i + n_items]


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


def bisect_spectrum(mzs, mz_value, tol):
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
