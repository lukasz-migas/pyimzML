import zlib


class NoCompression(object):
    """No compression"""

    name = "no compression"

    def __init__(self):
        pass

    @staticmethod
    def rounding(data):
        return data

    @staticmethod
    def compress(bytes):
        return bytes

    @staticmethod
    def decompress(bytes):
        return bytes


class ZlibCompression(object):
    """Zlib compression with optional rounding of values. Rounding helps the compression, but is lossy.

    :param round_amt:
        Number of digits after comma. None means no rounding.
    """
    name = "zlib compression"

    def __init__(self, round_amt=None):
        self.round_amt = round_amt

    def rounding(self, data):
        if self.round_amt is not None:
            return [round(x, self.round_amt) for x in data]
        return data

    @staticmethod
    def compress(bytes):
        return zlib.compress(bytes)

    @staticmethod
    def decompress(bytes):
        return zlib.decompress(bytes)
