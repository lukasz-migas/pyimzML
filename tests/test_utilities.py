from pyimzml.utilities import bisect_spectrum


class TestBisect:
    @staticmethod
    def test_bisect():
        mz_x = [100.0, 201.89, 201.99, 202.0, 202.01, 202.10000001, 400.0]
        test_mz = 202.0
        test_tol = 0.1
        ix_l, ix_u = bisect_spectrum(mz_x, test_mz, test_tol)
        assert ix_l == 2
        assert ix_u == 4
        assert ix_l <= ix_u
        assert mz_x[ix_l] >= test_mz - test_tol
        assert mz_x[ix_u] <= test_mz + test_tol
