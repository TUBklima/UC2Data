import unittest
from uc2data.Dataset import *
from pathlib import Path


class TestCheckResult(unittest.TestCase):

    file_dir = Path(__file__).parent / "test_files"

    def test_add(self):
        a = CheckResult()

        a.add(ResultCode.OK)
        self.assertTrue(a)

        a.add(ResultCode.ERROR, "blöd")
        self.assertFalse(a)

    def test_ok_files_pass(self):
        files = ["grid", "timeSeries", "timeSeriesProfile", "trajectory"]
        for fn in files:
            fn = self.file_dir / (fn + ".nc")

            data = Dataset(fn)
            data.uc2_check()
            self.assertTrue(data.check_result)
            self.assertTrue(len(data.check_result.errors) == 0)
            warn = data.check_result.warnings
            # TODO: Check that no warnings (currently there is one in origin lon/lat vs x/y in timeSeries.nc

            self.assertTrue(type(data.filename) == str)

    def test_nonsense_fails(self):
        fn = self.file_dir / "nonsense.nc"
        data = Dataset(fn)
        data.uc2_check()
        self.assertFalse(data.check_result)

if __name__ == '__main__':
    unittest.main()