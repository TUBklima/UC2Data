import unittest
import os
from UC2Data import *


class TestCheckResult(unittest.TestCase):

    file_dir = "test/test_files/"
    def test_add(self):
        a = CheckResult()

        a.add(ResultCode.OK)
        self.assertTrue(a)

        a.add(ResultCode.ERROR, "blöd")
        self.assertFalse(a)

    def test_ok_files_pass(self):
        files = ["grid", "timeSeries", "timeSeriesProfile", "trajectory"]
        for fn in files:
            fn = self.file_dir + fn + ".nc"

            data = UC2Data(fn)
            data.uc2_check()
            self.assertTrue(data.check_result)
            self.assertFalse(data.check_result.contains_warnings())

    def test_nonsense_fails(self):
        fn = self.file_dir + "nonsense.nc"
        data = UC2Data(fn)
        data.uc2_check()
        self.assertFalse(data.check_result)

if __name__ == '__main__':
    unittest.main()