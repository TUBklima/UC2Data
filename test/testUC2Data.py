import unittest
from UC2Data import *


class TestCheckResult(unittest.TestCase):

    def test_add(self):
        a = CheckResult()

        a.add(ResultCode.OK)
        self.assertTrue(a)

        a.add(ResultCode.ERROR, "blöd")
        self.assertFalse(a)


if __name__ == '__main__':
    unittest.main()