import xarray
import sys
from warnings import warn
import numpy
import datetime
import enum

is_win = sys.platform in ['win32', 'win64']
if not is_win:
    from cfchecker import cfchecks

class UC2Data(xarray.Dataset):

    def __init__(self, path):
        self.path = path
        super().__init__()
        tmp = xarray.open_dataset(self.path, decode_cf=False)
        self.update(tmp, inplace=True)
        self.attrs = tmp.attrs

    def has_var(self, varname):
        return varname in self.variables.keys()

    def has_dim(self, dimname):
        return dimname in self.dims.keys()

    def has_glob_attr(self, attrname):
        return attrname in self.attrs.keys()

    def uc2_check(self):

        CheckRules = get_check_rules(self)


        for ikey in CheckRules.keys():
            CheckRules[ikey].check(CheckRules)
            print(CheckRules[ikey].result)

class Variable:

    def __init__(self, variables, name):
        self.var = variables[name]
        self.name = name

    def has_attr(self, attrname):
        return attrname in self.var.attrs.keys()


class ResultCode(enum.Enum):
    OK = 1
    WARNING = 2
    ERROR = 3


class CheckResult:

    def __init__(self, result: ResultCode, message, group):
        self.result = result
        self.message = message
        self.group = group

    def __str__(self):
        return self.message

class Check:

    def __init__(self, do_check, func, param, msg, severity, group):
        self.do_check = do_check
        self.func = func
        self.param = param
        self.msg = msg
        self.severity = severity
        self.group = group
        self.result = None

    def check(self, other_checks):

        if type(self.do_check) == bool:
            do_check = self.do_check
        else:
            do_check = other_checks[self.do_check].check(other_checks)

        if not do_check:
            return False

        if self.func(*self.param):
            self.result = CheckResult(ResultCode.OK, "Test passed.", self.group)
        else:
            self.result = CheckResult(self.severity, self.msg, self.group)
        return self.result.result != ResultCode.ERROR

    def __str__(self):
        return self.result


def type_okay(var, allowed_types):
    return type(var) in allowed_types


def get_check_rules(data):
    CheckRules = dict()
    CheckRules["title_exists"] = Check(True, UC2Data.has_glob_attr, (data, "title"),
                                        "Global attribute 'title' expected. Not found.", ResultCode.ERROR,
                                        "obligatory global attributes")
    CheckRules["title_is_string"] = Check("title_exists", type_okay, (data.attrs["title"], [str]),
                                        "Global attribute 'title' must be string.", ResultCode.ERROR,
                                        "obligatory global attributes")
    return CheckRules