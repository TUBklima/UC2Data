import xarray
import numpy
import sys
import enum
import csv

is_win = sys.platform in ['win32', 'win64']
if not is_win:
    from cfchecker import cfchecks

data_content_file = "data_content.txt"
variables_file = "variables.txt"


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
            print(str(CheckRules[ikey].result) + ' (' + ikey + ')')


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

    def __init__(self, result: ResultCode, message):
        self.result = result
        self.message = message

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

        if type(self.func) == bool:
            if self.func:
                self.result = CheckResult(ResultCode.OK, "Test passed.", self.group)
            else:
                self.result = CheckResult(self.severity, self.msg, self.group)
        else:
            if self.func(*self.param):
                self.result = CheckResult(ResultCode.OK, "Test passed.", self.group)
            else:
                self.result = CheckResult(self.severity, self.msg, self.group)
        return self.result.result != ResultCode.ERROR

    def __str__(self):
        return self.result


def type_okay(var, allowed_types):
    return type(var) in allowed_types


def data_content_okay(value):
    with open(data_content_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            if value == row[1]:
                return True

    with open(variables_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            if value == row[3]:
                return True


def check_glob_attr(data, attrname, must_exist, allowed_types, exact_value):
    exists = attrname in data.attrs.keys()
    if not exists:
        if must_exist:
            return CheckResult(ResultCode.ERROR, "Required global attribute '" + attrname + "' not found.")
        else:
            return CheckResult(ResultCode.OK, "Global attribute '" + attrname + "' not found.")

    if len(allowed_types) > 0:
        if not type(data.attrs[attrname]) in allowed_types:
            return CheckResult(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong type. Should be " + \
                "one of the following: " + str(allowed_types))

    if len(exact_value) > 0:
        if data.attrs[attrname] != exact_value:
            return CheckResult(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong value. Should be " + \
                str(exact_value))

    return CheckResult(ResultCode.OK, "Test passed.")


def get_check_rules(data):
    CheckRules = dict()
    CheckRules["title_exists"] = Check(True, UC2Data.has_glob_attr, (data, "title"),
                                       "Global attribute 'title' expected. Not found.", ResultCode.ERROR,
                                       "obligatory global attributes")
    CheckRules["title_is_string"] = Check("title_exists", type_okay, (data.attrs["title"], [str]),
                                          "Global attribute 'title' must be string.", ResultCode.ERROR,
                                          "obligatory global attributes")
    CheckRules["data_content_exists"] = Check(True, UC2Data.has_glob_attr, (data, "data_content"),
                                              "Global attribute 'data_content' expected. Not found.", ResultCode.ERROR,
                                              "obligatory global attributes")
    CheckRules["data_content_okay"] = Check("data_content_exists", data_content_okay, (data.attrs["data_content"],),
                                            "Global attribute 'data_content' has unsupported value.", ResultCode.ERROR,
                                            "obligatory global attributes")
    CheckRules["source_exists"] = Check(True, UC2Data.has_glob_attr, (data, "source"),
                                        "Global attribute 'source' expected. Not found.", ResultCode.ERROR,
                                        "obligatory global attributes")
    CheckRules["source_is_string"] = Check("source_exists", type_okay, (data.attrs["source"], [str]),
                                           "Global attribute 'source' must be string.", ResultCode.ERROR,
                                           "obligatory global attributes")
    CheckRules["version_exists"] = Check(True, UC2Data.has_glob_attr, (data, "version"),
                                         "Global attribute 'version' expected. Not found.", ResultCode.ERROR,
                                         "obligatory global attributes")
    CheckRules["version_is_int"] = Check("version_exists", type_okay, (data.attrs["version"], [int, numpy.int16]),
                                           "Global attribute 'version' must be integer.", ResultCode.ERROR,
                                           "obligatory global attributes")
    CheckRules["version_range"] = Check("version_is_int", data.attrs["version"] in range(1,1000), None,
                                        "Global attribute 'version' out of allowed range (1,999)", ResultCode.ERROR,
                                        "obligatory global attributes")
    CheckRules["Conventions_exists"] = Check(True, UC2Data.has_glob_attr, (data, "Conventions"),
                                         "Global attribute 'Conventions' expected. Not found.", ResultCode.ERROR,
                                         "obligatory global attributes")
    CheckRules["Conventions_okay"] = Check("Conventions_exists", data.attrs["Conventions"] == "CF-1.7", None,
                                             "Global attribute 'Conventions' must be set to 'CF-1.7'.", ResultCode.ERROR,
                                             "obligatory global attributes")
    return CheckRules
