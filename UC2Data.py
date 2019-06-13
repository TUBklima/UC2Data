import xarray
import numpy
import sys
import enum
import csv
import re

is_win = sys.platform in ['win32', 'win64']
if not is_win:
    from cfchecker import cfchecks

data_content_file = "data_content.txt"
variables_file = "variables.txt"
institutions_file = "institutions.txt"


class UC2Data(xarray.Dataset):

    def __init__(self, path):
        self.path = path
        super().__init__()
        tmp = xarray.open_dataset(self.path, decode_cf=False)
        self.update(tmp, inplace=True)
        self.attrs = tmp.attrs

    def uc2_check(self):

        allowed_data_contents = get_allowed_data_contents(data_content_file, variables_file)
        allowed_institutions = get_allowed_institutions(institutions_file)
        allowed_licences = get_allowed_licences()

        ###
        # Check global attributes
        ###
        result = dict()
        result["title"] = self.check_glob_attr("title", True, str)
        result["data_content"] = self.check_glob_attr("data_content", True, str, allowed_values=allowed_data_contents) # TODO: Redo this test when variable is checked
        result["source"] = self.check_glob_attr("source", True, str)
        result["version"] = self.check_glob_attr("version", True, [int, numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.int64], allowed_values=list(range(1,1000))) # TODO: This ist going to be checked in DMS
        result["Conventions"] = self.check_glob_attr("Conventions", True, str, allowed_values=["CF-1.7"])
        result["dependencies"] = self.check_glob_attr("dependencies", True, str) # TODO: This is going to be checked by DMS
        result["history"] = self.check_glob_attr("history", True, str)
        result["institution"] = self.check_glob_attr("institution", True, str) # TODO: Check with acronym
        result["acronym"] = self.check_glob_attr("acronym", True, str)  # TODO: Check with institution
        result["author"] = self.check_glob_attr("author", True, str)  # TODO: regex
        result["contact_person"] = self.check_glob_attr("contact_person", True, str)  # TODO: regex
        result["references"] = self.check_glob_attr("references", True, str)
        result["comment"] = self.check_glob_attr("comment", True, str)
        result["keywords"] = self.check_glob_attr("keywords", True, str)
        result["licence"] = self.check_glob_attr("licence", True, str, allowed_values=allowed_licences)

        return result
        # CheckRules = get_check_rules(self)
        #
        # for ikey in CheckRules.keys():
        #     CheckRules[ikey].check(CheckRules)
        #     print(str(CheckRules[ikey].result) + ' (' + ikey + ')')

    def check_glob_attr(self, attrname, must_exist, allowed_types, allowed_values=[], regex=None):
        exists = attrname in self.attrs.keys()
        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR, "Required global attribute '" + attrname + "' not found.")
            else:
                return CheckResult(ResultCode.OK, "Global attribute '" + attrname + "' not found.")

        if not type(allowed_types) == list:
            allowed_types = [allowed_types]
        if len(allowed_types) > 0:
            if not type(self.attrs[attrname]) in allowed_types:
                return CheckResult(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong type. Should be " + \
                                   "one of the following: " + str(allowed_types))

        if not type(allowed_values) == list:
            allowed_values = [allowed_values]
        if len(allowed_values) > 0:
            if not self.attrs[attrname] in allowed_values:
                if len(allowed_values) == 1:
                    return CheckResult(ResultCode.ERROR,
                                       "Global attribute '" + attrname + "' has wrong value. Should be " + \
                                       str(allowed_values[0]))
                else:
                    return CheckResult(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong value")

        if not regex is None:
            pass # TODO: check regular expression (e.g. author)

        return CheckResult(ResultCode.OK, "Test passed.")


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


def get_allowed_data_contents(data_content_file, variables_file):
    out = []
    with open(data_content_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            out.append(row[1])

    with open(variables_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            out.append(row[3])

    return out


def get_allowed_institutions(institutions_file):
    pass


def get_allowed_licences():
    return ["[UC]2 MOSAIK Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
           "[UC]2 3DO Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
           "[UC]2 KliMoPrax Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
           "[UC]2 UseUClim Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
           "[UC]2 Restriced Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
           "[UC]2 Research Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
           "[UC]2 Open Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
           ]
