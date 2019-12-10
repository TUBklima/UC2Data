from numpy.core.defchararray import add as str_add
from collections import OrderedDict
import enum


class ResultCode(enum.Enum):

    """
    This object represents whether a check results in OK, warning or error.

    """

    OK = 1
    WARNING = 2
    ERROR = 3


class ResultItem:

    """
    This object combines the ResultCode class with a message string.

    Attributes
    ----------
    result : ResultCode
        Whether the check was OK or resulted in a warning or error
    message : str
        A message for the user

    """

    def __init__(self, result: ResultCode = ResultCode.OK, message: str = ""):

        """
        Creates a ResultItem object

        Parameters
        ----------
        result : ResultCode
            Whether the check was OK or resulted in a warning or error. Default: OK
        message : str
            A message for the user. Default: "Test passed."
        """

        self.result = result
        if result == ResultCode.OK:
            if message != "":
                raise Exception("cannot handle user message for ResultCode.OK")
            self.message = "Test passed."
        else:
            self.message = message

    def __bool__(self):

        """
        Checks whether the test passed

        Returns
        -------
        True if ResultCode is OK or WARNING, False otherwise

        """

        return self.result != ResultCode.ERROR


class CheckResult(OrderedDict):

    """
    This object collects ResultItems in a structured way.

    The ResultItems are structured in a dict-like structure.
    The tags can be nested. TODO: Keep on documenting here!!!

    Attributes
    ----------
    result : ResultCode
        Whether the check was OK or resulted in a warning or error
    message : str
        A message for the user

    """

    def __init__(self, *args, **kwargs):
        self.result = list()
        if args or kwargs:
            self.add(*args, **kwargs)

    def __getitem__(self, item):
        if item not in super().keys():
            self[item] = CheckResult()
        return super().__getitem__(item)

    def __bool__(self):
        if len(self.result) == 0:
            # an empty thing is not True (otherwise you couldnt check for bool(result["var"]) if "var" is not in result
            ok = len(self.keys()) != 0
        else:
            ok = all(i_ok for i_ok in self.result)

        if not ok:
            return ok

        for val in self.values():
            if not val:
                return False

        return True

    def contains_warnings(self):
        if len(self.result) == 0:
            has_warn = False
        else:
            has_warn = any([ir.result == ResultCode.WARNING for ir in self.result])

        if has_warn:
            return has_warn

        for val in self.values():
            if val.contains_warnings():
                return True

        return False

    def add(self, result, message=""):

        if isinstance(result, ResultCode):
            other = ResultItem(result, message)
        else:
            other = result

        if isinstance(other, ResultItem):
            if other.result == ResultCode.OK:
                if len(self.result) == 0:
                    self.result.append(other)
                else:
                    return  # There are ERRORs in result. => not OK, If OK => stays OK
            else:
                for i in self.result:
                    if i.result == ResultCode.OK:
                        self.result.remove(i)  # remove OK from result because ERROR is added.
                self.result.append(other)

        elif isinstance(other, CheckResult):
            for i in other.result:
                self.add(i)
            for key, value in other.items():
                self[key].add(value)
        else:
            raise Exception("unexpected type of other")

    def __repr__(self):

        out = list()
        for i in self.result:
            out.append(i.message + " (" + str(i.result) + ")")

        for k, v in self.items():
            out.append("[ " + k + " ]")
            out.extend(list(str_add("    ", v.__repr__().split("\n"))))

        out = "\n".join(out)
        return out

    def to_file(self, file, full=False):
        with open(file, "w") as outfile:
            if full:
                outfile.write(self.__repr__())
            else:
                outfile.write(str(self.warnings))
                outfile.write("\n")
                outfile.write(str(self.errors))

    @property
    def warnings(self):
        out = CheckResult()
        for i in self.result:
            if i.result == ResultCode.WARNING:
                out.add(i)

        for k, v in self.items():
            if v.contains_warnings():
                out[k].add(v.warnings)

        return out

    @property
    def errors(self):
        out = CheckResult()
        for i in self.result:
            if not i:
                out.add(i)

        for k, v in self.items():
            if not v:
                out[k].add(v.errors)

        return out
