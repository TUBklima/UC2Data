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
            A message for the user. Only applies if result is not OK. Default: "Test passed."
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

    The object is iterable which yields the unnested results.
    The object also contains tags like a dictionary.
    Each tag contains a deeper (nested) CheckResult object.

    Attributes
    ----------
    result : ResultItem
        Whether the test was passed or generated a warning/error. Includes a message for the user.

    """

    def __init__(self, *args, **kwargs):

        """
        Creates a CheckResult object.

        Parameters
        ----------
        *args
            passed on to add method for initialization with a specific result
        *kwargs
            passed on to add method for initialization with a specific result
        """

        self.result = list()
        if args or kwargs:
            self.add(*args, **kwargs)

    def __getitem__(self, item):

        """
        Get the result(s) of a tag

        Parameters
        ----------
        item
            name of the tag to get the results for
        """

        if item not in super().keys():
            self[item] = CheckResult()  # if there is no tag yet: Create it with an empty CheckResult
        return super().__getitem__(item)

    def __bool__(self):

        """
        Returns True if all results are OK or WARNING.

        If no results are present, False is returned.
        """

        if len(self.result) == 0:
            # an empty thing is not True (otherwise you couldnt check for bool(result["var"]) if "var" is not in result
            ok = len(self.keys()) != 0
        else:
            ok = all(i_ok for i_ok in self.result)

        if not ok:
            return ok

        for val in self.values():  # loop through nested tags
            if not val:
                return False

        return True

    def contains_warnings(self):

        """
        Returns True if there are WARNINGs within the object.
        """

        if len(self.result) == 0:
            has_warn = False
        else:
            has_warn = any([ir.result == ResultCode.WARNING for ir in self.result])  # This checks the unnested results

        if has_warn:
            return has_warn

        for val in self.values():
            if val.contains_warnings():
                return True

        return False

    def add(self, result, message=""):

        """
        Adds a ResultItem to the CheckResult object

        If an error is added to a tag which was OK before, the OK result is removed.
        If there was already an error, the new error is added to the tag.
        If an OK result is added to a tag which has an error or a warning, the OK result is actually not added.

        Parameters
        ----------
        result : Union[ResultCode, ResultItem, CheckResult]
            If result is of type ResultItem then message parameter is ignored
        message : str
            passed on to add method for initialization with a specific result

        Examples
        --------
        >>> d = uc2data.CheckResult()
        >>> d["thisTest"].add(uc2data.ResultCode.OK)  # add a nested OK result
        >>> print(d)
        [ thisTest ]
            Test passed. (ResultCode.OK)
        >>> d["thisTest"].add(uc2data.ResultCode.ERROR, "something went wrong")  # add an ERROR in the same tag (OK is removed)
        >>> print(d)
        [ thisTest ]
            something went wrong (ResultCode.ERROR)
        >>> d.add(uc2data.ResultItem(uc2data.ResultCode.OK))  # add an unnsted OK ResultItem
        >>> print(d)
        Test passed. (ResultCode.OK)
        [ thisTest ]
            something went wrong (ResultCode.ERROR)
        """

        if isinstance(result, ResultCode):
            other = ResultItem(result, message)
        else:
            other = result

        if isinstance(other, ResultItem):
            if other.result == ResultCode.OK:
                if len(self.result) == 0:
                    self.result.append(other)  # no result yet? => add this one
                else:
                    return  # There are ERRORs or WARNINGs in result. => not OK, If OK => stays OK
            else:
                for i in self.result:
                    if i.result == ResultCode.OK:
                        self.result.remove(i)  # remove OK from result because ERROR is added.
                self.result.append(other)

        elif isinstance(other, CheckResult):
            for i in other.result:
                self.add(i)  # add each unnested ResultItem of the other object
            for key, value in other.items():
                self[key].add(value)  # add each nested results of the other object
        else:
            raise Exception("unexpected type of other")

    def __repr__(self):

        """
        Prints the CheckResult in a nice way
        """

        out = list()
        for i in self.result:
            out.append(i.message + " (" + str(i.result) + ")")

        for k, v in self.items():
            out.append("[ " + k + " ]")
            out.extend(list(str_add("    ", v.__repr__().split("\n"))))

        out = "\n".join(out)
        return out

    def to_file(self, file, full=False):
        """
        Writes a file with results

        Parameters
        ----------
        file : str
            the filename
        full : bool
            If set, also OK results will be printed. Default: only ERRORs and WARNINGs
        """
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
