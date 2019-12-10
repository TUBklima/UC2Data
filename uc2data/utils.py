import re
import numpy
from .Result import ResultCode, CheckResult


def check_type(var, allowed_types):

    """
    Checks whether a variable is of an allowed type

    Parameters
    ----------
    var : any
        the variable to check
    allowed_types : list
        list of types that are allowed

    Returns
    -------
    True if the variable is of an allowed type. False otherwise

    """

    all_floats = [float, numpy.float, numpy.float16, numpy.float32, numpy.float64]
    all_ints = [int, numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.int64]

    try:
        this_type = var.dtype
    except AttributeError:
        this_type = type(var)

    if type(allowed_types) != list:
        allowed_types = [allowed_types]

    # allow all floats?
    if any(i in all_floats for i in allowed_types):
        for this in all_floats:
            if this not in allowed_types:
                allowed_types.append(this)

    # allow all ints?
    if any(i in all_ints for i in allowed_types):
        for this in all_ints:
            if this not in allowed_types:
                allowed_types.append(this)

    return this_type in allowed_types


def compare_utms(e1, n1, e2, n2):
    """
    Checks whether pairs of UTM coordinates refer to (roughly) the same location

    A warning is given if coordinate pairs differ by a small distance

    Parameters
    ----------
    e1 : float
        UTM easting(s) of the first point(s). Can be scalar of numpy.array
    n1 : float
        UTM northing(s) of the first point(s). Can be scalar of numpy.array
    e2 : float
        UTM easting(s) of the second point(s). Can be scalar of numpy.array
    n2 : float
        UTM northing(s) of the second point(s). Can be scalar of numpy.array

    Returns
    -------
    Dataset.CheckResult: The result of this check

    """

    if not isinstance(e1, numpy.ndarray):
        e1 = [e1]
    if not isinstance(n1, numpy.ndarray):
        n1 = [n1]
    if not isinstance(e2, numpy.ndarray):
        e2 = [e2]
    if not isinstance(n2, numpy.ndarray):
        n2 = [n2]

    max_diff = max(max(abs(numpy.subtract(e1, e2))),
                   max(abs(numpy.subtract(n1, n2))))

    out = CheckResult()

    if max_diff < 0.1:
        out.add(ResultCode.OK)
    elif max_diff < 1:
        out.add(ResultCode.WARNING, "UTM coordinates in file " +
                "differ from UTM coordinates calculated from lat/lon " +
                "by up to " + str(max_diff) + " m.")
    else:
        out.add(ResultCode.ERROR, "UTM coordinates in file " +
                "do not match calculated UTM coordinates from lat/lon in file")

    return out


def check_person_field(string, attrname):

    """
    Checks whether a string matches the definition of personal details

    personal detail strings must be formated like this:
    "LastName, FirstName[, email]"
    Multiple persons are separated with ";"

    Parameters
    ----------
    string : str
        The string to check
    attrname : str
        Name of the attribute which contains the string

    Returns
    -------
    Dataset.CheckResult: The result of this check

    """

    s = string.split(';')
    for i in s:
        i_s = i.split(',')
        if not len(i_s) in [2, 3]:
            return CheckResult(ResultCode.ERROR,
                               "Global attribute '" + attrname + "': Persons must be given as last_name, first_name[, email]")
        if len(i_s) == 3:
            if re.fullmatch(r"[^@]+@[^@]+\.[^@]+", i_s[2]) is None:
                return CheckResult(ResultCode.ERROR, "Global attribute '" + attrname + "': " + i_s[
                    2] + " is not a valid email address.")
    return CheckResult(ResultCode.OK)
