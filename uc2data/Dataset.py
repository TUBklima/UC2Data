from __future__ import annotations
import pyproj
import xarray
import numpy
import csv
import re
import calendar
import netCDF4
import importlib
import pathlib
import urllib.request
from cached_property import cached_property
from .utils import check_type, check_person_field, compare_utms
from .Result import ResultCode, CheckResult

libpath = pathlib.Path(importlib.import_module("uc2data").__file__)
respath = libpath.parent / "resources"
aggregations_file = respath / "aggregations.csv"
variables_file = respath / "uc2_table_A1.csv"
institutions_file = respath / "uc2_table_A3.csv"
data_content_file = respath / "uc2_table_A2.csv"
sites_file = respath / "uc2_table_A4.csv"


class Dataset:
    """
    This object represents data of NetCDF files following the UC2 data standard.

    Attributes
    ----------
    path : str or pathlib.Path
        The path of the file that is to be/was read
    check_result : Dataset.CheckResult
        The results of a call to Dataset.uc2_check
    ds : xarray.Dataset
        The representation of the data

    Methods
    -------
    uc2_check()
        runs the check for conformity with UC2 data standard
    _check_coordinates()
        checks coordinates according to UC2 data standard (only called internally)
    _check_geo_vars(lon_name, lat_name, eutm_name, nutm_name)
        checks coordinates according to UC2 data standard (only called internally)
    check_xy(xy)
        checks horizontal spatial reference variables according to UC2 data standard
    check_var(varname, varname, must_exist, allowed_types=None, allowed_range=None, dims=None, must_be_sorted_along=None, decrease_sort_allowed=True, fill_allowed=True)
        checks NetCDF variable for requested properties
    check_var_attr(varname, attrname, must_exist, allowed_types=None, allowed_values=None, regex=None, must_not_exist=None, allowed_range=None)
        checks NetCDF variable attribute for requested properties
    check_glob_attr(attrname, must_exist, allowed_types=None, allowed_values=None, max_strlen=None, regex=None, allowed_range=None)
        checks NetCDF global attribute for requested properties
    check_dims()
        checks_NetCDF dimensions according to UC2 data standard
    _check_all_vars()
        combines the variable checks of uc2_check method (only called internally)
    check_all_glob_attr()
        combines the global attribute checks of uc2_check method
    _check_cell_methods_agg_varname(varname)
        checks whether the varname matches the implications of cell_methods attribute (only called internally)
    _check_cell_methods_attribute(varname, is_agg_name)
        checks the cell_methods variable attribute (only called internally)
    geo2utm(x,y)
        converts geographic longitude / latitude to UTM coordinates
    """

    allowed_featuretypes = ["timeSeries", "timeSeriesProfile", "trajectory"]

    # try to download newest versions of the table files
    try:
        urllib.request.urlretrieve('http://www.uc2-program.org/uc2_table_A1.csv', variables_file)
    except Exception:
        pass
    try:
        urllib.request.urlretrieve('http://www.uc2-program.org/uc2_table_A2.csv', data_content_file)
    except Exception:
        pass
    try:
        urllib.request.urlretrieve('http://www.uc2-program.org/uc2_table_A3.csv', institutions_file)
    except Exception:
        pass
    try:
        urllib.request.urlretrieve('http://www.uc2-program.org/uc2_table_A4.csv', sites_file)
    except Exception:
        pass

    allowed_aggregations = dict()
    with open(aggregations_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for row in spamreader:
            allowed_aggregations[row[1]] = row[2]

    allowed_data_contents = list()
    with open(data_content_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for row in spamreader:
            allowed_data_contents.append(row[1])

    allowed_variables = dict()
    with open(variables_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for row in spamreader:
            allowed_variables[row[3]] = {
                "long_name": row[0],
                "standard_name": row[1]
            }

    allowed_institutions = []
    allowed_acronyms = []
    with open(institutions_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for row in spamreader:
            allowed_institutions.append(row[0])  # german name
            allowed_acronyms.append(row[1])
            allowed_institutions.append(row[2])  # english name
            allowed_acronyms.append(row[1])

    allowed_locations = []
    allowed_sites = []
    with open(sites_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for row in spamreader:
            allowed_locations.append(row[0])
            allowed_sites.append(row[1])

    allowed_licences = [
        "",
        "[UC]2 MOSAIK Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
        "[UC]2 3DO Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
        "[UC]2 KliMoPrax Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
        "[UC]2 UseUClim Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
        "[UC]2 Restricted Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
        "[UC]2 Research Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
        "[UC]2 Open Licence; see [UC]2 data policy available at www.uc2-program.org/uc2_data_policy.pdf",
    ]

    def __init__(self, path):
        """
        returns a Dataset object

        Parameters
        ----------
        path : str or pathlib.Path
            The path of the file to be read
        """

        self.path = path
        self.check_result = None

        # decode and mask are False for checking file without xarray's interpretation
        self.ds = xarray.open_dataset(self.path, decode_cf=False, mask_and_scale=False)

    @cached_property
    def is_ts(self):
        return self.featuretype == "timeSeries"

    @cached_property
    def is_tsp(self):
        return self.featuretype == "timeSeriesProfile"

    @cached_property
    def is_traj(self):
        return self.featuretype == "trajectory"

    @cached_property
    def is_grid(self):
        return self.featuretype == "None"

    @cached_property
    def is_iop(self):
        if "campaign" in self.ds.attrs:
            return self.ds.campaign[:3] == "IOP"
        else:
            return False

    @cached_property
    def is_lto(self):
        if "campaign" in self.ds.attrs:
            return self.ds.campaign == "LTO"
        else:
            return False

    @cached_property
    def featuretype(self):
        if "featureType" in self.ds.attrs:
            return self.ds.featureType
        else:
            return "None"

    @cached_property
    def filename(self):
        attrs = ["campaign", "location", "site", "acronym", "data_content", "origin_time", "version"]
        vals = list()

        if self.check_result is None:
            self.uc2_check()

        for i in attrs:
            if not self.check_result[i]:
                raise Exception(
                    "Cannot parse filename. Global attribute '" + i + "' did not pass UC2 conformity tests.")

            if i == "origin_time":
                vals.append(self.ds.attrs[i][: 10].replace("-", ""))
            elif i == "version":
                vals.append(str(self.ds.attrs[i]).zfill(3))
            else:
                vals.append(self.ds.attrs[i].replace("-", "_"))
                vals.append(self.ds.attrs[i].replace(".", "_"))
                vals.append(self.ds.attrs[i].replace("/", "_"))
                vals.append(self.ds.attrs[i].replace("\\", "_"))

        filename = "_".join(vals) + ".nc"
        return filename

    def uc2_check(self):
        """
        Performs all checks of conformity to the UC2 data standard.

        The results will be stored in the attribute check_result.

        Returns
        -------
        None

        """

        self.check_result = CheckResult()

        ###
        # Check global attributes
        ###

        self.check_all_glob_attr()
        if not self.check_result["featureType"]:
            return  # doesnt make sense to check file with unknown featureType

        ###
        # Check dims
        ###

        self.check_dims()

        ###
        # Check variables
        ###

        self._check_all_vars()

        # TODO: If all variables have cell_methods with time:point then no time_bounds (and bounds attribute)
        # TODO: If all variables have cell_methods with z:point then no z_bounds (and bounds attribute)
        # TODO: parse cell_methods more nicely: allow "lat: time: lon: mean over years z: sum"
        #       - check if dims have bounds
        # TODO: Make fast check without reading complete arrays?

        ###
        # Check geo vars
        ###

        if self.check_result["crs"]:
            self._check_coordinates()
        else:
            self.check_result["coordinate_transform"].add(ResultCode.ERROR, "Cannot check geographic coordinates " +
                                                          "because of error in 'crs' variable.")

    def _check_coordinates(self):
        """
        updates the attribute check_result with results of coordinates.

        Checks whether lat/lon matches UTM coordinates.

        Returns
        -------
        None

        """

        # Check if origin_lon/origin_lat matches origin_x/origin_y
        if all([self.check_result["origin_lon"], self.check_result["origin_lat"], self.check_result["origin_x"],
                self.check_result["origin_y"]]):
            e_orig_ll, n_orig_ll = self.geo2utm(self.ds.origin_lon, self.ds.origin_lat)

            self.check_result["origin_coords_match"].add(
                compare_utms(e_orig_ll, n_orig_ll, self.ds.origin_x, self.ds.origin_y))
        else:
            self.check_result["origin_coords_match"].add(
                ResultCode.ERROR, "Cannot check if origin_lon/lat matches origin_x/y because of error "
                                  "in one of these global attributes"
            )

        # Check if lon/lat matches E_UTM/N_UTM

        coord_list = [["lon", "lat", "E_UTM", "N_UTM"],  # standard coordinates
                      ["lonu", "latu", "Eu_UTM", "Nu_UTM"],  # coordinates of u-point in grid box
                      ["lonv", "latv", "Ev_UTM", "Nv_UTM"],  # coordinates of v-point in grid box
                      ["lons", "lats", "Es_UTM", "Ns_UTM"]]  # coordinates of surfaces

        for i_coord in coord_list:
            if all(elem in self.ds.variables for elem in i_coord):
                if all(self.check_result[x] for x in i_coord):
                    self.check_result["_".join(i_coord)].add(self._check_geo_vars(*i_coord))

    def _check_geo_vars(self, lon_name, lat_name, eutm_name, nutm_name):
        """
        Checks consistency of horizontal spatial variables

        Parameters
        ----------
        lon_name : str
            variable name for the longitudes
        lat_name : str
            variable name for the latitudes
        eutm_name : str
            variable name for the UTM eastings
        nutm_name : str
            variable name for the UTM northings


        Returns
        -------
        Dataset.CheckResult: The check results concerning these checks.

        """

        x = self.ds[lon_name].values.flatten()
        y = self.ds[lat_name].values.flatten()
        if self.ds[lon_name].dims != self.ds[eutm_name].dims:  # "inflate" array to y,x dims
            # this case is for un-rotated grid with E_UTM(x), N_UTM(y), lon(y,x), lat(y,x)
            e_utm = numpy.tile(self.ds[eutm_name].values, (self.ds[nutm_name].shape[0], 1))
            n_utm = numpy.tile(self.ds[nutm_name].values, (self.ds[eutm_name].shape[0], 1))
            n_utm = numpy.transpose(n_utm)
        else:
            e_utm = self.ds[eutm_name].values
            n_utm = self.ds[nutm_name].values
        e_utm = e_utm.flatten()
        n_utm = n_utm.flatten()

        # Check if fill values are at the same spot and remove them prior to comparison
        xfill = x == -9999
        yfill = y == -9999
        e_ut_mfill = e_utm == -9999
        n_ut_mfill = n_utm == -9999
        if not numpy.array_equal(xfill, yfill) or \
                not numpy.array_equal(xfill, e_ut_mfill) or \
                not numpy.array_equal(xfill, n_ut_mfill):
            return CheckResult(ResultCode.ERROR, "Coordinates have fill values at different indices: " +
                               ", ".join([lon_name, lat_name, eutm_name, nutm_name]) +
                               ". They should be parallel.")

        eutm, nutm = self.geo2utm(x[~xfill], y[~yfill])

        return compare_utms(eutm, nutm, e_utm[~e_ut_mfill], n_utm[~n_ut_mfill])

    def check_xy(self, xy):
        """
        Checks the spatial reference variable for consistency with the UC2 data standard.

        Parameters
        ----------
        xy : str
            the name of the variable to be checked

        Returns
        -------
        Dataset.CheckResult: The result of these checks

        """

        # is xy a box border?
        if xy in ["xu", "yv", "zw", "Eu_UTM", "Ev_UTM", "Nu_UTM", "Nv_UTM",  # TODO: this function is never called with "zw". Do we check it anywhere?
                  "lonu", "lonv", "latu", "latv"]:

            fill_allowed = False  # grid coordinate variables may not have fill values
            if xy in ["xu", "yv", "zw"]:
                dims = xy
                sort_along = xy
            elif xy in ["lonu", "lonv", "latu", "latv"]:
                dims = ("yv", "xu")
                sort_along = None
            else:
                if xy in ["Eu_UTM", "Ev_UTM"]:
                    dims = ["xu", ("yv", "xu")]  # can be 1-dim or 2-dim (if rotated)
                    sort_along = "xu"
                else:
                    dims = ["yv", ("yv", "xu")]  # can be 1-dim or 2-dim (if rotated)
                    sort_along = "yv"

        # surface coordinate?
        elif xy in ["xs", "ys", "lons", "lats", "Es_UTM", "Ns_UTM"]:
            dims = "s"
            sort_along = None
            fill_allowed = False

        # standard corrdinate?
        elif xy in ["x", "y", "lon", "lat", "E_UTM", "N_UTM"]:
            if self.is_grid:
                if "ncol" in self.ds.dims:  # pixel-based surfaces
                    dims = ("nrow", "ncol")
                    sort_along = None
                else:
                    if xy in ["x", "y"]:
                        dims = xy
                        sort_along = xy
                    elif xy in ["lon", "lat"]:
                        dims = ("y", "x")
                        sort_along = None
                    elif xy in ["E_UTM", "N_UTM"]:
                        if xy == "E_UTM":
                            dims = ["x", ("y", "x")]
                            sort_along = "x"
                        else:
                            dims = ["y", ("y", "x")]
                            sort_along = "y"
                    fill_allowed = False
            elif self.is_ts or self.is_tsp:
                dims = "station"
                sort_along = None
                fill_allowed = True
            elif self.is_traj:
                dims = ("traj", "ntime")
                sort_along = None
                fill_allowed = True
            else:
                raise Exception("Unexpected featureType")
        else:
            raise Exception('Unexpected variable: ' + xy)

        out = CheckResult(ResultCode.OK)
        out["variable"].add(self.check_var(xy, True,
                                           allowed_types=[int, float],
                                           dims=dims, must_be_sorted_along=sort_along, decrease_sort_allowed=True,
                                           fill_allowed=fill_allowed))

        if out["variable"]:  # if no error yet

            if xy in ["x", "xs", "y", "ys", "xu", "yv"]:
                long_n = "distance to origin in " + xy[0] + "-direction"
                standard_n = None
                axis = xy[0].upper()
                units = "m"
            elif xy in ["lon", "lons", "lat", "lats", "lonu", "lonv", "latu", "latv"]:
                if xy.startswith("lon"):
                    long_n = "longitude"
                    units = "degrees_east"
                else:
                    long_n = "latitude"
                    units = "degrees_north"
                standard_n = long_n
                axis = None
            elif xy in ["E_UTM", "Es_UTM", "N_UTM", "Ns_UTM", "Eu_UTM", "Ev_UTM", "Nu_UTM", "Nv_UTM"]:
                if xy.startswith("E"):
                    long_n = "easting"
                    standard_n = "projection_x_coordinate"
                else:
                    long_n = "northing"
                    standard_n = "projection_y_coordinate"
                axis = None
                units = "m"

            out["standard_name"].add(self.check_var_attr(xy, "standard_name", xy not in ["x", "xs", "y", "ys"],
                                                         must_not_exist=xy in ["x", "xs", "y", "ys"],
                                                         allowed_values=standard_n))
            out["long_name"].add(self.check_var_attr(xy, "long_name", True, allowed_types=str,
                                                     allowed_values=long_n))
            out["units"].add(self.check_var_attr(xy, "units", True, allowed_types=str,
                                                 allowed_values=units))
            out["axis"].add(self.check_var_attr(xy, "axis", axis is not None, allowed_types=str, allowed_values=axis))

        self.check_result[xy].add(out)

    def check_var(self, varname, must_exist, allowed_types=None, allowed_range: list = None, dims=None,
                  must_be_sorted_along=None, decrease_sort_allowed=True, fill_allowed=True):
        """
        Checks a NetCDF variable for requested properties

        Parameters
        ----------
        varname : str
            name of the variable to check
        must_exist : bool
            whether the variable needs to be present
        allowed_types : Union[str, Iterable], optional
            str or list of types that are allowed. If float then any Python/numpy float allowed. If int then any
            Python/numpy int allowed
        allowed_range : list, optional
            two-element list as [min, max] of allowed range
        dims : Union[str, list, tuple], optional
            The name(s) of the required dimension(s).
            If str: only this dimension allowed
            If list of strings: only these dimensions in the specified order allowed
            If list of tuples: Any of the dimension combinations as given by tuples allowed
            If tuple of strings: Only these dimensions allowed
            Empty tuple: Must have no dimensions (=scalar)
        must_be_sorted_along : str, optional
            Name of a dimension along which the data must be sorted (ignoring fill values)
        decrease_sort_allowed : bool, optional
            Whether variable may be sorted in descending order
        fill_allowed : bool, optional
            whether the variable may contain fill values

        Returns
        -------
        Dataset.CheckResult: the results of the variable check

        """

        exists = varname in self.ds.variables
        result = CheckResult(ResultCode.OK)

        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR, "Required variable '" + varname + "' not found.")
            else:
                return result

        this_var = self.ds[varname]

        if allowed_types is not None:

            if not check_type(this_var, allowed_types):
                result.add(ResultCode.ERROR, "Variable '" + varname + "' has wrong type. " +
                           "Should be one of the following: " + str(allowed_types) + ". " +
                           "Found type: " + str(this_var.dtype))

        if allowed_range is not None:
            if "_FillValue" in this_var.attrs:
                any_value = this_var.values[this_var.values != this_var._FillValue][0]  # the first occurrence of non-FillValue
                this_tmp = numpy.where(this_var.values == this_var._FillValue, any_value, this_var.values)
                this_var_min = this_tmp.min()
                this_var_max = this_tmp.max()
            else:
                this_var_min = this_var.min()
                this_var_max = this_var.max()
            if (this_var_min < allowed_range[0]) or (this_var_max > allowed_range[1]):
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "' is outside allowed range" + str(allowed_range) + ". " +
                           "Found range: [" + str(this_var_min) + "," + str(this_var_max) + "]")

        if dims is not None:
            if type(dims) == str:  # dims can be scalar string
                dims = (dims,)

            if type(dims) == tuple:  # make list of tuples of dims
                dims = [dims]

            for i_dim in range(len(dims)):
                if type(dims[i_dim]) == str:  # if dims is str within list make tuple in place
                    dims[i_dim] = (dims[i_dim],)

            if this_var.dims not in dims:
                result.add(ResultCode.ERROR, "Variable '" + varname + "' has wrong dimensions. Expected: " +
                           str(dims) + ". Found: " + str(this_var.dims))

        if must_be_sorted_along is not None:
            if must_be_sorted_along in this_var.dims:
                me = this_var.values.copy()
                me[me == -9999] = numpy.max(
                    me) + 1  # fill values must be in the end of array for coordinate variables (sort to end)
                sorted_arr = numpy.sort(me, axis=this_var.dims.index(must_be_sorted_along))

                if decrease_sort_allowed:
                    anti_sorted_arr = -numpy.sort(-me, axis=this_var.dims.index(must_be_sorted_along))
                    if not (numpy.array_equal(me, sorted_arr) or numpy.array_equal(me, anti_sorted_arr)):
                        result.add(ResultCode.ERROR,
                                   "Variable '" + varname + "' must be sorted along dimension '" + must_be_sorted_along + "'")
                else:
                    if not numpy.array_equal(me, sorted_arr):
                        result.add(ResultCode.ERROR,
                                   "Variable '" + varname + "' must be sorted along dimension '" + must_be_sorted_along + "'")
            else:
                result.add(ResultCode.ERROR, "Variable should be sorted along " + str(must_be_sorted_along) +
                           " but dim not found in variable.")

        if not fill_allowed:
            if "_FillValue" in this_var.attrs:
                result.add(ResultCode.WARNING, "Variable '" + varname + "' must not contain fill values but has " +
                           "the variable attribute '_FillValue'.")

            if -9999 in this_var:  # -9999 must always be the fill value
                result.add(ResultCode.ERROR, "Variable '" + varname + "' contains -9999. No fill values " +
                           "are allowed for this variables. -9999 is the fixed fill value in UC2 data standard.")

        return result

    def check_var_attr(self, varname, attrname, must_exist, allowed_types=None, allowed_range=None,
                       allowed_values=None, regex=None, must_not_exist=None):

        """
        Checks a NetCDF variable attributes for requested properties

        Parameters
        ----------
        varname : str
            name of the variable to check
        attrname : str
            name of the variable attribute to check
        must_exist : bool
            whether the variable attribute needs to be present
        allowed_types : Union[str, Iterable], optional
            str or list of types that are allowed. If float then any Python/numpy float allowed. If int then any
            Python/numpy int allowed
        allowed_range : list, optional
            two-element list as [min, max] of allowed range
        allowed_values : list
            list of values that are allowed. The actual attribute value must be within the list.
        regex : str
            regular expression that the attribute must match
        must_not_exist : bool
            if True, the variable attribute must not be defined

        Returns
        -------
        Dataset.CheckResult: the results of the variable attribute check

        """

        exists = attrname in self.ds[varname].attrs
        result = CheckResult(ResultCode.OK)
        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR,
                                   "Variable '" + varname + "': Required variable attribute '" + attrname + "' not found.")
            else:
                return result
        else:
            if must_not_exist:
                return CheckResult(ResultCode.ERROR,
                                   "Variable '" + varname + "' has attribute '" + attrname + "' defined. Not allowed.")

        this_value = self.ds[varname].attrs[attrname]

        if allowed_types is not None:
            if not check_type(this_value, allowed_types):
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "': Required variable attribute '" + attrname + "' has wrong type. Should be " +
                           "one of the following: " + str(allowed_types) + ". Found type: " + str(type(this_value)))
                return result

        if allowed_values is not None:
            if type(allowed_values) != list:
                allowed_values = [allowed_values]
            if this_value not in allowed_values:
                if len(allowed_values) == 1:
                    result.add(ResultCode.ERROR,
                               "Variable '" + varname + "': Required variable attribute '" + attrname + "'  has wrong value. Should be " +
                               str(allowed_values[0]) + ". Found value: " + str(this_value))
                else:
                    result.add(ResultCode.ERROR,
                               "Variable '" + varname + "': Required variable attribute '" + attrname + "' has wrong value. " +
                               "Found value: " + str(this_value))

        if allowed_range is not None:
            if this_value < allowed_range[0] or this_value > allowed_range[1]:
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "': Attribute '" + attrname + "' outside range. " +
                           "Expected: " + str(allowed_range) + ". " +
                           "Found: [" + str(numpy.min(this_value)) + "," + str(numpy.max(this_value)) + "]")

        if regex is not None:
            if re.fullmatch(regex, this_value) is None:
                result.add(ResultCode.ERROR,
                           "Global attribute '" + attrname + "' does not match regular expression " + regex + ". " +
                           "Found value: " + str(this_value))
        return result

    def check_glob_attr(self, attrname, must_exist, allowed_types=None, allowed_range=None, allowed_values=None,
                        max_strlen=None, regex=None):

        """
        Checks a NetCDF global attributes for requested properties

        Parameters
        ----------
        attrname : str
            name of the global attribute to check
        must_exist : bool
            whether the global attribute needs to be present
        allowed_types : Union[str, Iterable], optional
            str or list of types that are allowed. If float then any Python/numpy float allowed. If int then any
            Python/numpy int allowed
        allowed_range : list, optional
            two-element list as [min, max] of allowed range
        allowed_values : list
            list of values that are allowed. The actual attribute value must be within the list.
        max_strlen : int
            if the attribute is of type str then the length must not be larger than max_strlen
        regex : str
            regular expression that the attribute must match

        Returns
        -------
        Dataset.CheckResult: the results of the global attribute check

        """

        exists = attrname in self.ds.attrs
        result = CheckResult(ResultCode.OK)

        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR, "Required global attribute '" + attrname + "' not found.")
            else:
                return result

        this_value = self.ds.attrs[attrname]

        if allowed_types is not None:
            if not check_type(this_value, allowed_types):
                result.add(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong type. " +
                           "Should be one of the following: " + str(allowed_types) + ". " +
                           "Found type: " + str(type(this_value)))
                return result

        if not numpy.isscalar(this_value):
            result.add(ResultCode.ERROR, "Global attribute '" + attrname + "' must be scalar but is " +
                       str(type(this_value)))
            return result

        if allowed_values is not None:
            if type(allowed_values) != list:
                allowed_values = [allowed_values]
            if this_value not in allowed_values:
                if len(allowed_values) == 1:
                    result.add(ResultCode.ERROR,
                               "Global attribute '" + attrname + "' has wrong value. " +
                               "Should be " + str(allowed_values[0]) + ". " +
                               "Found value: " + str(this_value))
                else:
                    result.add(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong value. " +
                               "Found value: " + str(this_value))

        if regex is not None:
            if re.fullmatch(regex, this_value) is None:
                result.add(ResultCode.ERROR,
                           "Global attribute '" + attrname + "' does not match regular expression " + regex + ". " +
                           "Found value: " + str(this_value))

        if max_strlen is not None:
            if len(this_value) > max_strlen:
                result.add(ResultCode.ERROR,
                           "Global attribute '" + attrname + "' is too long. Must be max. " +
                           str(max_strlen) + " characters.")

        if allowed_range is not None:
            if (this_value < allowed_range[0]) or (this_value > allowed_range[1]):
                result.add(ResultCode.ERROR,
                           "Global attribute '" + attrname + "' is outside allowed range " + str(allowed_range) + ". " +
                           "Found value: " + str(this_value))

        return result

    def check_dims(self):

        """
        Checks dimensions within the NetCDF file for consistency

        This check is called internally by the uc2_check method.
        The check_result attribute of the Dataset object are updated.

        Returns
        -------
        None

        """

        tmp = netCDF4.Dataset(self.path)
        if any([x.isunlimited() for k, x in tmp.dimensions.items()]):
            self.check_result["unlimited_dim"].add(ResultCode.ERROR, "Unlimited dimensions not supported.")
        tmp.close()

        if "nv" in self.ds.dims:
            if self.ds.dims["nv"] != 2:
                self.check_result["nv_is_2"].add(ResultCode.ERROR, "Dimension 'nv' must have size of 2.")
        if "max_name_len" in self.ds.dims:
            if self.ds.dims["max_name_len"] != 32:
                self.check_result["max_name_len_is_32"].add(ResultCode.ERROR,
                                                            "Dimension 'max_name_len' must have size of 32.")

    def _check_all_vars(self):

        """
        Checks variables within the NetCDF file for consistency

        This check is called internally by the uc2_check method and only makes sense within that
        workflow. This is because check results of previous checks are required within this method.
        The check_result attribute of the Dataset object are updated.

        Returns
        -------
        None

        """

        # vrs

        self.check_result["vrs"]["variable"].add(self.check_var("vrs", True, dims=()))
        if self.check_result["vrs"]["variable"]:
            self.check_result["vrs"]["long_name"].add(self.check_var_attr("vrs", "long_name", True, allowed_types=str,
                                                                          allowed_values="vertical reference system"))
            self.check_result["vrs"]["system_name"].add(
                self.check_var_attr("vrs", "system_name", True, allowed_types=str, allowed_values="DHHN2016"))
            self.check_result["vrs"]["standard_name"].add(
                self.check_var_attr("vrs", "standard_name", False, must_not_exist=True))

        # time

        allowed_range = None
        if self.is_ts or self.is_tsp:
            time_dims = ("station", "ntime")
            time_dim_name = "ntime"
        elif self.is_traj:
            time_dims = ("traj", "ntime")
            time_dim_name = "ntime"
        else:
            if "ncol" in self.ds.dims:  # pixel-based surfaces
                pass  # TODO: do anything?
            else:  # is grid
                time_dims = ("time",)
                time_dim_name = "time"

        if self.is_iop:
            allowed_range = [.001, 86400]
        elif self.is_lto:
            if self.check_result["origin_time"]:
                ndays = calendar.monthrange(int(self.ds.origin_time[0:4]), int(self.ds.origin_time[5:7]))[1]
                allowed_range = [.01, ndays * 24 * 60 * 60]

        self.check_result["time"]["variable"].add(
            self.check_var("time", True, allowed_types=[int, float],
                           allowed_range=allowed_range, dims=time_dims,
                           must_be_sorted_along=time_dim_name,
                           decrease_sort_allowed=False,
                           fill_allowed=not self.is_grid))

        # Check that LTO time series have minimum time step of 30 min.
        if self.is_lto:
            if self.check_result["time"]:
                diff_ok = self.ds["time"].diff(time_dim_name) >= 1800  # is difference ok?
                # add 1 column to diff_ok because diff is one column shorter than time variable
                if self.ds["time"].ndim > 1:
                    add_to = numpy.ones((diff_ok.shape[0], 1), dtype=bool)
                    diff_ok = numpy.concatenate((add_to, diff_ok), axis=1)
                else:
                    add_to = numpy.array([True])
                    diff_ok = numpy.concatenate((add_to, diff_ok))
                is_valid = numpy.not_equal(self.ds["time"].values, -9999)  # -9999 is excluded from diff check
                if numpy.any(
                        numpy.logical_and(numpy.logical_not(diff_ok), is_valid)):  # if diff not okay and not -9999 -> error
                    self.check_result.add(ResultCode.ERROR, "Minimum time step in LTO must be 30 minutes")
            else:
                self.check_result["time"]["variable"].add(ResultCode.ERROR, "Cannot check time steps because of previous error in time variable.")

        if self.check_result["time"]["variable"]:
            self.check_result["time"]["long_name"].add(
                self.check_var_attr("time", "long_name", True, allowed_types=str, allowed_values="time"))
            self.check_result["time"]["standard_name"].add(
                self.check_var_attr("time", "standard_name", True, allowed_types=str,
                                    allowed_values="time"))
            self.check_result["time"]["calendar"].add(self.check_var_attr("time", "calendar", True, allowed_types=str,
                                                                          allowed_values="proleptic_gregorian"))
            self.check_result["time"]["axis"].add(
                self.check_var_attr("time", "axis", True, allowed_types=str, allowed_values="T"))
            self.check_result["time"]["fill_values"].add(
                self.check_var_attr("time", "_FillValue", False, allowed_types=self.ds["time"].dtype,
                                    must_not_exist=self.is_grid))
            self.check_result["time"]["units"].add(self.check_var_attr("time", "units", True, allowed_types=str,
                                                                       regex="seconds since [0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00"))
            if self.check_result["origin_time"] and self.check_result["time"]["units"]:
                if not self.ds["time"].units.endswith(self.ds.origin_time):
                    self.check_result["time"]["origin_time"].add(ResultCode.ERROR,
                                                                 "Global attribute 'origin_time' does not match units of variable 'time'.")
            # bounds attributes are checked below together with other variables.

        # z
        # TODO: z does not have to be there (e.g. T2 model output)

        if self.is_grid:
            z_dims = ("z",)
            must_be_sorted_along = "z"
        elif self.is_ts:
            z_dims = ("station",)
            must_be_sorted_along = None
        elif self.is_tsp:
            z_dims = ("station", "ntime", "nz")
            must_be_sorted_along = "nz"
        elif self.is_traj:
            z_dims = ("traj", "ntime")
            must_be_sorted_along = None
        else:
            raise Exception("unexpected featureType.")

        self.check_result["z"]["variable"].add(
            self.check_var("z", True, allowed_types=[int, float], dims=z_dims,
                           must_be_sorted_along=must_be_sorted_along,
                           fill_allowed=not self.is_grid))
        if self.check_result["z"]["variable"]:
            self.check_result["z"]["long_name"].add(
                self.check_var_attr("z", "long_name", True, allowed_types=str, allowed_values="height above origin"))
            self.check_result["z"]["axis"].add(
                self.check_var_attr("z", "axis", True, allowed_types=str, allowed_values="Z"))
            self.check_result["z"]["positive"].add(
                self.check_var_attr("z", "positive", True, allowed_types=str, allowed_values="up"))
            # Bounds will be checked below with all other variables.

            if self.check_result["z"] and self.check_result["origin_z"]:
                self.check_result["z"]["standard_name"].add(
                    self.check_var_attr("z", "standard_name", self.ds.origin_z == 0, allowed_types=str,
                                        allowed_values="height_above_mean_sea_level",
                                        must_not_exist=self.ds.origin_z != 0))

        # x, y
        self.check_xy("x")
        self.check_xy("y")
        self.check_xy("lon")
        self.check_xy("lat")
        self.check_xy("E_UTM")
        self.check_xy("N_UTM")

        if "s" in self.ds.dims:
            self.check_xy("xs")
            self.check_xy("ys")
            self.check_xy("lons")
            self.check_xy("lats")
            self.check_xy("Es_UTM")
            self.check_xy("Ns_UTM")

        if self.is_grid:
            # if one u is there, all are nedded
            if any(elem in self.ds.variables for elem in ["xu", "Eu_UTM", "Nu_UTM", "lonu", "latu"]):
                self.check_xy("xu")
                self.check_xy("Eu_UTM")
                self.check_xy("Nu_UTM")
                self.check_xy("latu")
                self.check_xy("lonu")
            # if one v is there, all are nedded
            if any(elem in self.ds.variables for elem in ["xv", "Ev_UTM", "Nv_UTM", "lonv", "latv"]):
                self.check_xy("yv")
                self.check_xy("Ev_UTM")
                self.check_xy("Nv_UTM")
                self.check_xy("lonv")
                self.check_xy("latv")

        # crs
        self.check_result["crs"]["variable"].add(self.check_var("crs", True, dims=()))
        if self.check_result["crs"]:
            self.check_result["crs"]["standard_name"].add(
                self.check_var_attr("crs", "standard_name", False, must_not_exist=True))
            self.check_result["crs"]["long_name"].add(self.check_var_attr("crs", "long_name", True,
                                                                          allowed_values="coordinate reference system"))
            self.check_result["crs"]["grid_mapping_name"].add(self.check_var_attr("crs", "grid_mapping_name", True,
                                                                                  allowed_values="transverse_mercator"))
            self.check_result["crs"]["semi_major_axis"].add(self.check_var_attr("crs", "semi_major_axis", True,
                                                                                allowed_types=[int, float],
                                                                                allowed_values=6378137))
            self.check_result["crs"]["inverse_flattening"].add(self.check_var_attr("crs", "inverse_flattening", True,
                                                                                   allowed_types=float,
                                                                                   allowed_range=[
                                                                                       298.257222101 - 0.0001,
                                                                                       298.257222101 + 0.0001]))
            self.check_result["crs"]["longitude_of_prime_meridian"].add(
                self.check_var_attr("crs", "longitude_of_prime_meridian",
                                    True,
                                    allowed_types=[int, float],
                                    allowed_values=0))
            self.check_result["crs"]["longitude_of_central_meridian"].add(
                self.check_var_attr("crs", "longitude_of_central_meridian",
                                    True, allowed_types=[int, float],
                                    allowed_values=[3, 9,
                                                    15]))  # TODO: this could lead to comparison error if 3.000001 is set as float
            self.check_result["crs"]["scale_factor_at_central_meridian"].add(
                self.check_var_attr("crs", "scale_factor_at_central_meridian",
                                    True, allowed_types=[int, float],
                                    allowed_range=[0.9996 - 0.0001, 0.9996 + 0.0001]))
            self.check_result["crs"]["latitude_of_projection_origin"].add(
                self.check_var_attr("crs", "latitude_of_projection_origin",
                                    True, allowed_types=[int, float],
                                    allowed_values=0))
            self.check_result["crs"]["false_easting"].add(
                self.check_var_attr("crs", "false_easting", True, allowed_values=500000, allowed_types=[int, float]))
            self.check_result["crs"]["false_northing"].add(
                self.check_var_attr("crs", "false_northing", True, allowed_values=0, allowed_types=[int, float]))
            self.check_result["crs"]["units"].add(self.check_var_attr("crs", "units", True, allowed_values="m"))
            self.check_result["crs"]["epsg_code"].add(self.check_var_attr("crs", "epsg_code", True,
                                                                          allowed_values=["EPSG:25831", "EPSG:25832",
                                                                                          "EPSG:25833"]))

        #
        # other (auxiliary) coordinate variables
        #

        check_platform = False
        if self.is_ts or self.is_tsp:
            check_platform = True
            name = "station_name"
            long_name = "station name"
            dim = "station"
            cf_id = "timeseries_id"
        elif self.is_traj:
            check_platform = True
            name = "traj_name"
            long_name = "trajectory name"
            dim = "traj"
            cf_id = "trajectory_id"

        if check_platform:
            self.check_result[name]["variable"].add(self.check_var(name, True, allowed_types=numpy.dtype("S1"),
                                                                   dims=(dim, "max_name_len")))
            if self.check_result[name]["variable"]:
                self.check_result[name]["long_name"].add(self.check_var_attr(name, "long_name", True, allowed_types=str,
                                                                             allowed_values=long_name))
                self.check_result[name]["standard_name"].add(
                    self.check_var_attr(name, "standard_name", True, allowed_types=str,
                                        allowed_values="platform_name"))
                self.check_result[name]["cf_role"].add(self.check_var_attr(name, "cf_role", True, allowed_types=str,
                                                                           allowed_values=cf_id))

        if self.is_ts or self.is_tsp:
            self.check_result["station_h"]["variable"].add(self.check_var("station_h", True,
                                                                          allowed_types=[int, float], dims="station",
                                                                          fill_allowed=False))
            if self.check_result["station_h"]["variable"]:
                self.check_result["station_h"]["long_name"].add(
                    self.check_var_attr("station_h", "long_name", True, allowed_types=str,
                                        allowed_values="surface altitude"))
                self.check_result["station_h"]["standard_name"].add(
                    self.check_var_attr("station_h", "standard_name", True, allowed_types=str,
                                        allowed_values="surface_altitude"))
                self.check_result["station_h"]["units"].add(self.check_var_attr("station_h", "units", True,
                                                                                allowed_types=str,
                                                                                allowed_values="m"))
        if self.is_traj:
            self.check_result["height"]["variable"].add(self.check_var("height", True, dims=[(), ("traj", "ntime")],
                                                                       allowed_types=[int, float]))
            if self.check_result["height"]["variable"]:
                self.check_result["height"]["long_name"].add(self.check_var_attr("height", "long_name", True,
                                                                                 allowed_types=str,
                                                                                 allowed_values="height above surface"))
                self.check_result["height"]["standard_name"].add(self.check_var_attr("height", "standard_name", True,
                                                                                     allowed_types=str,
                                                                                     allowed_values="height"))
                self.check_result["height"]["units"].add(self.check_var_attr("height", "units", True, allowed_types=str,
                                                                             allowed_values="m"))

        ###
        # Data variables
        ###

        data_content_var_names = list()
        dont_check = ["station_h", "crs", "vrs", "height"]
        known_coordinates = ["station_name", "traj_name",
                             "z", "zw", "zs",
                             "x", "xu", "xs",
                             "y", "yv", "ys",
                             "lon", "lonu", "lonv", "lons",
                             "lat", "latu", "latv", "lats",
                             "E_UTM", "Eu_UTM", "Ev_UTM", "Es_UTM",
                             "N_UTM", "Nu_UTM", "Nv_UTM", "Ns_UTM",
                             "s",
                             "time",
                             "azimuth", "azimuths", "zenith", "zeniths"]
        if self.is_ts:
            data_dims = ("station", "ntime")
        elif self.is_tsp:
            data_dims = ("station", "ntime", "nz")
        elif self.is_traj:
            data_dims = ("traj", "ntime")
        else:
            if "ncol" in self.ds.dims:  # pixel-based surfaces
                data_dims = ("time", "nrow", "ncol")  # TODO: Wird es erlaubt werden, pixel ohne time abzulegen?
            else:
                data_dims = None

        # get all coordinates that appear in this file
        existing_coordinates = list()
        for ikey in self.ds.variables:
            if (ikey in known_coordinates) or ikey.startswith("bands_"):
                existing_coordinates.append(ikey)
        existing_coordinates.sort()

        for ikey in self.ds.variables:
            if ikey in dont_check:
                continue
            is_normal = ikey in self.allowed_variables
            is_agg_name = ikey in [a + "_" + b for a in self.allowed_variables for b in
                                   self.allowed_aggregations]
            is_bounds = ikey.endswith("_bounds")
            is_bands = ikey.startswith("bands_") and not is_bounds
            is_ancillary = ikey.startswith("ancillary_") and not is_bounds
            is_coordinate = ikey in existing_coordinates

            if not any([is_normal, is_agg_name, is_bounds, is_bands, is_ancillary, is_coordinate]):
                self.check_result[ikey].add(ResultCode.ERROR, "'" + ikey + "' is not a supported variable name.")
                continue

            if is_bands:
                self.check_result[ikey].add(
                    self.check_var(ikey, True, dims=ikey, fill_allowed=False, must_be_sorted_along=ikey))
            elif is_bounds:

                main_key = ikey.replace("_bounds", "")
                if main_key not in self.ds.variables:
                    self.check_result[ikey].add(ResultCode.ERROR,
                                                "Variable '" + ikey + "' seems to be a bounds variable " 
                                                                      "but there is no main variable (expected '" +
                                                main_key + "')")
                else:
                    self.check_result[main_key].add(self.check_var_attr(main_key, "bounds", True,
                                                                        allowed_types=str, allowed_values=ikey))
                    self.check_result[ikey].add(self.check_var(ikey, True, allowed_types=self.ds[main_key].dtype,
                                                               dims=self.ds[main_key].dims + ("nv",)))
                    if len(self.ds[ikey].attrs) != 0:
                        self.check_result[ikey]["attributes"].add(ResultCode.ERROR,
                                                                  "Variable '" + ikey + "' must not have any attributes.")
                # Time must be end of time period
                if ikey == "time_bounds":
                    if self.check_result[ikey]:
                        if not self.ds[main_key].equals(self.ds[ikey][..., 1]):
                            self.check_result[ikey]["variable"].add(ResultCode.ERROR,
                                                                    "second column of 'time_bounds' must equal data of variable 'time'")
                    else:
                        self.check_result[ikey]["variable"].add(ResultCode.ERROR,
                                                                "Could not check values of variable '" + ikey + "'" +
                                                                " because of previous error with this variable.")
                # z must be in middle of z bounds
                if ikey == "z_bounds":
                    if self.check_result[ikey]:
                        z_bound_lower = self.ds[ikey][..., 0]
                        z_bound_upper = self.ds[ikey][..., 1]
                        z_bound_mid = z_bound_lower + (z_bound_upper - z_bound_lower) * 0.5
                        if not numpy.allclose(self.ds[main_key].values, z_bound_mid.values):
                            self.check_result[ikey]["variable"].add(ResultCode.ERROR,
                                                                    "values of z must be in the middle between z_bounds.")
                    else:
                        self.check_result[ikey]["variable"].add(ResultCode.ERROR,
                                                                "Could not check values of variable '" + ikey + "'" +
                                                                " because of previous error with this variable.")

            elif is_ancillary:
                # Check ancillary
                # This is an inner loop over all variables again, to find the one that references this ancillary variable.
                for tmpKey in self.ds.variables:
                    if "ancillary_variables" in self.ds[tmpKey].attrs:
                        if ikey in self.ds[tmpKey].ancillary_variables.split(" "):
                            main_var = self.ds[tmpKey]
                            if main_var.dims != self.ds[ikey].dims:
                                self.check_result.add(ResultCode.ERROR, "Dimensions of ancillary variable '" +
                                                      ikey + "' (" + str(self.ds[ikey].dims) + ") must be the same " +
                                                      "as the referencing variable '" + tmpKey + "' (" +
                                                      str(main_var.dims) + ")")

            elif is_coordinate:
                # TODO: actually we may not pass! E.g., time must go through bounds check below!
                pass  # TODO: Check coordinates? azimuth etc could be 0-360
            else:
                if is_normal:
                    expected_data_content = ikey
                elif is_agg_name:
                    expected_data_content = "_".join(ikey.split("_")[:-1])

                if expected_data_content not in data_content_var_names:
                    data_content_var_names.append(expected_data_content)

                # Check var
                self.check_result[ikey]["variable"].add(self.check_var(ikey, True, dims=data_dims))

                # Check obligatory attributes
                self.check_result[ikey]["long_name"].add(self.check_var_attr(ikey, "long_name", True, allowed_types=str,
                                                                             allowed_values=self.allowed_variables[
                                                                                 expected_data_content][
                                                                                 "long_name"]))
                self.check_result[ikey]["units"].add(
                    self.check_var_attr(ikey, "units", True, allowed_types=str))  # TODO: check conversion
                self.check_result[ikey]["_FillValue"].add(
                    self.check_var_attr(ikey, "_FillValue", True, allowed_types=self.ds[ikey].dtype,
                                        allowed_values=-9999))
                self.check_result[ikey]["coordinates"].add(
                    self.check_var_attr(ikey, "coordinates", True, allowed_types=str))
                if self.check_result[ikey]["coordinates"]:
                    this_coords = self.ds[ikey].coordinates.split(" ")
                    this_coords.sort()

                    coords_in_var_not_in_file = set(this_coords).difference(set(existing_coordinates))
                    coords_in_file_not_in_var = set(existing_coordinates).difference(set(this_coords))

                    if len(coords_in_file_not_in_var) != 0:
                        self.check_result[ikey]["coordinates"].add(ResultCode.WARNING,
                                                                   "variable attribute 'coordinates' does not " +
                                                                   "contain all (auxiliary) coordinates. These are missing: " +
                                                                   str(coords_in_file_not_in_var))
                    if len(coords_in_var_not_in_file) != 0:
                        self.check_result[ikey]["coordinates"].add(ResultCode.ERROR,
                                                                   "variable attribute 'coordinates' contains a reference " +
                                                                   "to a coordinate that is not found in file: " +
                                                                   str(coords_in_var_not_in_file))
                self.check_result[ikey]["grid_mapping"].add(
                    self.check_var_attr(ikey, "grid_mapping", True, allowed_types=str,
                                        allowed_values="crs"))
                # other attributes
                self.check_result[ikey]["standard_name"].add(self.check_var_attr(ikey, "standard_name",
                                                                                 self.allowed_variables[
                                                                                     expected_data_content][
                                                                                     "standard_name"] != "",
                                                                                 allowed_types=str,
                                                                                 allowed_values=
                                                                                 self.allowed_variables[
                                                                                     expected_data_content][
                                                                                     "standard_name"],
                                                                                 must_not_exist=
                                                                                 self.allowed_variables[
                                                                                     expected_data_content][
                                                                                     "standard_name"] == ""))
                self.check_result[ikey]["units_alt"].add(
                    self.check_var_attr(ikey, "units_alt", False, allowed_types=str))  # TODO: check conversion
                self.check_result[ikey]["uncertainty_rel"].add(self.check_var_attr(ikey, "uncertainty_rel", False,
                                                                                   allowed_types=float))
                self.check_result[ikey]["processing_level"].add(self.check_var_attr(ikey, "processing_level", False,
                                                                                    allowed_types=int,
                                                                                    allowed_range=[0, 3]))
                self.check_result[ikey]["processing_info"].add(self.check_var_attr(ikey, "processing_info", False,
                                                                                   allowed_types=str))
                self.check_result[ikey]["instrument_name"].add(self.check_var_attr(ikey, "instrument_name", False,
                                                                                   allowed_types=str))
                self.check_result[ikey]["instrument_nr"].add(self.check_var_attr(ikey, "instrument_nr", False,
                                                                                 allowed_types=str))
                self.check_result[ikey]["instrument_sn"].add(self.check_var_attr(ikey, "instrument_sn", False,
                                                                                 allowed_types=str))

                # check cell_methods if variable has name xyz_method
                if is_agg_name:
                    self.check_result[ikey]["cell_methods"].add(self._check_cell_methods_agg_varname(ikey))

                # check cell_methods if cell_methods in variable attributes
                if "cell_methods" in self.ds[ikey].attrs:
                    self.check_result[ikey]["cell_methods"].add(self._check_cell_methods_attribute(ikey, is_agg_name))

                # Check ancillary_variables attribute
                if "ancillary_variables" in self.ds[ikey].attrs:
                    anc_var = self.ds[ikey].ancillary_variables.split(" ")
                    self.check_result[ikey]["ancillary_variables"].add(self.check_var_attr(ikey, "ancillary_variables",
                                                                                           True, allowed_types=str))
                    for i in anc_var:
                        if i not in self.ds.variables:
                            self.check_result[ikey]["ancillary_variables"].add(ResultCode.ERROR,
                                                                               "Expected ancillary variable '" +
                                                                               i + "' not found in file.")

                # Check bounds attribute
                if "bounds" in self.ds[ikey].attrs:
                    if ikey + "_bounds" not in self.ds.variables:
                        self.check_result[ikey]["bounds"].add(ResultCode.ERROR,
                                                              "Expected variable '" + ikey + "_bounds' not found.")
                    self.check_result[ikey]["bounds"].add(self.check_var_attr(ikey, "bounds", True,
                                                                              allowed_types=str,
                                                                              allowed_values=ikey + "_bounds"))

        if len(data_content_var_names) == 0:
            self.check_result.add(ResultCode.ERROR, "No data variable found.")
        elif len(data_content_var_names) == 1:
            if self.check_result["data_content"] and self.ds.data_content != data_content_var_names[0]:
                self.check_result["data_content"].add(ResultCode.ERROR, "Only one data variable found. '" +
                                                      data_content_var_names[0] +
                                                      "'. Expected global attribute 'data_content'" +
                                                      " to be '" + data_content_var_names[0] + "'.")
        else:
            if self.check_result["data_content"] and self.ds.data_content not in self.allowed_data_contents:
                self.check_result["data_content"].add(ResultCode.ERROR, "Multiple data variables found. "
                                                                        "In this case only one of the allowed"
                                                                        "data_content variable categories must be"
                                                                        "used. You used '" + self.ds.data_content + "'.")

    def check_all_glob_attr(self):
        """
        Checks all global attributes within the NetCDF file for consistency

        This check is called internally by the uc2_check method.
        The check_result attribute of the Dataset object are updated.

        Returns
        -------
        None

        """

        self.check_result["title"].add(self.check_glob_attr("title", True, str))
        self.check_result["data_content"].add(
            self.check_glob_attr("data_content", True, str,
                                 allowed_values=Dataset.allowed_data_contents + list(Dataset.allowed_variables.keys()),
                                 max_strlen=16))
        self.check_result["source"].add(self.check_glob_attr("source", True, str))
        self.check_result["version"].add(self.check_glob_attr("version", True, int, allowed_values=list(
            range(1, 1000))))  # TODO: This is going to be checked in DMS
        self.check_result["Conventions"].add(self.check_glob_attr("Conventions", True, str, allowed_values=["CF-1.7"]))
        self.check_result["dependencies"].add(
            self.check_glob_attr("dependencies", True, str))  # TODO: This is going to be checked by DMS
        self.check_result["history"].add(self.check_glob_attr("history", True, str))
        self.check_result["references"].add(self.check_glob_attr("references", True, str))
        self.check_result["comment"].add(self.check_glob_attr("comment", True, str))
        self.check_result["keywords"].add(self.check_glob_attr("keywords", True, str))
        self.check_result["licence"].add(
            self.check_glob_attr("licence", True, str, allowed_values=Dataset.allowed_licences))
        self.check_result["creation_time"].add(
            self.check_glob_attr("creation_time", True, str,
            regex="[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00"))
        self.check_result["origin_time"].add(
            self.check_glob_attr("origin_time", True, str,
            regex="[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00"))
        self.check_result["origin_lon"].add(self.check_glob_attr("origin_lon", True, float,
                                                                 allowed_range=[-180, 180]))
        self.check_result["origin_lat"].add(self.check_glob_attr("origin_lat", True, float,
                                                                 allowed_range=[-90, 90]))
        self.check_result["origin_x"].add(self.check_glob_attr("origin_x", True, float))
        self.check_result["origin_y"].add(self.check_glob_attr("origin_y", True, float))
        self.check_result["rotation_angle"].add(self.check_glob_attr("rotation_angle", True,
                                                                     float,
                                                                     allowed_range=[0, 360]))

        # non-standard checks

        if self.is_grid:
            self.check_result["featureType"].add(ResultCode.OK)
        else:
            self.check_result["featureType"].add(self.check_glob_attr("featureType", False, str,
                                                                      allowed_values=self.allowed_featuretypes))
            if not self.check_result["featureType"]:
                return

        self.check_result["origin_z"].add(self.check_glob_attr("origin_z", True, float,
                                                               allowed_values=0 if (not self.is_grid) else None))
        self.check_result["location"].add(
            self.check_glob_attr("location", True, str, allowed_values=Dataset.allowed_locations))
        self.check_result["site"].add(self.check_glob_attr("site", True, str, allowed_values=Dataset.allowed_sites,
                                                           max_strlen=12)) # TODO: max_strlen gilt nur für UC2 Projekt?
        if self.check_result["location"] and self.check_result["site"]:
            if Dataset.allowed_locations[Dataset.allowed_sites.index(self.ds.site)] != self.ds.location:
                self.check_result["site"].add(ResultCode.ERROR, "site '" + self.ds.site +
                                              "' does not match location '" + self.ds.location + "'")

        self.check_result["institution"].add(
            self.check_glob_attr("institution", True, str, allowed_values=Dataset.allowed_institutions))
        self.check_result["acronym"].add(
            self.check_glob_attr("acronym", True, str, allowed_values=Dataset.allowed_acronyms,
                                 max_strlen=12)) # TODO: max_strlen gilt nur für UC2 Projekt?
        if self.check_result["institution"] and self.check_result["acronym"]:
            if Dataset.allowed_acronyms[Dataset.allowed_institutions.index(self.ds.institution)] != \
                    self.ds.acronym:
                self.check_result["institution"].add(ResultCode.ERROR, "institution '" + self.ds.institution +
                                                     "' does not match acronym '" + self.ds.acronym + "'")

        self.check_result["author"].add(self.check_glob_attr("author", True, str))
        if self.check_result["author"]:
            if self.ds.author != "":
                self.check_result["author"].add(check_person_field(self.ds.author, "author"))

        self.check_result["contact_person"].add(self.check_glob_attr("contact_person", True, str))
        if self.check_result["contact_person"]:
            self.check_result["contact_person"].add(check_person_field(self.ds.contact_person, "contact_person"))

        self.check_result["campaign"].add(self.check_glob_attr("campaign", True, str, regex="^[A-Za-z0-9\._-]+$",
                                                               max_strlen=12)) # TODO: max_strlen gilt nur für UC2 Projekt?
        if self.check_result["campaign"]:
            if self.is_iop:
                try:
                    if (len(self.ds.campaign) != 5) or (not int(self.ds.campaign[3:]) in range(1, 100)):
                        self.check_result["campaign"].add(ResultCode.ERROR,
                                                          "Global attribute 'campaign': If IOP then string must be IOPxx")
                except ValueError:
                    self.check_result["campaign"].add(ResultCode.ERROR, "If global attribute 'campaign' starts with " +
                                                      "'IOP' then numbers must follow.")
            elif self.ds.campaign.startswith("VALR") or self.ds.campaign.startswith("VALM"):
                try:
                    if (len(self.ds.campaign) != 6) or (not int(self.ds.campaign[4:]) in range(1, 100)):
                        self.check_result["campaign"].add(ResultCode.ERROR,
                                                          "Global attribute 'campaign': If VALM/VALR then string must be VALMxx/VALRxx")
                except ValueError:
                    self.check_result["campaign"].add(ResultCode.ERROR, "If global attribute 'campaign' starts with " +
                                                      "'VALM' or 'VALR' then numbers must follow.")

    def _check_cell_methods_agg_varname(self, varname):

        """
        Checks whether variable names match their cell_methods attribute

        This check is called internally by the uc2_check method and only makes sense within that
        workflow. This is because check results of previous checks are required within this method.
        The check_result attribute of the Dataset object are updated.

        Returns
        -------
        None

        """

        out = CheckResult(ResultCode.OK)
        out[varname]["cell_methods"].add(
            self.check_var_attr(varname, "cell_methods", True, allowed_types=str))
        if out[varname]["cell_methods"]:
            this_agg_short = varname.split("_")[-1]
            this_agg_cf = self.allowed_aggregations[this_agg_short]
            if not re.match(r".*?\btime\b( )?:( )?" + re.escape(this_agg_cf) + r"\b",
                            self.ds[varname].cell_methods):
                out[varname]["cell_methods"].add(ResultCode.ERROR,
                                                 "The variable name indicates a " +
                                                 "temporal aggregation. This must be given by cell_methods: " +
                                                 "'time: " + this_agg_cf + "'.")
            if "time_bounds" not in self.ds.variables:
                out[varname]["cell_methods"].add(ResultCode.ERROR,
                                                 "The variable name indicates a " +
                                                 "temporal aggregation. Therefore the variable " +
                                                 "'time_bounds' is needed.")
        return out

    def _check_cell_methods_attribute(self, varname, is_agg_name):
        """
        Checks the cell_method attribute of variables within the NetCDF file for consistency

        This check is called internally by the uc2_check method and only makes sense within that
        workflow. This is because check results of previous checks are required within this method.
        The check_result attribute of the Dataset object are updated.

        Returns
        -------
        None

        """

        out = CheckResult(ResultCode.OK)

        this_cm = self.ds[varname].cell_methods
        if re.match(r".*?\btime\b( )?:", this_cm):  # contains "time:"?
            method = re.match(r".*?\btime\b ?: ?([a-zA-Z]+)", this_cm)  # get method string after "time:"
            if method:
                method = method.groups()[0].lower()
                if method == "point" or method == "mean" or method == "sum":
                    if is_agg_name:
                        out[varname]["cell_methods"].add(
                            ResultCode.ERROR, "Variable attribute 'cell_methods' specifies temporal "
                                              "aggregation '" + method + "'. Variable name must not contain"
                                                                         " agg_method in this case."
                        )
                else:
                    short_method = None
                    for iag, ival in self.allowed_aggregations.items():  # get short version of method for variable name
                        if ival == method:
                            short_method = iag
                            break
                    if not short_method:
                        out[varname]["cell_methods"].add(
                            ResultCode.ERROR, "cell_methods contain 'time:" + method + "'. Method is unsupported."
                        )
                    else:
                        if not varname.endswith("_" + short_method):
                            out[varname]["cell_methods"].add(
                                ResultCode.ERROR, "cell_method contains 'time:" + method + "'. Variable name should"
                                                                                           " end with '_" + short_method + "'."
                            )

                if method != "point":
                    if "bounds" not in self.ds["time"].attrs:
                        out[varname]["cell_methods"].add(
                            ResultCode.ERROR, "Variable '" + varname + "' contains cell methods 'time:...'. A "
                                                                       "variable 'time_bounds' is needed and variable 'time' must contain attribute "
                                                                       "'bounds' (='time_bounds')."
                        )

            else:
                out[varname]["cell_methods"].add(
                    ResultCode.ERROR, "Variable attribute 'cell_methods' has unsupported format after 'time: '"
                )
        return out

    def geo2utm(self, x, y):
        """
        Transforms geographical coordinates (EPSG:4258) to UTM coordinates in the crs of the NetCDF file

        Parameters
        ----------
        x : float
            longitude. Can be numpy or regular python array, python
            list/tuple or scalar
        y : float
            latitude. Can be numpy or regular python array, python
            list/tuple or scalar

        Returns
        -------
        tuple: UTM coordinates

        """

        utm = pyproj.CRS(self.ds["crs"].epsg_code.lower())
        geo = pyproj.CRS("epsg:4258")

        return pyproj.transform(geo, utm, y, x)
