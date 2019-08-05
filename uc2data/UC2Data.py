from __future__ import annotations
import pyproj
import xarray
import numpy
import enum
import csv
import re
import calendar
from numpy.core.defchararray import add as str_add
from collections import OrderedDict
import netCDF4


aggregations_file = "resources/aggregations.txt"
data_content_file = "resources/data_content.txt"
variables_file = "resources/variables.txt"
institutions_file = "resources/institutions.txt"
sites_file = "resources/sites.txt"

class UC2Data(xarray.Dataset):
    all_floats = [float, numpy.float, numpy.float16, numpy.float32, numpy.float64]
    all_ints = [int, numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.int64]

    allowed_featuretypes = ["timeSeries", "timeSeriesProfile", "trajectory"]

    allowed_aggregations = dict()
    with open(aggregations_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter='|', quotechar='"')
        for row in spamreader:
            allowed_aggregations[row[1]] = row[2]

    allowed_data_contents = list()
    with open(data_content_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            allowed_data_contents.append(row[1])

    allowed_variables = dict()
    with open(variables_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            allowed_variables[row[3]] = {
                "long_name": row[0],
                "standard_name": row[1]
            }
    allowed_data_contents.extend(allowed_variables.keys())

    allowed_institutions = []
    allowed_acronyms = []
    with open(institutions_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            allowed_institutions.append(row[0])
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

        self.path = path
        self.featuretype = None
        self.filename = None
        self.is_grid = None
        self.is_ts = None
        self.is_tsp = None
        self.is_traj = None
        self.is_iop = None
        self.is_lto = None
        self.check_result = None
        super().__init__()
        # decode and mask are False for checking file without xarray's interpretation
        tmp = xarray.open_dataset(self.path, decode_cf=False, mask_and_scale=False)
        self.update(tmp, inplace=True)
        self.attrs = tmp.attrs
        tmp.close()

    def uc2_check(self):

        if "featureType" in self.attrs.keys():
            self.featuretype = self.attrs["featureType"]
        else:
            self.featuretype = "None"

        self.is_ts = self.featuretype == "timeSeries"
        self.is_tsp = self.featuretype == "timeSeriesProfile"
        self.is_traj = self.featuretype == "trajectory"
        self.is_grid = self.featuretype == "None"
        self.is_iop = False
        self.is_lto = False
        if "campaign" in self.attrs:
            self.is_iop = self.attrs["campaign"][:3] == "IOP"
            self.is_lto = self.attrs["campaign"] == "LTO"

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

        tmp = netCDF4.Dataset(self.path)
        if any([x.isunlimited() for k, x in tmp.dimensions.items()]):
            self.check_result["unlimited_dim"].add(ResultCode.ERROR, "Unlimited dimensions not supported.")

        if "nv" in self.dims:
            if self.dims["nv"] != 2:
                self.check_result["nv_is_2"].add(ResultCode.ERROR, "Dimension 'nv' must have size of 2.")
        if "max_name_len" in self.dims:
            if self.dims["max_name_len"] != 32:
                self.check_result["max_name_len_is_32"].add(ResultCode.ERROR,
                                                            "Dimension 'max_name_len' must have size of 32.")

        ###
        # Check variables
        ###

        self.check_all_vars()

        # TODO: If all variables have cell_methods with time:point then no time_bounds (and bounds attribute)
        # TODO: If all variables have cell_methods with z:point then no z_bounds (and bounds attribute)

        ###
        # Check geo vars
        ###

        if self.check_result["crs"]:
            self.check_coordinates()
        else:
            self.check_result["coordinate_transform"].add(ResultCode.ERROR, "Cannot check geographic coordinates " +
                                                          "because of error in 'crs' variable.")

    def check_coordinates(self):

        # Check if origin_lon/origin_lat matches origin_x/origin_y
        if all([self.check_result["origin_lon"], self.check_result["origin_lat"], self.check_result["origin_x"],
                self.check_result["origin_y"]]):
            lon_orig = self.attrs["origin_lon"]
            lat_orig = self.attrs["origin_lat"]
            e_orig_ll, n_orig_ll = self.geo2UTM(lon_orig, lat_orig)

            self.check_result["origin_coords_match"].add(
                compare_UTMs(e_orig_ll, n_orig_ll, self.attrs["origin_x"], self.attrs["origin_y"]))

        # Check if lon/lat matches E_UTM/N_UTM
        if all(elem in self.keys() for elem in ["lon", "lat", "E_UTM", "N_UTM"]):
            if all([self.check_result["lon"], self.check_result["lat"], self.check_result["E_UTM"],
                    self.check_result["N_UTM"]]):
                self.check_result["lon_lat_E_UTM_N_UTM"].add(self.check_geo_vars("lon", "lat", "E_UTM", "N_UTM"))
        if all(elem in self.keys() for elem in ["lonu", "latu", "Eu_UTM", "Nu_UTM"]):
            if all([self.check_result["lonu"], self.check_result["latu"], self.check_result["Eu_UTM"],
                    self.check_result["Nu_UTM"]]):
                self.check_result["lonu_latu_Eu_UTM_Nu_UTM"].add(
                    self.check_geo_vars("lonu", "latu", "Eu_UTM", "Nu_UTM"))
        if all(elem in self.keys() for elem in ["lonv", "latv", "Ev_UTM", "Nv_UTM"]):
            if all([self.check_result["lonv"], self.check_result["latv"], self.check_result["Ev_UTM"],
                    self.check_result["Nv_UTM"]]):
                self.check_result["lonv_latv_Ev_UTM_Nv_UTM"].add(
                    self.check_geo_vars("lonv", "latv", "Ev_UTM", "Nv_UTM"))
        if all(elem in self.keys() for elem in ["lons", "lats", "Es_UTM", "Ns_UTM"]):
            if all([self.check_result["lons"], self.check_result["lats"], self.check_result["Es_UTM"],
                    self.check_result["Ns_UTM"]]):
                self.check_result["lons_lats_Es_UTM_Ns_UTM"].add(
                    self.check_geo_vars("lons", "lats", "Es_UTM", "Ns_UTM"))

    def check_geo_vars(self, lon_name, lat_name, eutm_name, nutm_name):

        x = self[lon_name].values.flatten()
        y = self[lat_name].values.flatten()
        if self[lon_name].dims != self[eutm_name].dims:  # "inflate" array to y,x dims
            # must be un-rotated grid with E_UTM(x), N_UTM(y), lon(y,x), lat(y,x)
            E_UTM = numpy.tile(self[eutm_name].values, (self[nutm_name].shape[0], 1))
            N_UTM = numpy.tile(self[nutm_name].values, (self[eutm_name].shape[0], 1))
            N_UTM = numpy.transpose(N_UTM)
        else:
            E_UTM = self[eutm_name].values
            N_UTM = self[nutm_name].values
        E_UTM = E_UTM.flatten()
        N_UTM = N_UTM.flatten()

        # Check if fill values are at the same spot and remove them prior to comparison
        xfill = x == -9999
        yfill = y == -9999
        E_UTMfill = E_UTM == -9999
        N_UTMfill = N_UTM == -9999
        if not numpy.array_equal(xfill, yfill) or \
                not numpy.array_equal(xfill, E_UTMfill) or \
                not numpy.array_equal(xfill, N_UTMfill):
            return CheckResult(ResultCode.ERROR, "Coordinates have fill values at different indices: " +
                               ", ".join([lon_name, lat_name, eutm_name, nutm_name]) +
                               ". They should be parallel.")

        eutm, nutm = self.geo2UTM(x[~xfill], y[~yfill])

        return compare_UTMs(eutm, nutm, E_UTM[~E_UTMfill], N_UTM[~N_UTMfill])

    def check_xy(self, xy):

        if xy in ["xu", "yv", "zw", "Eu_UTM", "Ev_UTM", "Nu_UTM", "Nv_UTM",
                  "lonu", "lonv", "latu", "latv"]:
            fill_allowed = False
            if xy in ["xu", "yv", "zw"]:
                dims = xy
                sort_along = xy
            elif xy in ["lonu", "lonv", "latu", "latv"]:
                dims = ("yv", "xu")
                sort_along = None
            else:
                dims = "xu" if xy in ["Eu_UTM", "Ev_UTM", "Nu_UTM", "Nv_UTM"] else "yv"
                sort_along = dims
        elif xy in ["xs", "ys", "lons", "lats", "Es_UTM", "Ns_UTM"]:
            dims = "s"
            sort_along = None
            fill_allowed = False
        elif xy in ["x", "y", "lon", "lat", "E_UTM", "N_UTM"]:
            if self.is_grid:
                if "ncol" in self.dims:  # pixel-based surfaces
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
                        dims = "x" if xy == "E_UTM" else "y"
                        sort_along = dims
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
        if out["variable"]:

            if xy in ["x", "xs", "y", "ys", "xu", "yv"]:
                long_n = "distance to origin in " + xy[0] + "-direction"
                standard_n = None
                axis = xy.upper()
                units = "m"
            elif xy in ["lon", "lons", "lat", "lats", "lonu", "lonv", "latu", "latv"]:
                long_n = "longitude" if xy[:3] == "lon" else "latitude"
                standard_n = long_n
                axis = None
                units = "degrees_east" if xy[:3] == "lon" else "degrees_north"
            elif xy in ["E_UTM", "Es_UTM", "N_UTM", "Ns_UTM", "Eu_UTM", "Ev_UTM", "Nu_UTM", "Nv_UTM"]:
                long_n = "easting" if xy[0] == "E" else "northing"
                standard_n = "projection_x_coordinate" if xy[0] == "E" else "projection_y_coordinate"
                axis = None
                units = "m"

            out["standard_name"].add(self.check_var_attr(xy, "standard_name", not xy in ["x", "xs", "y", "ys"],
                                                         must_not_exist=xy in ["x", "xs", "y", "ys"],
                                                         allowed_values=standard_n))
            out["long_name"].add(self.check_var_attr(xy, "long_name", True, allowed_types=str,
                                                     allowed_values=long_n))
            out["units"].add(self.check_var_attr(xy, "units", True, allowed_types=str,
                                                 allowed_values=units))
            out["axis"].add(self.check_var_attr(xy, "axis", axis is not None, allowed_types=str, allowed_values=axis))

        self.check_result[xy].add(out)

    def check_var(self, varname, must_exist, allowed_types=None, allowed_range=None, dims=None,
                  must_be_sorted_along=None, decrease_sort_allowed=True, fill_allowed=True):
        exists = varname in self.variables.keys()
        result = CheckResult(ResultCode.OK)

        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR, "Required variable '" + varname + "' not found.")
            else:
                return result

        if allowed_types is not None:
            if type(allowed_types) == type or type(allowed_types) == numpy.dtype:
                allowed_types = [allowed_types]

            # allow all floats?
            if any(i in self.all_floats for i in allowed_types):
                for this in self.all_floats:
                    if this not in allowed_types:
                        allowed_types.append(this)

            # allow all ints?
            if any(i in self.all_ints for i in allowed_types):
                for this in self.all_ints:
                    if this not in allowed_types:
                        allowed_types.append(this)

            if not self[varname].dtype in allowed_types:
                result.add(ResultCode.ERROR, "Variable '" + varname + "' has wrong type. " +
                           "Should be one of the following: " + str(allowed_types) + ". " +
                           "Found type: " + str(self[varname].dtype))

        if allowed_range is not None:
            if (self[varname].min() < allowed_range[0]) or (self[varname].max() > allowed_range[1]):
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "' is outside allowed range" + str(allowed_range) + ". " +
                           "Found ramge: [" + str(self[varname]) + "," + str(self[varname]) + "]")

        if dims is not None:
            if type(dims) == list:
                dims = tuple(dims)
            elif type(dims) == str:
                dims = tuple([dims])
            if self[varname].dims != dims:
                result.add(ResultCode.ERROR, "Variable '" + varname + "' has wrong dimensions. Expected: " +
                           str(dims) + ". Found: " + str(self[varname].dims))

        if must_be_sorted_along is not None:
            if must_be_sorted_along in self[varname].dims:
                me = self[varname].values.copy()
                me[numpy.where(me == -9999)] = numpy.max(
                    me) + 1  # fill values must be in the end of array for coordinate variables (sort to end)
                sorted_arr = numpy.sort(me, axis=self[varname].dims.index(must_be_sorted_along))

                if decrease_sort_allowed:
                    anti_sorted_arr = -numpy.sort(-me, axis=self[varname].dims.index(must_be_sorted_along))
                    if not (numpy.array_equal(me, sorted_arr) or numpy.array_equal(me, anti_sorted_arr)):
                        result.add(ResultCode.ERROR,
                                   "Variable '" + varname + "' must be sorted along dimension '" + must_be_sorted_along + "'")
                else:
                    if not numpy.array_equal(me, sorted_arr):
                        result.add(ResultCode.ERROR,
                                   "Variable '" + varname + "' must be sorted along dimension '" + must_be_sorted_along + "'")

        if not fill_allowed:
            if "_FillValue" in self[varname].attrs.keys():
                result.add(ResultCode.WARNING, "Variable '" + varname + "' must not contain fill values but has " +
                           "the variable '_FillValue'.")

            if -9999 in self[varname]:  # -9999 must always be the fill value
                result.add(ResultCode.ERROR, "Variable '" + varname + "' contains -9999. No fill values " +
                           "are allowed for this variables. -9999 is the fixed fill value in UC2 data standard.")

        return result

    def check_var_attr(self, varname, attrname, must_exist, allowed_types=None, allowed_values=None, regex=None,
                       must_not_exist=None, allowed_range=None):
        exists = attrname in self[varname].attrs.keys()
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

        this_value = self[varname].attrs[attrname]

        if allowed_types is not None:
            if type(allowed_types) == type or type(allowed_types) == numpy.dtype:
                allowed_types = [allowed_types]
            this_type = type(this_value)

            # allow all floats?
            if any(i in self.all_floats for i in allowed_types):
                for this in self.all_floats:
                    if this not in allowed_types:
                        allowed_types.append(this)

            # allow all ints?
            if any(i in self.all_ints for i in allowed_types):
                for this in self.all_ints:
                    if this not in allowed_types:
                        allowed_types.append(this)

            if not this_type in allowed_types:
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "': Required variable attribute '" + attrname + "' has wrong type. Should be " +
                           "one of the following: " + str(allowed_types) + ". Found type: " + str(this_type))
                return result

        if allowed_values is not None:
            if type(allowed_values) != list:
                allowed_values = [allowed_values]
            if not this_value in allowed_values:
                if len(allowed_values) == 1:
                    result.add(ResultCode.ERROR,
                               "Variable '" + varname + "': Required variable attribute '" + attrname + "'  has wrong value. Should be " +
                               str(allowed_values[0]) + ". Found value: " + str(this_value))
                else:
                    result.add(ResultCode.ERROR,
                               "Variable '" + varname + "': Required variable attribute '" + attrname + "' has wrong value. " +
                               "Found value: " + str(this_value))

        if allowed_range is not None:
            if this_value < allowed_range[0] or \
                    this_value > allowed_range[1]:
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "': Attribute '" + attrname + "' outside range. " +
                           "Expected: " + str(allowed_range) + ". " +
                           "Found: [" + str(numpy.min(this_value)) + "," + str(numpy.max(this_value)) + "]")

        if regex is not None:
            if re.fullmatch(regex, self[varname].attrs[attrname]) is None:
                result.add(ResultCode.ERROR,
                           "Global attribute '" + attrname + "' does not match regular expression " + regex + ". " +
                           "Found value: " + str(this_value))
        return result

    def check_glob_attr(self, attrname, must_exist, allowed_types=None, allowed_values=None,
                        max_strlen=None, regex=None, allowed_range=None):
        exists = attrname in self.attrs.keys()
        result = CheckResult(ResultCode.OK)

        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR, "Required global attribute '" + attrname + "' not found.")
            else:
                return result

        this_value = self.attrs[attrname]

        if allowed_types is not None:
            if type(allowed_types) == type or type(allowed_types) == numpy.dtype:
                allowed_types = [allowed_types]

            this_type = type(this_value)

            # allow all floats?
            if any(i in self.all_floats for i in allowed_types):
                for this in self.all_floats:
                    if this not in allowed_types:
                        allowed_types.append(this)

            # allow all ints?
            if any(i in self.all_ints for i in allowed_types):
                for this in self.all_ints:
                    if this not in allowed_types:
                        allowed_types.append(this)

            if not this_type in allowed_types:
                result.add(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong type. " +
                           "Should be one of the following: " + str(allowed_types) + ". " +
                           "Found type: " + str(this_type))
                return result

        if not numpy.isscalar(this_value):
            result.add(ResultCode.ERROR, "Global attribute '" + attrname + "' must be scalar but is " +
                       str(this_type))
            return result

        if allowed_values is not None:
            if numpy.isscalar(allowed_values):
                allowed_values = [allowed_values]
            if not this_value in allowed_values:
                if len(allowed_values) == 1:
                    result.add(ResultCode.ERROR,
                               "Global attribute '" + attrname + "' has wrong value. " +
                               "Should be " + str(allowed_values[0]) + ". " +
                               "Found value: " + str(this_value))
                else:
                    result.add(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong value")

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

    def check_all_vars(self):

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
            if "ncol" in self.dims:  # pixel-based surfaces
                pass  # TODO: Wird es erlaubt, pixel ohne time anzulegen?
            else:
                time_dims = ("time")
                time_dim_name = "time"

        if self.is_iop:
            allowed_range = [.01, 86400]
        elif self.is_lto:
            if self.check_result["origin_time"]:
                ndays = calendar.monthrange(int(self.attrs["origin_time"][0:4]), int(self.attrs["origin_time"][5:7]))[1]
                allowed_range = [.01, ndays * 24 * 60 * 60]

        self.check_result["time"]["variable"].add(
            self.check_var("time", True, allowed_types=[int, float],
                           allowed_range=allowed_range, dims=time_dims,
                           must_be_sorted_along=time_dim_name,
                           decrease_sort_allowed=False,
                           fill_allowed=not self.is_grid))
        # TODO: check max dt in LTO is <= 30 min
        if self.is_lto:
            pass

        if self.check_result["time"]["variable"]:
            # bounds are checked below together with other variables.
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
                self.check_var_attr("time", "_FillValue", False, allowed_types=self["time"].dtype,
                                    must_not_exist=self.is_grid))
            self.check_result["time"]["units"].add(self.check_var_attr("time", "units", True, allowed_types=str,
                                                                       regex="seconds since [0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00"))
            if self.check_result["origin_time"] and self.check_result["time"]["units"]:
                if self.attrs["origin_time"] != self["time"].units[14:]:
                    self.check_result["time"]["origin_time"].add(ResultCode.ERROR,
                                                                 "Global attribute 'origin_time' does not match units of variable 'time'.")

        # z

        if self.is_grid:
            z_dims = ("z")
            must_be_sorted_along = "z"
        elif self.is_ts:
            z_dims = ("station")
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
            if self.check_result["z"]:
                if self.check_result["origin_z"]:
                    self.check_result["z"]["standard_name"].add(
                        self.check_var_attr("z", "standard_name", self.attrs["origin_z"] == 0, allowed_types=str,
                                            allowed_values="height_above_mean_sea_level",
                                            must_not_exist=self.attrs["origin_z"] != 0))

        if self.is_ts or self.is_tsp:
            self.check_result["station_h"].add(self.check_var("station_h", True,
                                                              allowed_types=[int, float],
                                                              dims=("station")))

        # x, y
        self.check_xy("x")
        self.check_xy("y")
        self.check_xy("lon")
        self.check_xy("lat")
        self.check_xy("E_UTM")
        self.check_xy("N_UTM")

        if "s" in self.dims:
            self.check_xy("xs")
            self.check_xy("ys")
            self.check_xy("lons")
            self.check_xy("lats")
            self.check_xy("Es_UTM")
            self.check_xy("Ns_UTM")

        if self.is_grid:
            # if one u is there, all are nedded
            if any(elem in self for elem in ["xu", "Eu_UTM", "Nu_UTM", "lonu", "latu"]):
                self.check_xy("xu")
                self.check_xy("Eu_UTM")
                self.check_xy("Nu_UTM")
                self.check_xy("latu")
                self.check_xy("lonu")
            # if one v is there, all are nedded
            if any(elem in self for elem in ["xv", "Ev_UTM", "Nv_UTM", "lonv", "latv"]):
                self.check_xy("yv")
                self.check_xy("Ev_UTM")
                self.check_xy("Nv_UTM")
                self.check_xy("lonv")
                self.check_xy("latv")

        # crs
        self.check_result["crs"]["variable"].add(self.check_var("crs", True))
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
                                                                                   allowed_range=[298.2572, 298.2573]))
            self.check_result["crs"]["longitude_of_prime_meridian"].add(
                self.check_var_attr("crs", "longitude_of_prime_meridian",
                                    True,
                                    allowed_types=[int, float],
                                    allowed_values=0))
            self.check_result["crs"]["longitude_of_central_meridian"].add(
                self.check_var_attr("crs", "longitude_of_central_meridian",
                                    True, allowed_types=[int, float],
                                    allowed_values=[3, 9, 15]))
            self.check_result["crs"]["scale_factor_at_central_meridian"].add(
                self.check_var_attr("crs", "scale_factor_at_central_meridian",
                                    True, allowed_types=[int, float],
                                    allowed_range=[0.9995, 0.9997]))
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
            id = "timeseries_id"
        elif self.is_traj:
            check_platform = True
            name = "traj_name"
            long_name = "trajectory name"
            dim = "traj"
            id = "trajectory_id"

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
                                                                           allowed_values=id))

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
            self.check_result["height"]["variable"].add(self.check_var("height", True,
                                                                       allowed_types=[int, float]))
            if self.check_result["height"]["variable"]:
                if self["height"].dims != () and self["height"].dims != ("traj", "ntime"):
                    self.check_result["height"]["variable"].add(ResultCode.ERROR,
                                                                "Variable 'height' must either be scalar " +
                                                                "or have dimensions (traj, ntime).")
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

        dv = dict()
        data_content_var_names = list()
        dont_check = ["station_h", "crs", "vrs", "height"]
        known_coordinates = ["station_name", "traj_name",
                             "height",
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
            if "ncol" in self.dims:  # pixel-based surfaces
                data_dims = ("time", "nrow", "ncol")  # TODO: Wird es erlaubt werden, pixel ohne time abzulegen?
            else:
                data_dims = None

        # get all coordinates that appear in this file
        existing_coordinates = list()
        for ikey in self.variables:
            if (ikey in known_coordinates) or ikey.startswith("bands_"):
                existing_coordinates.append(ikey)
        existing_coordinates.sort()

        for ikey in self.variables:
            if ikey in dont_check:
                continue
            is_normal = ikey in self.allowed_variables.keys()
            is_agg = ikey in [a + "_" + b for a in self.allowed_variables.keys() for b in
                              self.allowed_aggregations.keys()]
            is_bounds = ikey.endswith("_bounds")
            is_bands = ikey.startswith("bands_")
            is_ancillary = ikey.startswith("ancillary_")
            is_coordinate = ikey in existing_coordinates

            if not any([is_normal, is_agg, is_bands, is_bounds, is_bands, is_ancillary, is_coordinate]):
                self.check_result[ikey].add(ResultCode.ERROR, "'" + ikey + "' is not a supported variable name.")
                continue

            if is_bands and not is_bounds:  # if is_bands and is_bounds: that would mean, e.g., "bands_xyz_bounds" which is actually only bounds
                # Check bands (bands are coordinate variables => need dim of same name)
                self.check_result[ikey].add(
                    self.check_var(ikey, True, dims=ikey, fill_allowed=False, must_be_sorted_along=ikey))
            elif is_bounds:

                main_key = ikey[:-7]
                if main_key not in self:
                    self.check_result[ikey].add(ResultCode.ERROR,
                                                "Variable '" + ikey + "' seems to be a bounds variable " \
                                                                      "but there is no main variable (expected '" + \
                                                main_key + "')")
                else:
                    self.check_result[main_key].add(self.check_var_attr(main_key, "bounds", True,
                                                                        allowed_types=str, allowed_values=ikey))
                    self.check_result[ikey].add(self.check_var(ikey, True, allowed_types=self[main_key].dtype,
                                                               dims=self[main_key].dims + ("nv",)))
                    if len(self[ikey].attrs) != 0:
                        self.check_result[ikey]["attributes"].add(ResultCode.ERROR,
                                                                  "Variable '" + ikey + "' must not have any attributes.")
                # Time must be end of time period
                if ikey == "time_bounds":
                    if self.check_result[ikey]:
                        if not self[main_key].equals(self[ikey][..., 1]):
                            self.check_result[ikey]["variable"].add(ResultCode.ERROR,
                                                                    "second column of 'time_bounds' must equal data of variable 'time'")
                    else:
                        self.check_result[ikey]["variable"].add(ResultCode.ERROR,
                                                                "Could not check values of variable '" + ikey + "'" +
                                                                " because of previous error with this variable.")
                # z must be in middle of z bounds
                if ikey == "z_bounds":
                    if self.check_result[ikey]:
                        z_bound_lower = self[ikey][..., 0]
                        z_bound_upper = self[ikey][..., 1]
                        z_bound_mid = z_bound_lower + (z_bound_upper - z_bound_lower) * 0.5
                        if not numpy.allclose(self[main_key].values, z_bound_mid.values, equal_nan=True):
                            self.check_result[ikey]["variable"].add(ResultCode.ERROR,
                                                                    "values of z must be in the middle between z_bounds.")
                    else:
                        self.check_result[ikey]["variable"].add(ResultCode.ERROR,
                                                                "Could not check values of variable '" + ikey + "'" +
                                                                " because of previous error with this variable.")

            elif is_ancillary:
                # Check ancillary
                # This is an inner loop over all variables again, to find the one that references this ancillary variable.
                for tmpKey in self.variables:
                    if "ancillary_variables" in self[tmpKey].attrs.keys():
                        if ikey in self[tmpKey].attrs["ancillary_variables"].split(" "):
                            main_var = self[tmpKey]
                            self.check_result[ikey].add(self.check_var(ikey, True, dims=main_var.dims))

            elif is_coordinate:
                pass
            else:
                if is_normal:
                    expected_data_content = ikey
                elif is_agg:
                    expected_data_content = "_".join(ikey.split("_")[:-1])
                else:
                    raise Exception("Unexpected var type: " + ikey)

                if expected_data_content not in data_content_var_names:
                    data_content_var_names.append(expected_data_content)

                # Check var
                self.check_result[ikey]["variable"].add(self.check_var(ikey, True, dims=data_dims))

                # Check obligatory attributes
                self.check_result[ikey]["long_name"].add(self.check_var_attr(ikey, "long_name", True, allowed_types=str,
                                                                             allowed_values=
                                                                             self.allowed_variables[
                                                                                 expected_data_content][
                                                                                 "long_name"]))
                self.check_result[ikey]["units"].add(
                    self.check_var_attr(ikey, "units", True, allowed_types=str))  # TODO: check conversion
                self.check_result[ikey]["_FillValue"].add(
                    self.check_var_attr(ikey, "_FillValue", True, allowed_types=self[ikey].dtype,
                                        allowed_values=-9999))
                self.check_result[ikey]["coordinates"].add(
                    self.check_var_attr(ikey, "coordinates", True, allowed_types=str))
                if self.check_result[ikey]["coordinates"]:
                    this_coords = self[ikey].attrs["coordinates"].split(" ")
                    this_coords.sort()
                    if this_coords != existing_coordinates:
                        self.check_result[ikey]["coordinates"].add(ResultCode.WARNING,
                                                                   "variable attribute 'coordinates' does not " +
                                                                   "contain all (auxiliary) coordinates. These are missing: " +
                                                                   str(set(existing_coordinates).difference(
                                                                       set(this_coords))))

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

                # Check cell_methods
                if is_agg:
                    self.check_result[ikey]["cell_methods"].add(
                        self.check_var_attr(ikey, "cell_methods", True, allowed_types=str))
                    if self.check_result[ikey]["cell_methods"]:
                        this_agg_short = ikey.split("_")[-1]
                        this_agg_cf = self.allowed_aggregations[this_agg_short]
                        if not re.match(r".*?\btime\b( )?:( )?" + re.escape(this_agg_cf) + r"\b",
                                        self[ikey].attrs["cell_methods"]):
                            self.check_result[ikey]["cell_methods"].add(ResultCode.ERROR,
                                                                        "The variable name indicates a " +
                                                                        "temporal aggregation. This must be given by cell_methods: " +
                                                                        "'time: " + this_agg_cf + "'.")

                # Check ancillary_variables attribute
                if "ancillary_variables" in self[ikey].attrs.keys():
                    anc_var = self[ikey].attrs["ancillary_variables"].split(" ")
                    for i in anc_var:
                        if i not in self.keys():
                            self.check_result[ikey]["ancillary_variables"].add(ResultCode.ERROR,
                                                                               "Expected ancillary variable '" +
                                                                               i + "' not found in file.")
                    self.check_result[ikey]["ancillary_variables"].add(self.check_var_attr(ikey, "ancillary_variables",
                                                                                           True, allowed_types=str))

                # Check bounds attribute
                if "bounds" in self[ikey].attrs.keys():
                    if ikey + "_bounds" not in self.keys():
                        self.check_result[ikey]["bounds"].add(ResultCode.ERROR,
                                                              "Expected variable '" + ikey + "_bounds' not found.")
                    self.check_result[ikey]["bounds"].add(self.check_var_attr(ikey, "bounds", True,
                                                                              allowed_types=str,
                                                                              allowed_values=ikey + "_bounds"))

        if len(data_content_var_names) == 0:
            self.check_result.add(ResultCode.ERROR, "No data variable found.")
        elif len(data_content_var_names) == 1:
            if self.check_result["data_content"]:
                if self.attrs["data_content"] != data_content_var_names[0]:
                    self.check_result["data_content"].add(ResultCode.ERROR, "Only one data variable found. '" +
                                                          data_content_var_names[
                                                              0] + "'. Expected global attribute 'data_content'" +
                                                          " to be '" + data_content_var_names[0] + "'.")

    def check_all_glob_attr(self):
        self.check_result["title"].add(self.check_glob_attr("title", True, str))
        self.check_result["data_content"].add(self.check_glob_attr("data_content", True, str,
                                                                   allowed_values=UC2Data.allowed_data_contents,
                                                                   max_strlen=16))
        self.check_result["source"].add(self.check_glob_attr("source", True, str))
        self.check_result["version"].add(self.check_glob_attr("version", True,
                                                              int,
                                                              allowed_values=list(
                                                                  range(1,
                                                                        1000))))  # TODO: This is going to be checked in DMS
        self.check_result["Conventions"].add(self.check_glob_attr("Conventions", True, str, allowed_values=["CF-1.7"]))
        self.check_result["dependencies"].add(self.check_glob_attr("dependencies", True,
                                                                   str))  # TODO: This is going to be checked by DMS
        self.check_result["history"].add(self.check_glob_attr("history", True, str))
        self.check_result["references"].add(self.check_glob_attr("references", True, str))
        self.check_result["comment"].add(self.check_glob_attr("comment", True, str))
        self.check_result["keywords"].add(self.check_glob_attr("keywords", True, str))
        self.check_result["licence"].add(
            self.check_glob_attr("licence", True, str, allowed_values=UC2Data.allowed_licences))
        self.check_result["creation_time"].add(self.check_glob_attr("creation_time", True, str,
                                                                    regex="[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00"))
        self.check_result["origin_time"].add(self.check_glob_attr("origin_time", True, str,
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

        if not self.is_grid:
            self.check_result["featureType"].add(self.check_glob_attr("featureType", False, str,
                                                                      allowed_values=self.allowed_featuretypes))
            if not self.check_result["featureType"]:
                return
        else:
            self.check_result["featureType"].add(ResultCode.OK)

        self.check_result["origin_z"].add(self.check_glob_attr("origin_z", True, float,
                                                               allowed_values=0 if (not self.is_grid) else None))
        self.check_result["location"].add(
            self.check_glob_attr("location", True, str, allowed_values=UC2Data.allowed_locations))
        self.check_result["site"].add(self.check_glob_attr("site", True, str, allowed_values=UC2Data.allowed_sites,
                                                           max_strlen=12))
        if self.check_result["location"] and self.check_result["site"]:
            if UC2Data.allowed_locations[UC2Data.allowed_sites.index(self.attrs["site"])] != self.attrs["location"]:
                self.check_result["site"].add(ResultCode.ERROR, "site '" + self.attrs[
                    "site"] + "' does not match location '" + self.attrs["location"] + "'")
                self.check_result["location"].add(ResultCode.ERROR, "site '" + self.attrs[
                    "site"] + "' does not match location '" + self.attrs["location"] + "'")

        self.check_result["institution"].add(
            self.check_glob_attr("institution", True, str, allowed_values=UC2Data.allowed_institutions))
        self.check_result["acronym"].add(
            self.check_glob_attr("acronym", True, str, allowed_values=UC2Data.allowed_acronyms,
                                 max_strlen=12))
        if self.check_result["institution"] and self.check_result["acronym"]:
            if UC2Data.allowed_institutions.index(self.attrs["institution"]) != UC2Data.allowed_acronyms.index(
                    self.attrs["acronym"]):
                self.check_result["institution"].add(ResultCode.ERROR, "institution '" + self.attrs[
                    "institution"] + "' does not match acronym '" + self.attrs["acronym"] + "'")
                self.check_result["acronym"].add(ResultCode.ERROR, "institution '" + self.attrs[
                    "institution"] + "' does not match acronym '" + self.attrs["acronym"] + "'")

        self.check_result["author"].add(self.check_glob_attr("author", True, str))
        if self.check_result["author"]:
            if self.attrs["author"] != "":
                self.check_result["author"].add(check_person_field(self.attrs["author"], "author"))

        self.check_result["contact_person"].add(self.check_glob_attr("contact_person", True, str))
        if self.check_result["contact_person"]:
            self.check_result["contact_person"].add(check_person_field(self.attrs["contact_person"], "contact_person"))

        self.check_result["campaign"].add(self.check_glob_attr("campaign", True, str, regex="^[A-Za-z0-9\._-]+$",
                                                               max_strlen=12))
        if self.check_result["campaign"]:
            if self.is_iop:
                if (len(self.attrs["campaign"]) != 5) or (not int(self.attrs["campaign"][3:]) in range(1, 100)):
                    self.check_result["campaign"].add(ResultCode.ERROR,
                                                      "Global attribute 'campaign': If IOP then string must be IOPxx")
            elif self.attrs["campaign"][0:4] in ["VALR", "VALM"]:
                if (len(self.attrs["campaign"]) != 6) or (not int(self.attrs["campaign"][4:]) in range(1, 100)):
                    self.check_result["campaign"].add(ResultCode.ERROR,
                                                      "Global attribute 'campaign': If VALM/VALR then string must be VALMxx/VALRxx")

    def get_filename(self):
        attrs = ["campaign", "location", "site", "acronym", "data_content", "origin_time", "version"]
        vals = list()

        if self.check_result is None:
            self.uc2_check()

        for i in attrs:
            if not self.check_result[i]:
                raise Exception("Cannot parse filename. Global attribute '" + i + "' did not pass UC2 conformity tests.")

            if i == "origin_time":
                vals.append(self.attrs[i][: 10].replace("-", ""))
            elif i == "version":
                vals.append(str(self.attrs[i]).zfill(3))
            else:
                vals.append(self.attrs[i].replace("-", "_"))

        filename = "_".join(vals) + ".nc"
        self.filename = filename
        return filename

    def geo2UTM(self, x, y):

        utm = pyproj.CRS(self["crs"].attrs["epsg_code"].lower(), preserve_units=False)
        geo = pyproj.CRS("epsg:4258")

        return pyproj.transform(geo, utm, y, x)


def compare_UTMs(e1, n1, e2, n2):
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
    s = string.split(';')
    for i in s:
        i_s = i.split(',')
        if not len(i_s) in [2, 3]:
            return CheckResult(ResultCode.ERROR,
                               "Global attribute '" + attrname + "': Perons must be given as last_name, first_name[, email]")
        if len(i_s) == 3:
            if re.fullmatch(r"[^@]+@[^@]+\.[^@]+", i_s[2]) is None:
                return CheckResult(ResultCode.ERROR, "Global attribute '" + attrname + "': " + i_s[
                    2] + " is not a valid email address.")
    return CheckResult(ResultCode.OK)


class ResultCode(enum.Enum):
    OK = 1
    WARNING = 2
    ERROR = 3


class ResultItem:

    def __init__(self, result: ResultCode = ResultCode.OK, message: str = ""):
        self.result = result
        if result == ResultCode.OK:
            if message != "":
                raise Exception("cannot handle user message for ResultCode.OK")
            self.message = "Test passed."
        else:
            self.message = message

    def __bool__(self):
        return self.result != ResultCode.ERROR


class CheckResult(OrderedDict):

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
                elif ResultCode.OK in [i.result for i in self.result]:
                    return  # We need only one OK per result
                else:
                    return  # There are ERRORs in result. => not OK
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

    def warnings(self):
        out = CheckResult()
        for i in self.result:
            if i.result == ResultCode.WARNING:
                out.add(i)

        for k, v in self.items():
            if v.contains_warnings():
                out[k].add(v.warnings())

        return out

    def errors(self):
        out = CheckResult()
        for i in self.result:
            if not i:
                out.add(i)

        for k, v in self.items():
            if not v:
                out[k].add(v.errors())

        return out