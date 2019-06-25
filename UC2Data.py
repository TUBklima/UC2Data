import xarray
import numpy
import sys
import enum
import csv
import re
import calendar
from numpy.core.defchararray import add as str_add


is_win = sys.platform in ['win32', 'win64']
if not is_win:
    from cfchecker import cfchecks

data_content_file = "data_content.txt"
variables_file = "variables.txt"
institutions_file = "institutions.txt"
sites_file = "sites.txt"


class UC2Data(xarray.Dataset):

    allowed_data_contents = list()
    with open(data_content_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            allowed_data_contents.append(row[1])

    allowed_variables = list()
    with open(variables_file, encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='"')
        for row in spamreader:
            allowed_variables.append(row[3])
    allowed_data_contents.extend(allowed_variables)

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
        self.check_result = None
        super().__init__()
        tmp = xarray.open_dataset(self.path, decode_cf=False)
        self.update(tmp, inplace=True)
        self.attrs = tmp.attrs

        if "featureType" in self.attrs.keys():
            self.featuretype = self.attrs["featureType"]
        else:
            self.featuretype = "None"

        self.check_result = self.uc2_check()

        self.filename = self.get_filename()

    def uc2_check(self):

        is_ts = self.featuretype == "timeSeries"
        is_tsp = self.featuretype == "timeSeriesProfile"
        is_traj = self.featuretype == "trajectory"
        is_grid = self.featuretype == "None"

        ###
        # Check global attributes
        ###

        result = dict()

        result["title"] = self.check_glob_attr("title", True, str)
        result["data_content"] = self.check_glob_attr("data_content", True, str,
                                                      allowed_values=UC2Data.allowed_data_contents,
                                                      max_strlen=16)  # TODO: Redo this test when variable is checked
        result["source"] = self.check_glob_attr("source", True, str)
        result["version"] = self.check_glob_attr("version", True,
                                                 [int, numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.int64],
                                                 allowed_values=list(
                                                     range(1, 1000)))  # TODO: This is going to be checked in DMS
        result["Conventions"] = self.check_glob_attr("Conventions", True, str, allowed_values=["CF-1.7"])
        result["dependencies"] = self.check_glob_attr("dependencies", True,
                                                      str)  # TODO: This is going to be checked by DMS
        result["history"] = self.check_glob_attr("history", True, str)
        result["references"] = self.check_glob_attr("references", True, str)
        result["comment"] = self.check_glob_attr("comment", True, str)
        result["keywords"] = self.check_glob_attr("keywords", True, str)
        result["licence"] = self.check_glob_attr("licence", True, str, allowed_values=UC2Data.allowed_licences)
        result["creation_time"] = self.check_glob_attr("creation_time", True, str,
                                                       regex="[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00")
        result["origin_time"] = self.check_glob_attr("origin_time", True, str,
                                                     regex="[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00")  # TODO: Check later with time units.
        result["origin_lon"] = self.check_glob_attr("origin_lon", True, [numpy.float, numpy.float32, numpy.float64],
                                                    allowed_range=[-180, 180])
        result["origin_lat"] = self.check_glob_attr("origin_lat", True, [numpy.float, numpy.float32, numpy.float64],
                                                    allowed_range=[-90, 90])
        result["origin_x"] = self.check_glob_attr("origin_x", True, [numpy.float, numpy.float32, numpy.float64])
        result["origin_y"] = self.check_glob_attr("origin_y", True, [numpy.float, numpy.float32, numpy.float64])
        result["rotation_angle"] = self.check_glob_attr("rotation_angle", True,
                                                        [numpy.float, numpy.float32, numpy.float64],
                                                        allowed_range=[0, 360])

        # non-standard checks

        if not is_grid:
            result["featureType"] = self.check_glob_attr("featureType", False, str,
                                                         allowed_values=["timeSeries", "timeSeriesProfile",
                                                                         "trajectory"])
            if not result["featureType"]:
                return result

        result["origin_z"] = self.check_glob_attr("origin_z", True, [numpy.float, numpy.float32, numpy.float64],
                                                  allowed_values=0 if (not is_grid) else None)
        result["location"] = self.check_glob_attr("location", True, str, allowed_values=UC2Data.allowed_locations)
        result["site"] = self.check_glob_attr("site", True, str, allowed_values=UC2Data.allowed_sites,
                                              max_strlen=12)
        if result["location"] and result["site"]:
            if UC2Data.allowed_locations[UC2Data.allowed_sites.index(self.attrs["site"])] != self.attrs["location"]:
                result["site"].append(CheckResult(ResultCode.ERROR, "site '" + self.attrs[
                    "site"] + "' does not match location '" + self.attrs["location"] + "'"))
                result["location"].append(CheckResult(ResultCode.ERROR, "site '" + self.attrs[
                    "site"] + "' does not match location '" + self.attrs["location"] + "'"))

        result["institution"] = self.check_glob_attr("institution", True, str, allowed_values=UC2Data.allowed_institutions)
        result["acronym"] = self.check_glob_attr("acronym", True, str, allowed_values=UC2Data.allowed_acronyms,
                                                 max_strlen=12)
        if result["institution"] and result["acronym"]:
            if UC2Data.allowed_institutions.index(self.attrs["institution"]) != UC2Data.allowed_acronyms.index(self.attrs["acronym"]):
                result["institution"].append(CheckResult(ResultCode.ERROR, "institution '" + self.attrs[
                    "institution"] + "' does not match acronym '" + self.attrs["acronym"] + "'"))
                result["acronym"].append(CheckResult(ResultCode.ERROR, "institution '" + self.attrs[
                    "institution"] + "' does not match acronym '" + self.attrs["acronym"] + "'"))

        result["author"] = self.check_glob_attr("author", True, str)
        if result["author"]:
            if self.attrs["author"] != "":
                result["author"].append(check_person_field(self.attrs["author"], "author"))

        result["contact_person"] = self.check_glob_attr("contact_person", True, str)
        if result["contact_person"]:
            result["contact_person"].append(check_person_field(self.attrs["contact_person"], "contact_person"))

        is_iop = False
        is_lto = False
        result["campaign"] = self.check_glob_attr("campaign", True, str, regex="^[A-Za-z0-9\._-]+$",
                                                  max_strlen=12)
        if result["campaign"]:
            if self.attrs["campaign"][0:3] == "IOP":
                is_iop = True
                if (len(self.attrs["campaign"]) != 5) or (not int(self.attrs["campaign"][3:]) in range(1, 100)):
                    result["campaign"].append(CheckResult(ResultCode.ERROR,
                                                          "Global attribute 'campaign': If IOP then string must be IOPxx"))
            elif self.attrs["campaign"][0:4] in ["VALR", "VALM"]:
                is_lto = True
                if (len(self.attrs["campaign"]) != 6) or (not int(self.attrs["campaign"][4:]) in range(1, 100)):
                    result["campaign"].append(CheckResult(ResultCode.ERROR,
                                                          "Global attribute 'campaign': If VALM/VALR then string must be VALMxx/VALRxx"))

        ###
        # Check dims
        ###

        # TODO: There must not be any UNLIMITED dimension.
        if "nv" in self.dims.keys():
            if self.dims["nv"] != 2:
                result["nv_is_2"] = CheckResult(ResultCode.ERROR, "Dimension 'nv' must have size of 2.")
        if "max_name_len" in self.dims.keys():
            if self.dims["max_name_len"] != 32:
                result["max_name_len_is_32"] = CheckResult(ResultCode.ERROR, "Dimension 'max_name_len' must have size of 32.")

        ###
        # Check variables
        ###

        # vrs

        result["vrs"] = dict()
        result["vrs"]["variable"] = self.check_var("vrs", True, dims=())
        if result["vrs"]["variable"]:
            result["vrs"]["long_name"] = self.check_var_attr("vrs", "long_name", True, allowed_types=str,
                                                             allowed_values="vertical reference system")
            result["vrs"]["system_name"] = self.check_var_attr("vrs", "system_name", True, allowed_types=str, allowed_values="DHHN2016")
            result["vrs"]["standard_name"] = self.check_var_attr("vrs", "standard_name", False, must_not_exist=True)

        # time

        allowed_range = None
        if is_ts or is_tsp:
            time_dims = ("station", "ntime")
            time_bounds_dims = ("station", "ntime", "nv")
            time_dim_name = "ntime"
        elif is_traj:
            time_dims = ("trag", "ntime")
            time_bounds_dims = ("traj", "ntime", "nv")
            time_dim_name = "ntime"
        else:
            time_dims = ("time")
            time_bounds_dims = ("time", "nv")
            time_dim_name = "time"

        if is_iop:
            allowed_range = [.01, 86400]
        elif is_lto:
            if result["origin_time"]:
                ndays = calendar.monthrange(int(self.attrs["origin_time"][0:4]), int(self.attrs["origin_time"][5:7]))[1]
                allowed_range = [.01, ndays * 24 * 60 * 60]

        result["time"] = dict()
        result["time"]["variable"] = self.check_var("time", True, allowed_types=[numpy.int16, numpy.int32, numpy.float],
                                                    allowed_range=allowed_range, dims=time_dims,
                                                    must_be_sorted_along=time_dim_name,
                                                    decrease_sort_allowed=False,
                                                    fill_allowed=not is_grid)  # TODO: If coordinate var, then no missing values allowed.
        if result["time"]["variable"]:
            result["time"]["long_name"] = self.check_var_attr("time", "long_name", True, allowed_types=str, allowed_values="time")
            result["time"]["standard_name"] = self.check_var_attr("time", "standard_name", True, allowed_types=str,
                                                                  allowed_values="time")
            result["time"]["calendar"] = self.check_var_attr("time", "calendar", True, allowed_types=str,
                                                             allowed_values="proleptic_gregorian")
            result["time"]["axis"] = self.check_var_attr("time", "axis", True, allowed_types=str, allowed_values="T")
            result["time"]["fill_values"] = self.check_var_attr("time", "_FillValue", False, allowed_types=self.variables["time"].dtype,
                                                                must_not_exist=is_grid)
            result["time"]["units"] = self.check_var_attr("time", "units", True, allowed_types=str,
                                                          regex="seconds since [0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00")
            if result["origin_time"] and result["time"]["units"]:
                if self.attrs["origin_time"] != self.variables["time"].attrs["units"][14:]:
                    result["time"]["origin_time"] = CheckResult(ResultCode.ERROR,
                                                                "Global attribute 'origin_time' does not match units of variable 'time'.")
        result["time_bounds"] = dict()
        if "time_bounds" in self.variables.keys():
            result["time"]["bounds"] = self.check_var_attr("time", "bounds", True, allowed_types=str, allowed_values="time_bounds")
            if result["time"]["variable"]:
                result["time_bounds"]["variable"] = self.check_var("time_bounds", True, allowed_types=self.variables["time"].dtype,
                                                                   dims=time_bounds_dims)
                if not self.variables["time"][0].equals(self.variables["time_bounds"][0,:,1]):
                    result["time_bounds"]["variable"].append(CheckResult(ResultCode.ERROR, "second column of 'time_bounds' must equal data of variable 'time'"))

                if len(self.variables["time_bounds"].attrs) != 0:
                    result["time_bounds"]["attributes"] = CheckResult(ResultCode.ERROR, "Variable 'time_bounds' must not have any attributes.")
            else:
                result["time_bounds"]["variable"] = CheckResult(ResultCode.ERROR, "Cannot check variable 'time_bounds' because of error in variable 'time'")
        else:
            result["time"]["bounds"] = self.check_var_attr("time", "bounds", False, must_not_exist=True)

        # z

        if is_grid:
            z_dims = ("z")
            z_bounds_dims = ("z", "nv")
            must_be_sorted_along = "z"
        elif is_ts:
            z_dims = ("station")
            z_bounds_dims = ("station", "nv")
            must_be_sorted_along = None
        elif is_tsp:
            z_dims = ("station", "ntime", "nz")
            z_bounds_dims = ("station", "ntime", "nz", "nv")
            must_be_sorted_along = "nz"
        elif is_traj:
            z_dims = ("traj", "ntime")
            z_bounds_dims = ("traj", "ntime", "nv")
            must_be_sorted_along = None
        else: raise Exception("unexpected featureType.")

        result["z"] = dict()

        result["z"]["variable"] = self.check_var("z", True, allowed_types=[numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.float,
                                                                           numpy.float16, numpy.float32, numpy.float64], dims=z_dims,
                                                 must_be_sorted_along=must_be_sorted_along,
                                                 fill_allowed=not is_grid)
        result["z"]["long_name"] = self.check_var_attr("z", "long_name", True, allowed_types=str, allowed_values="height above origin")
        result["z"]["axis"] = self.check_var_attr("z", "axis", True, allowed_types=str, allowed_values="Z")
        result["z"]["positive"] = self.check_var_attr("z", "positive", True, allowed_types=str, allowed_values="up")
        if result["origin_z"]:
            result["z"]["standard_name"] = self.check_var_attr("z", "standard_name", self.attrs["origin_z"] == 0, allowed_types=str,
                                                               allowed_values="height_above_mean_sea_level",
                                                               must_not_exist=self.attrs["origin_z"] != 0)
        result["z_bounds"] = dict()
        if "z_bounds" in self.variables.keys():
            result["z"]["bounds"] = self.check_var_attr("z", "bounds", True, allowed_types=str, allowed_values="z_bounds")
            if result["z"]["variable"]:
                result["z_bounds"]["variable"] = self.check_var("z_bounds", True, allowed_types=self.variables["z"].dtype, dims=z_bounds_dims)

                z_bound_lower = self.variables["z_bounds"][0,:,0]
                z_bound_upper = self.variables["z_bounds"][0,:,1]
                z_bound_mid = z_bound_lower + (z_bound_upper - z_bound_lower) * 0.5
                if not numpy.allclose(self.variables["z"][0].values, z_bound_mid.values, equal_nan=True):
                    result["z_bounds"]["variable"].append(CheckResult(ResultCode.ERROR, "values of z must be in the middle between z_bounds."))

                if len(self.variables["z_bounds"].attrs) != 0:
                    result["z_bounds"]["attributes"] = CheckResult(ResultCode.ERROR, "Variable 'z_bounds' must not have any attributes.")
            else:
                result["z_bounds"]["variable"] = CheckResult(ResultCode.ERROR, "Cannot check variable 'z_bounds' because of error in variable 'z'")
        else:
            result["z"]["bounds"] = self.check_var_attr("z", "bounds", False, must_not_exist=True)

        if is_ts or is_tsp:
            result["station_h"] = self.check_var("station_h", True,
                                                 allowed_types=[numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.float,
                                                                numpy.float16, numpy.float32, numpy.float64], dims=("station"))

        # x, y
        result["x"] = self.check_xy("x")
        result["y"] = self.check_xy("y")
        result["lon"] = self.check_xy("lon")
        result["lat"] = self.check_xy("lat")
        result["E_UTM"] = self.check_xy("E_UTM")
        result["N_UTM"] = self.check_xy("N_UTM")
        # TODO: check xu, yv, zw, Eu_UTM, Nu_UTM, Ev_UTM, Nv_UTM, lonu, latu, lonv, latv

        # crs
        result["crs"] = dict()
        result["crs"]["variable"] = self.check_var("crs", True)
        result["crs"]["standard_name"] = self.check_var_attr("crs", "standard_name", False, must_not_exist=True)
        result["crs"]["long_name"] = self.check_var_attr("crs", "long_name", True,
                                                         allowed_values="coordinate reference system")
        result["crs"]["grid_mapping_name"] = self.check_var_attr("crs", "grid_mapping_name", True,
                                                                 allowed_values="transverse_mercator")
        result["crs"]["semi_major_axis"] = self.check_var_attr("crs", "semi_major_axis", True,
                                                               allowed_values=6378137)
        result["crs"]["inverse_flattening"] = self.check_var_attr("crs", "inverse_flattening", True,
                                                                  allowed_range=[298.2572,298.2573])
        result["crs"]["longitude_of_prime_meridian"] = self.check_var_attr("crs", "longitude_of_prime_meridian",
                                                                           True, allowed_values=0)
        result["crs"]["longitude_of_central_meridian"] = self.check_var_attr("crs", "longitude_of_central_meridian",
                                                                             True, allowed_values=[3,9,15])
        result["crs"]["scale_factor_at_central_meridian"] = self.check_var_attr("crs", "scale_factor_at_central_meridian",
                                                                                True, allowed_range=[0.9995,0.9997])
        result["crs"]["latitude_of_projection_origin"] = self.check_var_attr("crs", "latitude_of_projection_origin",
                                                                             True, allowed_values=0)
        result["crs"]["false_easting"] = self.check_var_attr("crs", "false_easting", True, allowed_values=500000)
        result["crs"]["false_northing"] = self.check_var_attr("crs", "false_northing", True, allowed_values=0)
        result["crs"]["units"] = self.check_var_attr("crs", "units", True, allowed_values="m")
        result["crs"]["epsg_code"] = self.check_var_attr("crs", "epsg_code", True,
                                                         allowed_values=["EPSG:25831","EPSG:25832","EPSG:25833"])

        #
        # other (auxiliary) coordinate variables
        #

        check_platform = False
        if is_ts or is_tsp:
            check_platform = True
            name = "station_name"
            long_name = "station name"
            dim = "station"
            id = "timeseries_id"
        elif is_traj:
            check_platform = True
            name = "traj_name"
            long_name = "trajectory name"
            dim = "traj"
            id = "trajectory_id"

        if check_platform:
            result[name] = dict()
            result[name]["variable"] = self.check_var(name, True, allowed_types=numpy.dtype("S1"),
                                                                dims=(dim, "max_name_len"))
            if result[name]["variable"]:
                result[name]["long_name"] = self.check_var_attr(name, "long_name", True, allowed_types=str,
                                                                allowed_values=long_name)
                result[name]["standard_name"] = self.check_var_attr(name, "standard_name", True, allowed_types=str,
                                                                    allowed_values="platform_name")
                result[name]["cf_role"] = self.check_var_attr(name, "cf_role", True, allowed_types=str,
                                                              allowed_values=id)

        if is_ts or is_tsp:
            result["station_h"] = dict()
            result["station_h"]["variable"] = self.check_var("station_h", True,
                                                             allowed_types=[numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.float,
                                                             numpy.float16, numpy.float32, numpy.float64], dims="station",
                                                             fill_allowed=False)
            if result["station_h"]["variable"]:
                result["station_h"]["long_name"] = self.check_var_attr("station_h", "long_name", True, allowed_types=str,
                                                                       allowed_values="surface altitude")
                result["station_h"]["standard_name"] = self.check_var_attr("station_h", "standard_name", True, allowed_types=str,
                                                                       allowed_values="surface_altitude")
                result["station_h"]["units"] = self.check_var_attr("station_h", "units", True,
                                                                   allowed_types=str,
                                                                   allowed_values="m")
        if is_traj:
            result["height"] = dict()
            result["height"]["variable"] = self.check_var("height", True, allowed_types=[numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.float,
                                                             numpy.float16, numpy.float32, numpy.float64])
            if result["height"]["variable"]:
                if self.variables["height"].dims != () or self.variables["height"].dims != ("traj", "ntime"):
                    result["height"]["variable"] = CheckResult(ResultCode.ERROR, "Variable 'height' must either be scalar "+
                                                               "or have dimensions (traj, ntime).")
            if result["height"]["variable"]:
                result["height"]["long_name"] = self.check_var_attr("height", "long_name", True,
                                                                  allowed_types=str, allowed_values="height above surface")
                result["height"]["standard_name"] = self.check_var_attr("height", "standard_name", True,
                                                                      allowed_types=str, allowed_values="height")
                result["height"]["units"] = self.check_var_attr("height", "units", True, allowed_types=str,
                                                                allowed_values="m")

        ###
        # Data variables
        ###
        dv = dict()
        for ikey in self.data_vars:
            if not ikey in ["station_name", "station_h", "z", "z_bounds",
                            "crs", "vrs", "x", "xu", "xs", "y", "yv", "ys", "z", "zw", "zs",
                            "lon", "lonu", "lonv", "lons", "lat", "latu", "latv", "lats",
                            "E_UTM", "Eu_UTM", "Ev_UTM", "Es_UTM", "N_UTM", "Nu_UTM", "Nv_UTM", "Ns_UTM",
                            "s",
                            "time", "time_bounds",
                            "azimuth", "azimuths", "zenith", "zeniths"]:
                dv[ikey] = self.data_vars[ikey]

        for ikey in dv:
            if ikey in self.allowed_variables:
                result[ikey] = dict()
                result[ikey]["variable"] = self.check_var()
            else:
                pass
            # TODO: Check for extension by agg_method
            # TODO: Check agg_method with cell_methods
            # TODO: Check flags, ancillaries and spectra
            pass

        # TODO: "all coordinate and auxiliary coordinate variables must be specified in the attribute coordinates of the respective data variable"
        # TODO: If only one variable, data_content MUST be variable name.
        # TODO: If all variables have cell_methods with time:point then no time_bounds (and bounds attribute)
        # TODO: If all variables have cell_methods with z:point then no z_bounds (and bounds attribute)
        # TODO: Need height variable for "trajectory"
        # TODO: Check azimuth, zenith
        # TODO: Check ranges of allowed values and "almost_equal"

        ###
        # TODO: Check geo between var and glob att
        ###


        return result

    def check_xy(self, xy):
        if not xy in ["x", "y", "lon", "lat", "E_UTM", "N_UTM"]:
            raise Exception('Unexpected variable: ' + xy)

        if self.featuretype == "None":
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
        elif self.featuretype in ["timeSeries", "timeSeriesProfile"]:
            dims = "station"
            sort_along = None
            fill_allowed = True
        elif self.featuretype == "trajectory":
            dims = ("traj", "ntime")
            sort_along = None
            fill_allowed = True
        else:
            raise Exception("Unexpected featureType")
            # TODO: In case of "pixel-based surfaces" x has dimensions (time, nrow, ncol)

        out = dict()
        out["variable"] = self.check_var(xy, True,
                                         allowed_types=[numpy.float, numpy.float16, numpy.float32, numpy.float64,
                                                        numpy.int, numpy.int8, numpy.int16, numpy.int32],
                                         dims=dims, must_be_sorted_along=sort_along, decrease_sort_allowed=True,
                                         fill_allowed=fill_allowed)
        if out["variable"]:

            if xy in ["x", "y"]:
                long_n = "distance to origin in "+xy+"-direction"
                standard_n = None
                axis = xy.upper()
                units = "m"
            elif xy in ["lon", "lat"]:
                long_n = "longitude" if xy == "lon" else "latitude"
                standard_n = long_n
                axis = None
                units = "degrees_east" if xy == "lon" else "degrees_north"
            elif xy in ["E_UTM", "N_UTM"]:
                long_n = "easting" if xy == "E_UTM" else "northing"
                standard_n = "projection_x_coordinate" if xy == "E_UTM" else "projection_y_coordinate"
                axis = None
                units = "m"

            out["standard_name"] = self.check_var_attr(xy, "standard_name", not xy in ["x", "y"],
                                                       must_not_exist=xy in ["x", "y"], allowed_values=standard_n)
            out["long_name"] = self.check_var_attr(xy, "long_name", True, allowed_types=str,
                                                   allowed_values=long_n)
            out["units"] = self.check_var_attr(xy, "units", True, allowed_types=str,
                                               allowed_values=units)
            out["axis"] = self.check_var_attr(xy, "axis", axis is not None, allowed_types=str, allowed_values=axis)

        return out

    def check_var(self, varname, must_exist, allowed_types=None, allowed_range=None, dims=None,
                  must_be_sorted_along=None, decrease_sort_allowed=True, fill_allowed=True):
        exists = varname in self.variables.keys()
        result = CheckResult(ResultCode.OK, "Test passed.")

        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR, "Required variable '" + varname + "' not found.")
            else:
                return result

        if allowed_types:
            if type(allowed_types) == numpy.dtype:
                allowed_types = [allowed_types]
            if not self.variables[varname].dtype in allowed_types:
                result.append(CheckResult(ResultCode.ERROR, "Variable '" + varname + "' has wrong type. Should be " +
                                          "one of the following: " + str(allowed_types)))

        if allowed_range:
            if (self.variables[varname].min() < allowed_range[0]) or (self.variables[varname].max() > allowed_range[1]):
                result.append(CheckResult(ResultCode.ERROR,
                                          "Variable '" + varname + "' is outside allowed range" + str(allowed_range)))

        if dims:
            if type(dims) == list:
                dims = tuple(dims)
            elif type(dims) == str:
                dims = tuple([dims])
            if self.variables[varname].dims != dims:
                result.append(CheckResult(ResultCode.ERROR, "Variable '" + varname + "' has wrong dimensions. Expected: " + str(dims)))

        if must_be_sorted_along:
            if must_be_sorted_along in self.variables[varname].dims:
                sorted_arr = numpy.sort(self.variables[varname], axis=self.variables[varname].dims.index(must_be_sorted_along))
                if decrease_sort_allowed:
                    anti_sorted_arr = -numpy.sort(-self.variables[varname], axis=self.variables[varname].dims.index(must_be_sorted_along))
                    if not (numpy.array_equal(self.variables[varname], sorted_arr) or numpy.array_equal(self.variables[varname], anti_sorted_arr)):
                        result.append(CheckResult(ResultCode.ERROR, "Variable '" + varname + "' must be sorted along dimension '" + must_be_sorted_along + "'"))
                else:
                    if not numpy.array_equal(self.variables[varname], sorted_arr):
                        result.append(CheckResult(ResultCode.ERROR, "Variable '" + varname + "' must be sorted along dimension '" + must_be_sorted_along + "'"))

        if fill_allowed:
            pass # TODO: do it.

        return result


    def check_var_attr(self, varname, attrname, must_exist, allowed_types=None, allowed_values=None, regex=None,
                       must_not_exist=None, allowed_range=None):
        exists = attrname in self.variables[varname].attrs.keys()
        result = CheckResult(ResultCode.OK, "Test passed.")
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

        if allowed_types:
            if type(allowed_types) == type:
                allowed_types = [allowed_types]
            if not type(self.variables[varname].attrs[attrname]) in allowed_types:
                result.append(CheckResult(ResultCode.ERROR,
                                          "Variable '" + varname + "': Required variable attribute '" + attrname + "' has wrong type. Should be " +
                                          "one of the following: " + str(allowed_types)))

        if allowed_values:
            if type(allowed_values) != list:
                allowed_values = [allowed_values]
            if not self.variables[varname].attrs[attrname] in allowed_values:
                if len(allowed_values) == 1:
                    result.append(CheckResult(ResultCode.ERROR,
                                              "Variable '" + varname + "': Required variable attribute '" + attrname + "'  has wrong value. Should be " +
                                              str(allowed_values[0])))
                else:
                    result.append(CheckResult(ResultCode.ERROR,
                                              "Variable '" + varname + "': Required variable attribute '" + attrname + "' has wrong value"))

        if allowed_range:
            if self.variables[varname].attrs[attrname] < allowed_range[0] or \
                    self.variables[varname].attrs[attrname] > allowed_range[1]:
                result.append(CheckResult(ResultCode.ERROR,
                                          "Variable '" + varname + "': Attribute '"+attrname+"' outside range. Expected: "+
                                          str(allowed_range)))

        if regex:
            if re.fullmatch(regex, self.variables[varname].attrs[attrname]) is None:
                result.append(CheckResult(ResultCode.ERROR,
                                          "Global attribute '" + attrname + "' does not match regular expression " + regex))
        return result

    def check_glob_attr(self, attrname, must_exist, allowed_types=None, allowed_values=None,
                        max_strlen=None, regex=None, allowed_range=None):
        exists = attrname in self.attrs.keys()
        result = CheckResult(ResultCode.OK, "Test passed.")
        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR, "Required global attribute '" + attrname + "' not found.")
            else:
                return result

        if allowed_types:
            if type(allowed_types) == type:
                allowed_types = [allowed_types]
            if not type(self.attrs[attrname]) in allowed_types:
                result.append(CheckResult(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong type. Should be " +
                                          "one of the following: " + str(allowed_types)))

        if allowed_values:
            if not self.attrs[attrname] in allowed_values:
                if len(allowed_values) == 1:
                    result.append(CheckResult(ResultCode.ERROR,
                                              "Global attribute '" + attrname + "' has wrong value. Should be " +
                                              str(allowed_values[0])))
                else:
                    result.append(CheckResult(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong value"))

        if regex:
            if re.fullmatch(regex, self.attrs[attrname]) is None:
                result.append(CheckResult(ResultCode.ERROR,
                                          "Global attribute '" + attrname + "' does not match regular expression " + regex))

        if max_strlen:
            if len(self.attrs[attrname]) > max_strlen:
                result.append(CheckResult(ResultCode.ERROR,
                                          "Global attribute '" + attrname + "' is too long. Must be max. "+
                                          str(max_strlen) + " characters."))

        if allowed_range:
            if (self.attrs[attrname] < allowed_range[0]) or (self.attrs[attrname] > allowed_range[1]):
                return CheckResult(ResultCode.ERROR,
                                   "Global attribute '" + attrname + "' is outside allowed range " + str(allowed_range))

        return result

    def get_filename(self):
        attrs = ["campaign", "location", "site", "acronym", "data_content", "origin_time", "version"]
        vals = list()
        for i in attrs:
            if not self.check_result[i]:
                raise Exception("Cannot parse filename. Global attribute '" + i + "' did not pass UC2 conformity test.")

            if i == "origin_time":
                vals.append(self.attrs[i][: 10].replace("-", ""))
            elif i == "version":
                vals.append(str(self.attrs[i]).zfill(3))
            else:
                vals.append(self.attrs[i].replace("-", "_"))
        return "_".join(vals) + ".nc"


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
    return CheckResult()


class ResultCode(enum.Enum):
    OK = 1
    WARNING = 2
    ERROR = 3


class CheckResult:

    def __init__(self, result=ResultCode.OK, message="Test passed."):
        self.message = [message]
        self.result = [result]

    def __bool__(self):
        if  ResultCode.ERROR in self.result:
            return False
        return True

    def append(self, value):
        if value.result == ResultCode.OK:
            return
        else:
            while ResultCode.OK in self.result:
                idx = self.result.index(ResultCode.OK)
                self.result.pop(idx)
                self.message.pop(idx)
        self.result += value.result
        self.message += value.message


def printCheckResult(value):
    out = list()
    for k in value.keys():
        if isinstance(value[k], dict):
            out.append(k)
            out.extend(list(str_add("  ", printCheckResult(value[k]))))
        else:
            out.append(k)
            for idx, i in enumerate(value[k].result):
                out.extend(list(str_add("  ", [value[k].message[idx] + ' (' + str(value[k].result[idx]) + ')'])))
    return out