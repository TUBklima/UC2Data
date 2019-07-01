from __future__ import annotations
import xarray
import numpy
import sys
import enum
import csv
import re
import calendar
from numpy.core.defchararray import add as str_add
from collections import OrderedDict
from typing import Union

is_win = sys.platform in ['win32', 'win64']
if not is_win:
    from cfchecker import cfchecks


aggregations_file = "aggregations.txt"
data_content_file = "data_content.txt"
variables_file = "variables.txt"
institutions_file = "institutions.txt"
sites_file = "sites.txt"


class UC2Data(xarray.Dataset):

    all_floats = [float, numpy.float, numpy.float16, numpy.float32, numpy.float64]
    all_ints = [int, numpy.int, numpy.int8, numpy.int16, numpy.int32, numpy.int64]

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

        result = CheckResult()

        result["title"].add(self.check_glob_attr("title", True, str))
        result["data_content"].add(self.check_glob_attr("data_content", True, str,
                                                        allowed_values=UC2Data.allowed_data_contents,
                                                        max_strlen=16))  # TODO: Redo this test when variable is checked
        result["source"].add(self.check_glob_attr("source", True, str))
        result["version"].add(self.check_glob_attr("version", True,
                                                   int,
                                                   allowed_values=list(
                                                       range(1, 1000))))  # TODO: This is going to be checked in DMS
        result["Conventions"].add(self.check_glob_attr("Conventions", True, str, allowed_values=["CF-1.7"]))
        result["dependencies"].add(self.check_glob_attr("dependencies", True,
                                                        str))  # TODO: This is going to be checked by DMS
        result["history"].add(self.check_glob_attr("history", True, str))
        result["references"].add(self.check_glob_attr("references", True, str))
        result["comment"].add(self.check_glob_attr("comment", True, str))
        result["keywords"].add(self.check_glob_attr("keywords", True, str))
        result["licence"].add(self.check_glob_attr("licence", True, str, allowed_values=UC2Data.allowed_licences))
        result["creation_time"].add(self.check_glob_attr("creation_time", True, str,
                                                         regex="[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00"))
        result["origin_time"].add(self.check_glob_attr("origin_time", True, str,
                                                       regex="[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00"))  # TODO: Check later with time units.
        result["origin_lon"].add(self.check_glob_attr("origin_lon", True, float,
                                                      allowed_range=[-180, 180]))
        result["origin_lat"].add(self.check_glob_attr("origin_lat", True, float,
                                                      allowed_range=[-90, 90]))
        result["origin_x"].add(self.check_glob_attr("origin_x", True, float))
        result["origin_y"].add(self.check_glob_attr("origin_y", True, float))
        result["rotation_angle"].add(self.check_glob_attr("rotation_angle", True,
                                                          float,
                                                          allowed_range=[0, 360]))

        # non-standard checks

        if not is_grid:
            result["featureType"].add(self.check_glob_attr("featureType", False, str,
                                                           allowed_values=["timeSeries", "timeSeriesProfile",
                                                                           "trajectory"]))
            if not result["featureType"]:
                return result

        result["origin_z"].add(self.check_glob_attr("origin_z", True, float,
                                                    allowed_values=0 if (not is_grid) else None))
        result["location"].add(self.check_glob_attr("location", True, str, allowed_values=UC2Data.allowed_locations))
        result["site"].add(self.check_glob_attr("site", True, str, allowed_values=UC2Data.allowed_sites,
                                                max_strlen=12))
        if result["location"] and result["site"]:
            if UC2Data.allowed_locations[UC2Data.allowed_sites.index(self.attrs["site"])] != self.attrs["location"]:
                result["site"].add(ResultCode.ERROR, "site '" + self.attrs[
                    "site"] + "' does not match location '" + self.attrs["location"] + "'")
                result["location"].add(ResultCode.ERROR, "site '" + self.attrs[
                    "site"] + "' does not match location '" + self.attrs["location"] + "'")

        result["institution"].add(
            self.check_glob_attr("institution", True, str, allowed_values=UC2Data.allowed_institutions))
        result["acronym"].add(self.check_glob_attr("acronym", True, str, allowed_values=UC2Data.allowed_acronyms,
                                                   max_strlen=12))
        if result["institution"] and result["acronym"]:
            if UC2Data.allowed_institutions.index(self.attrs["institution"]) != UC2Data.allowed_acronyms.index(
                    self.attrs["acronym"]):
                result["institution"].add(ResultCode.ERROR, "institution '" + self.attrs[
                    "institution"] + "' does not match acronym '" + self.attrs["acronym"] + "'")
                result["acronym"].add(ResultCode.ERROR, "institution '" + self.attrs[
                    "institution"] + "' does not match acronym '" + self.attrs["acronym"] + "'")

        result["author"].add(self.check_glob_attr("author", True, str))
        if result["author"]:
            if self.attrs["author"] != "":
                result["author"].add(check_person_field(self.attrs["author"], "author"))

        result["contact_person"].add(self.check_glob_attr("contact_person", True, str))
        if result["contact_person"]:
            result["contact_person"].add(check_person_field(self.attrs["contact_person"], "contact_person"))

        is_iop = False
        is_lto = False
        result["campaign"].add(self.check_glob_attr("campaign", True, str, regex="^[A-Za-z0-9\._-]+$",
                                                    max_strlen=12))
        if result["campaign"]:
            if self.attrs["campaign"][0:3] == "IOP":
                is_iop = True
                if (len(self.attrs["campaign"]) != 5) or (not int(self.attrs["campaign"][3:]) in range(1, 100)):
                    result["campaign"].add(ResultCode.ERROR,
                                           "Global attribute 'campaign': If IOP then string must be IOPxx")
            elif self.attrs["campaign"][0:4] in ["VALR", "VALM"]:
                is_lto = True
                if (len(self.attrs["campaign"]) != 6) or (not int(self.attrs["campaign"][4:]) in range(1, 100)):
                    result["campaign"].add(ResultCode.ERROR,
                                           "Global attribute 'campaign': If VALM/VALR then string must be VALMxx/VALRxx")

        ###
        # Check dims
        ###

        # TODO: There must not be any UNLIMITED dimension.
        if "nv" in self.dims.keys():
            if self.dims["nv"] != 2:
                result["nv_is_2"].add(ResultCode.ERROR, "Dimension 'nv' must have size of 2.")
        if "max_name_len" in self.dims.keys():
            if self.dims["max_name_len"] != 32:
                result["max_name_len_is_32"].add(ResultCode.ERROR, "Dimension 'max_name_len' must have size of 32.")

        ###
        # Check variables
        ###

        # vrs

        result["vrs"]["variable"].add(self.check_var("vrs", True, dims=()))
        if result["vrs"]["variable"]:
            result["vrs"]["long_name"].add(self.check_var_attr("vrs", "long_name", True, allowed_types=str,
                                                               allowed_values="vertical reference system"))
            result["vrs"]["system_name"].add(
                self.check_var_attr("vrs", "system_name", True, allowed_types=str, allowed_values="DHHN2016"))
            result["vrs"]["standard_name"].add(self.check_var_attr("vrs", "standard_name", False, must_not_exist=True))

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

        result["time"]["variable"].add(
            self.check_var("time", True, allowed_types=[int, float],
                           allowed_range=allowed_range, dims=time_dims,
                           must_be_sorted_along=time_dim_name,
                           decrease_sort_allowed=False,
                           fill_allowed=not is_grid))  # TODO: If coordinate var, then no missing values allowed.
        if result["time"]["variable"]:
            # bounds are checked below together with other variables.
            result["time"]["long_name"].add(
                self.check_var_attr("time", "long_name", True, allowed_types=str, allowed_values="time"))
            result["time"]["standard_name"].add(self.check_var_attr("time", "standard_name", True, allowed_types=str,
                                                                    allowed_values="time"))
            result["time"]["calendar"].add(self.check_var_attr("time", "calendar", True, allowed_types=str,
                                                               allowed_values="proleptic_gregorian"))
            result["time"]["axis"].add(self.check_var_attr("time", "axis", True, allowed_types=str, allowed_values="T"))
            result["time"]["fill_values"].add(
                self.check_var_attr("time", "_FillValue", False, allowed_types=self["time"].dtype,
                                    must_not_exist=is_grid))
            result["time"]["units"].add(self.check_var_attr("time", "units", True, allowed_types=str,
                                                            regex="seconds since [0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} \+00"))
            if result["origin_time"] and result["time"]["units"]:
                if self.attrs["origin_time"] != self["time"].units[14:]:
                    result["time"]["origin_time"].add(ResultCode.ERROR,
                                                      "Global attribute 'origin_time' does not match units of variable 'time'.")

        # z

        if is_grid:
            z_dims = ("z")
            must_be_sorted_along = "z"
        elif is_ts:
            z_dims = ("station")
            must_be_sorted_along = None
        elif is_tsp:
            z_dims = ("station", "ntime", "nz")
            must_be_sorted_along = "nz"
        elif is_traj:
            z_dims = ("traj", "ntime")
            must_be_sorted_along = None
        else:
            raise Exception("unexpected featureType.")

        result["z"]["variable"].add(
            self.check_var("z", True, allowed_types=[int, float], dims=z_dims,
                           must_be_sorted_along=must_be_sorted_along,
                           fill_allowed=not is_grid))
        result["z"]["long_name"].add(
            self.check_var_attr("z", "long_name", True, allowed_types=str, allowed_values="height above origin"))
        result["z"]["axis"].add(self.check_var_attr("z", "axis", True, allowed_types=str, allowed_values="Z"))
        result["z"]["positive"].add(self.check_var_attr("z", "positive", True, allowed_types=str, allowed_values="up"))
        # Bounds will be checked below with all other variables.
        if result["z"]:
            if result["origin_z"]:
                result["z"]["standard_name"].add(
                    self.check_var_attr("z", "standard_name", self.attrs["origin_z"] == 0, allowed_types=str,
                                        allowed_values="height_above_mean_sea_level",
                                        must_not_exist=self.attrs["origin_z"] != 0))

        if is_ts or is_tsp:
            result["station_h"].add(self.check_var("station_h", True,
                                                   allowed_types=[int, float],
                                                   dims=("station")))

        # x, y
        result["x"].add(self.check_xy("x"))
        result["y"].add(self.check_xy("y"))
        result["lon"].add(self.check_xy("lon"))
        result["lat"].add(self.check_xy("lat"))
        result["E_UTM"].add(self.check_xy("E_UTM"))
        result["N_UTM"].add(self.check_xy("N_UTM"))

        if "s" in self.dims:
            result["xs"].add(self.check_xy("xs"))
            result["ys"].add(self.check_xy("ys"))
            result["lons"].add(self.check_xy("lons"))
            result["lats"].add(self.check_xy("lats"))
            result["Es_UTM"].add(self.check_xy("Es_UTM"))
            result["Ns_UTM"].add(self.check_xy("Ns_UTM"))

        if is_grid:
            # if one u is there, all are nedded
            if any(elem in self for elem in ["xu", "Eu_UTM", "Nu_UTM", "lonu", "latu"]):
                result["xu"].add(self.check_xy("xu"))
                result["Eu_UTM"].add(self.check_xy("Eu_UTM"))
                result["Nu_UTM"].add(self.check_xy("Nu_UTM"))
                result["latu"].add(self.check_xy("latu"))
                result["lonu"].add(self.check_xy("lonu"))
            # if one v is there, all are nedded
            if any(elem in self for elem in ["xv", "Ev_UTM", "Nv_UTM", "lonv", "latv"]):
                result["yv"].add(self.check_xy("yv"))
                result["Ev_UTM"].add(self.check_xy("Ev_UTM"))
                result["Nv_UTM"].add(self.check_xy("Nv_UTM"))
                result["lonv"].add(self.check_xy("lonv"))
                result["latv"].add(self.check_xy("latv"))

        # crs
        result["crs"]["variable"].add(self.check_var("crs", True))
        result["crs"]["standard_name"].add(self.check_var_attr("crs", "standard_name", False, must_not_exist=True))
        result["crs"]["long_name"].add(self.check_var_attr("crs", "long_name", True,
                                                           allowed_values="coordinate reference system"))
        result["crs"]["grid_mapping_name"].add(self.check_var_attr("crs", "grid_mapping_name", True,
                                                                   allowed_values="transverse_mercator"))
        result["crs"]["semi_major_axis"].add(self.check_var_attr("crs", "semi_major_axis", True,
                                                                 allowed_values=6378137))
        result["crs"]["inverse_flattening"].add(self.check_var_attr("crs", "inverse_flattening", True,
                                                                    allowed_range=[298.2572, 298.2573]))
        result["crs"]["longitude_of_prime_meridian"].add(self.check_var_attr("crs", "longitude_of_prime_meridian",
                                                                             True, allowed_values=0))
        result["crs"]["longitude_of_central_meridian"].add(self.check_var_attr("crs", "longitude_of_central_meridian",
                                                                               True, allowed_values=[3, 9, 15]))
        result["crs"]["scale_factor_at_central_meridian"].add(
            self.check_var_attr("crs", "scale_factor_at_central_meridian",
                                True, allowed_range=[0.9995, 0.9997]))
        result["crs"]["latitude_of_projection_origin"].add(self.check_var_attr("crs", "latitude_of_projection_origin",
                                                                               True, allowed_values=0))
        result["crs"]["false_easting"].add(self.check_var_attr("crs", "false_easting", True, allowed_values=500000))
        result["crs"]["false_northing"].add(self.check_var_attr("crs", "false_northing", True, allowed_values=0))
        result["crs"]["units"].add(self.check_var_attr("crs", "units", True, allowed_values="m"))
        result["crs"]["epsg_code"].add(self.check_var_attr("crs", "epsg_code", True,
                                                           allowed_values=["EPSG:25831", "EPSG:25832", "EPSG:25833"]))

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
            result[name]["variable"].add(self.check_var(name, True, allowed_types=numpy.dtype("S1"),
                                                        dims=(dim, "max_name_len")))
            if result[name]["variable"]:
                result[name]["long_name"].add(self.check_var_attr(name, "long_name", True, allowed_types=str,
                                                                  allowed_values=long_name))
                result[name]["standard_name"].add(self.check_var_attr(name, "standard_name", True, allowed_types=str,
                                                                      allowed_values="platform_name"))
                result[name]["cf_role"].add(self.check_var_attr(name, "cf_role", True, allowed_types=str,
                                                                allowed_values=id))

        if is_ts or is_tsp:
            result["station_h"]["variable"].add(self.check_var("station_h", True,
                                                               allowed_types=[int, float], dims="station",
                                                               fill_allowed=False))
            if result["station_h"]["variable"]:
                result["station_h"]["long_name"].add(
                    self.check_var_attr("station_h", "long_name", True, allowed_types=str,
                                        allowed_values="surface altitude"))
                result["station_h"]["standard_name"].add(
                    self.check_var_attr("station_h", "standard_name", True, allowed_types=str,
                                        allowed_values="surface_altitude"))
                result["station_h"]["units"].add(self.check_var_attr("station_h", "units", True,
                                                                     allowed_types=str,
                                                                     allowed_values="m"))
        if is_traj:
            result["height"]["variable"].add(self.check_var("height", True,
                                                            allowed_types=[int, float]))
            if result["height"]["variable"]:
                if self["height"].dims != () or self["height"].dims != ("traj", "ntime"):
                    result["height"]["variable"].add(ResultCode.ERROR, "Variable 'height' must either be scalar " +
                                                     "or have dimensions (traj, ntime).")
            if result["height"]["variable"]:
                result["height"]["long_name"].add(self.check_var_attr("height", "long_name", True,
                                                                      allowed_types=str,
                                                                      allowed_values="height above surface"))
                result["height"]["standard_name"].add(self.check_var_attr("height", "standard_name", True,
                                                                          allowed_types=str, allowed_values="height"))
                result["height"]["units"].add(self.check_var_attr("height", "units", True, allowed_types=str,
                                                                  allowed_values="m"))

        ###
        # Data variables
        ###

        dv = dict()
        data_content_var_names = list()
        known_coordinates = ["station_name", "station_h", "z", "z_bounds",
                             "crs", "vrs", "x", "xu", "xs", "y", "yv", "ys", "z", "zw", "zs",
                             "lon", "lonu", "lonv", "lons", "lat", "latu", "latv", "lats",
                             "E_UTM", "Eu_UTM", "Ev_UTM", "Es_UTM", "N_UTM", "Nu_UTM", "Nv_UTM", "Ns_UTM",
                             "s",
                             "time", "time_bounds",
                             "azimuth", "azimuths", "zenith", "zeniths"]
        if is_ts:
            data_dims = ("station", "ntime")
        elif is_tsp:
            data_dims = ("station", "ntime", "nz")
        elif is_traj:
            data_dims = ("traj", "ntime")
        else:
            data_dims = None

        for ikey in self.variables:

            is_normal = ikey in self.allowed_variables.keys()
            is_agg = ikey in [a + "_" + b for a in self.allowed_variables.keys() for b in self.allowed_aggregations.keys()]
            is_bounds = ikey.endswith("_bounds")
            is_bands = ikey.startswith("bands_")
            is_ancillary = ikey.startswith("ancillary_")
            is_coordinate = ikey in known_coordinates or is_bands

            if not any([is_normal, is_agg, is_bands, is_bounds, is_bands, is_ancillary, is_coordinate]):
                result[ikey].add(ResultCode.ERROR, "'" + ikey + "' is not a supported variable name.")

            if is_bands and not is_bounds:
                # Check bands (bands are coordinate variables => need dim of same name)
                result[ikey].add(self.check_var(ikey, True, dims=ikey, fill_allowed=False, must_be_sorted_along=ikey))
            else:
                if is_bounds:


                    main_key = ikey[:-7]
                    if main_key not in self:
                        result[ikey].add(ResultCode.ERROR,
                                         "Variable '" + ikey + "' seems to be a bounds variable " \
                                                               "but there is no main variable (expected '" + \
                                         main_key + "')")
                    else:
                        result[main_key].add(self.check_var_attr(main_key, "bounds", True,
                                                                 allowed_types=str, allowed_values=ikey))
                        result[ikey].add(self.check_var(ikey, True, allowed_types=self[main_key].dtype,
                                                        dims=self[main_key].dims + ("nv",)))
                        if len(self[ikey].attrs) != 0:
                            result[ikey]["attributes"].add(ResultCode.ERROR,
                                                           "Variable '" + ikey + "' must not have any attributes.")
                    # Time must be end of time period
                    if ikey == "time_bounds":
                        if not self[main_key][0].equals(self[ikey][0, :, 1]):
                            result[ikey]["variable"].add(ResultCode.ERROR,
                                                         "second column of 'time_bounds' must equal data of variable 'time'")
                    # z must be in middle of z bounds
                    if ikey == "z_bounds":
                        z_bound_lower = self[ikey][0, :, 0]
                        z_bound_upper = self[ikey][0, :, 1]
                        z_bound_mid = z_bound_lower + (z_bound_upper - z_bound_lower) * 0.5
                        if not numpy.allclose(self[main_key][0].values, z_bound_mid.values, equal_nan=True):
                            result[ikey]["variable"].add(ResultCode.ERROR,
                                                         "values of z must be in the middle between z_bounds.")


                elif is_ancillary:
                    # Check ancillary
                    pass
                elif is_coordinate:
                    # What to do?
                    pass
                else:
                    if is_normal:
                        expected_data_content = ikey
                    elif is_agg:
                        expected_data_content = "_".join(ikey[0].split("_")[:-1])
                    else:
                        raise Exception("unknown variable type")
                    if expected_data_content not in data_content_var_names:
                        data_content_var_names.append(expected_data_content)

                    # Check var
                    result[ikey]["variable"].add(self.check_var(ikey, True, dims=data_dims))

                    # Check obligatory attributes
                    result[ikey]["long_name"].add(self.check_var_attr(ikey, "long_name", True, allowed_types=str,
                                                                      allowed_values=self.allowed_variables[ikey]["long_name"]))
                    result[ikey]["units"].add(self.check_var_attr(ikey, "units", True, allowed_types=str)) # TODO: check conversion
                    result[ikey]["_FillValue"].add(self.check_var_attr(ikey, "_FillValue", True, allowed_types=self[ikey].dtype,
                                                                       allowed_values=-9999))
                    result[ikey]["coordinates"].add(self.check_var_attr(ikey, "coordinates", True, allowed_types=str)) # TODO: Check consistency
                    result[ikey]["grid_mapping"].add(self.check_var_attr(ikey, "grid_mapping", True, allowed_types=str,
                                                                         allowed_values="crs"))
                    # other attributes
                    result[ikey]["standard_name"].add(self.check_var_attr(ikey, "standard_name",
                                                                          self.allowed_variables[ikey]["standard_name"] != "",
                                                                          allowed_types=str,
                                                                          allowed_values=self.allowed_variables[ikey]["standard_name"],
                                                                          must_not_exist=self.allowed_variables[ikey]["standard_name"] == ""))
                    result[ikey]["units_alt"].add(self.check_var_attr(ikey, "units_alt", False, allowed_types=str)) # TODO: check conversion
                    result[ikey]["uncertainty_rel"].add(self.check_var_attr(ikey, "uncertainty_rel", False,
                                                                            allowed_types=float))
                    result[ikey]["processing_level"].add(self.check_var_attr(ikey, "processing_level", False,
                                                                             allowed_types=int, allowed_range=[0,3]))
                    result[ikey]["processing_info"].add(self.check_var_attr(ikey, "processing_info", False,
                                                                            allowed_types=str))
                    result[ikey]["instrument_name"].add(self.check_var_attr(ikey, "instrument_name", False,
                                                                            allowed_types=str))
                    result[ikey]["instrument_nr"].add(self.check_var_attr(ikey, "instrument_nr", False,
                                                                          allowed_types=str))
                    result[ikey]["instrument_sn"].add(self.check_var_attr(ikey, "instrument_sn", False,
                                                                          allowed_types=str))

                    # Check ancillary_variables attribute
                    if "ancillary_variables" in self[ikey].attrs.keys():
                        anc_var = self[ikey].attrs["ancillary_variables"].split(" ")
                        for i in anc_var:
                            if i not in self.keys():
                                result[ikey]["ancillary_variables"].add(ResultCode.ERROR,
                                                                        "Expected ancillary variable '"+
                                                                        i+"' not found in file.")
                        result[ikey]["ancillary_variables"].add(self.check_var_attr(ikey, "ancillary_variables",
                                                                                    True, allowed_types=str))


                    # Check bounds attribute
                    if "bounds" in self[ikey].attrs.keys():
                        if ikey + "_bounds" not in self.keys():
                            result[ikey]["bounds"].add(ResultCode.ERROR,
                                                       "Expected variable '" + ikey + "_bounds' not found.")
                        result[ikey]["bounds"].add(self.check_var_attr(ikey, "bounds", True,
                                                                       allowed_types=str,
                                                                       allowed_values=ikey + "_bounds"))



        if len(data_content_var_names) == 0:
            result.add(ResultCode.ERROR, "No data variable found.")
        elif len(data_content_var_names) == 1:
            if result["data_content"]:
                if self.attrs["data_content"] != data_content_var_names[0]:
                    result["data_content"].add(ResultCode.ERROR, "Only one data variable found. '"+
                                               data_content_var_names[0] + "'. Expected global attribute 'data_content'" +
                                               " to be '"+data_content_var_names[0]+"'.")




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
        if not xy in ["x", "y", "lon", "lat", "E_UTM", "N_UTM",
                      "xs", "ys", "lons", "lats", "Es_UTM", "Ns_UTM",
                      "xu", "yv", "zw", "Eu_UTM", "Ev_UTM", "Nu_UTM", "Nv_UTM",
                      "lonu", "lonv", "latu", "latv"]:
            raise Exception('Unexpected variable: ' + xy)

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
        else:
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

        return out

    def check_var(self, varname, must_exist, allowed_types=None, allowed_range=None, dims=None,
                  must_be_sorted_along=None, decrease_sort_allowed=True, fill_allowed=True):
        exists = varname in self.variables.keys()
        result = CheckResult(ResultCode.OK)

        if not exists:
            if must_exist:
                return CheckResult(ResultCode.ERROR, "Required variable '" + varname + "' not found.")
            else:
                return result

        if allowed_types:
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
                result.add(ResultCode.ERROR, "Variable '" + varname + "' has wrong type. Should be " +
                           "one of the following: " + str(allowed_types))

        if allowed_range:
            if (self[varname].min() < allowed_range[0]) or (self[varname].max() > allowed_range[1]):
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "' is outside allowed range" + str(allowed_range))

        if dims:
            if type(dims) == list:
                dims = tuple(dims)
            elif type(dims) == str:
                dims = tuple([dims])
            if self[varname].dims != dims:
                result.add(ResultCode.ERROR, "Variable '" + varname + "' has wrong dimensions. Expected: " + str(dims))

        if must_be_sorted_along:
            if must_be_sorted_along in self[varname].dims:
                sorted_arr = numpy.sort(self[varname], axis=self[varname].dims.index(must_be_sorted_along))
                if decrease_sort_allowed:
                    anti_sorted_arr = -numpy.sort(-self[varname], axis=self[varname].dims.index(must_be_sorted_along))
                    if not (numpy.array_equal(self[varname], sorted_arr) or numpy.array_equal(self[varname],
                                                                                              anti_sorted_arr)):
                        result.add(ResultCode.ERROR,
                                   "Variable '" + varname + "' must be sorted along dimension '" + must_be_sorted_along + "'")
                else:
                    if not numpy.array_equal(self[varname], sorted_arr):
                        result.add(ResultCode.ERROR,
                                   "Variable '" + varname + "' must be sorted along dimension '" + must_be_sorted_along + "'")

        if fill_allowed:
            pass  # TODO: do it.

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

        if allowed_types:
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

            if not type(self[varname].attrs[attrname]) in allowed_types:
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "': Required variable attribute '" + attrname + "' has wrong type. Should be " +
                           "one of the following: " + str(allowed_types))

        if allowed_values:
            if type(allowed_values) != list:
                allowed_values = [allowed_values]
            if not self[varname].attrs[attrname] in allowed_values:
                if len(allowed_values) == 1:
                    result.add(ResultCode.ERROR,
                               "Variable '" + varname + "': Required variable attribute '" + attrname + "'  has wrong value. Should be " +
                               str(allowed_values[0]))
                else:
                    result.add(ResultCode.ERROR,
                               "Variable '" + varname + "': Required variable attribute '" + attrname + "' has wrong value")

        if allowed_range:
            if self[varname].attrs[attrname] < allowed_range[0] or \
                    self[varname].attrs[attrname] > allowed_range[1]:
                result.add(ResultCode.ERROR,
                           "Variable '" + varname + "': Attribute '" + attrname + "' outside range. Expected: " +
                           str(allowed_range))

        if regex:
            if re.fullmatch(regex, self[varname].attrs[attrname]) is None:
                result.add(ResultCode.ERROR,
                           "Global attribute '" + attrname + "' does not match regular expression " + regex)
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

        if allowed_types:
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

            if not type(self.attrs[attrname]) in allowed_types:
                result.add(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong type. Should be " +
                           "one of the following: " + str(allowed_types))

        if allowed_values:
            if not self.attrs[attrname] in allowed_values:
                if len(allowed_values) == 1:
                    result.add(ResultCode.ERROR,
                               "Global attribute '" + attrname + "' has wrong value. Should be " +
                               str(allowed_values[0]))
                else:
                    result.add(ResultCode.ERROR, "Global attribute '" + attrname + "' has wrong value")

        if regex:
            if re.fullmatch(regex, self.attrs[attrname]) is None:
                result.add(ResultCode.ERROR,
                           "Global attribute '" + attrname + "' does not match regular expression " + regex)

        if max_strlen:
            if len(self.attrs[attrname]) > max_strlen:
                result.add(ResultCode.ERROR,
                           "Global attribute '" + attrname + "' is too long. Must be max. " +
                           str(max_strlen) + " characters.")

        if allowed_range:
            if (self.attrs[attrname] < allowed_range[0]) or (self.attrs[attrname] > allowed_range[1]):
                result.add(ResultCode.ERROR,
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
            ok = True
        else:
            ok = all(i_ok for i_ok in self.result)

        if not ok:
            return ok

        for val in self.values():
            if not val:
                return False

        return True

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
