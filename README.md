# UC2Data
Handles NetCDF files in the [UC]² data standard. Includes checks for conformity.

## Installation

The following line should do the job

`pip install git+https://github.com/TUBklima/UC2Data.git`

## Open a file

`import uc2data  # import the package`

`my_dataset = uc2data.Dataset(filename)`

## Access data in the file

The data itself is represented as an xarray dataset. Access it like this:

`data = my_dataset.ds`

## Check for conformity with [UC]² data standard

`my_dataset.uc2_check()  # do the check`

`print(my_dataset.check_result)  # the results are saved within the dataset's attribute check_result`

`print(my_dataset.check_result.errors)  # get just the errors`

`print(my_dataset.check_result.warnings)  # get just the warnings`
