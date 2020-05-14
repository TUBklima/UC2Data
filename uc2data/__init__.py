from .Dataset import Dataset
from .Result import ResultCode, ResultItem, CheckResult
from .helpers import check_multi
from warnings import warn

data_standard_version = (1,4,1)

msg = "\nThe uc2data repository migrated to https://gitlab.klima.tu-berlin.de/klima/uc2data.git \n"+ \
      "You are using an old version from https://github.com/TUBklima/UC2Data.git \n"+\
      "Please checkout the latest version from the new repository location"
warn(msg)
