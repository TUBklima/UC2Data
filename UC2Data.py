import xarray
import sys
from warnings import warn

is_win = sys.platform in ['win32', 'win64']
if not is_win:
    from cfchecker import cfchecks


class UC2Data(xarray.Dataset):

    def __init__(self, path):
        self.path = path
        super().__init__()

    def read(self):
        tmp = xarray.open_dataset(self.path)
        self.update(tmp, inplace=True)

    def check(self):
        self.cf_check()
        self.uc2_check()

    def cf_check(self):
        if is_win:
            warn('CFchecker does not work on Windows systems. Sorry')

    def uc2_check(self):
        pass
