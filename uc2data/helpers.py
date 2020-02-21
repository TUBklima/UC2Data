from .Dataset import Dataset
from pathlib import Path


def check_multi(folder):
    pathlist = Path(folder).glob('**/*.nc')

    for path in pathlist:
        outfile = Path(str(path).replace(".nc", ".check"))

        if not outfile.exists():

            try:
                i_data = Dataset(path)
                i_data.uc2_check()
                print_me = str(i_data.check_result.errors) + str(i_data.check_result.warnings)
            except Exception:
                print_me = "Could not read file: "+str(path)

            text_file = open(str(outfile), "w")
            text_file.write(print_me)
            text_file.close()
