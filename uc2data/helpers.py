from .Dataset import Dataset
from pathlib import Path


def check_multi(folder):
    pathlist = Path(folder).glob('**/*.nc')

    for path in pathlist:
        outfile = Path(str(path).replace(".nc", ".check"))

        try:
            with Dataset(path) as i_data:
                i_data.uc2_check()
                print_me_err = str(i_data.check_result.errors)
                print_me_warn = str(i_data.check_result.warnings)
                if print_me_err != "" and print_me_warn != "":
                    connect = "\n"
                else:
                    connect = ""
                print_me = print_me_err + connect + print_me_warn
        except Exception:
            print_me = "Could not read file: "+str(path)

        text_file = open(str(outfile), "w")
        text_file.write(print_me)
        text_file.close()
