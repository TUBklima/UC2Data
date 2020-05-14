from .Dataset import Dataset
from pathlib import Path


def check_multi(folder):
    pathlist = Path(folder).glob('**/*.nc')

    for path in pathlist:
        outfile = Path(str(path).replace(".nc", ".check"))

        try:
            with Dataset(path) as i_data:
                i_data.uc2_check()
                i_data.check_result.to_file(outfile, full=False)
        except Exception:
            text_file = open(str(outfile), "w")
            text_file.write("Could not read file: "+str(path))
            text_file.close()

