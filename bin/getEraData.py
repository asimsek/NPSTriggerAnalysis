import re
import subprocess
from collections import OrderedDict

def get_run_era_ranges_dict(years_and_eras):
    """
    For each year+era in `years_and_eras`, list
      /ExpressPhysics/Run{year}{era}-Express-v*/FEVT
    via dasgoclient, then for each matching dataset:
      - if only one vX existskey = "{year}{era}"
      - if multiple key = "{year}{era}v{X}"
    finally query all runs in each dataset and return min/max.
    """
    result = {}

    for year, era_list in years_and_eras.items():
        for era in era_list:
            base_pattern = f"/ExpressPhysics/Run{year}{era}-Express-v*/FEVT"
            # find all versions
            ds_proc = subprocess.run(
                ["dasgoclient", "-query", f"dataset={base_pattern}"],
                capture_output=True, text=True, check=False
            )
            ds_names = [l.strip() for l in ds_proc.stdout.splitlines() if l.strip()]
            if not ds_names:
                raise RuntimeError(f"No ExpressPhysics datasets found for era {year}{era}")

            # decide whether to split by version
            split_by_version = len(ds_names) > 1
            for ds in ds_names:
                # extract version number
                m = re.search(r"-v(\d+)/FEVT$", ds)
                ver = m.group(1) if m else "1"

                # build the key
                if split_by_version:
                    key = f"{year}{era}v{ver}"
                else:
                    key = f"{year}{era}"

                # query runs for this dataset
                run_proc = subprocess.run(
                    ["dasgoclient", "-query", f"run dataset={ds}"],
                    capture_output=True, text=True, check=False
                )
                runs = [int(x) for x in run_proc.stdout.split() if x.isdigit()]
                if not runs:
                    raise RuntimeError(f"No runs found in dataset {ds}")
                result[key] = [str(min(runs)), str(max(runs))]

    # sort by year then by min_run
    return OrderedDict(
        sorted(
            result.items(),
            key=lambda it: (int(it[0][:4]), int(it[1][0]))
        )
    )

