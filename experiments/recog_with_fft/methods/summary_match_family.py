import pandas as pd
from glob import glob
from sys import argv as args
from tqdm import tqdm
from tools import eprint
import numpy as np
import re

res = pd.DataFrame(columns=[
    "Window_Size",
    "N_Peaks",
    "Overlap",
    "Skip_First_K_Freqs",
    "Selection_Method",
    "Significance",
    "Average_F_Score",
    "Average_Precision",
    "Average_Sharpness",
    "avg(Family_Members/Match_Count)",
    "Average_Hash_Intersection_Size",
    "Runtime_match-family"
    ]
)
for file in tqdm(args[1:]):
    try:
        summary = pd.read_csv(file)
        logfile = file.split("/")
        logfile.insert(-1, "_logs")
        logfile = glob("%s*.err" % "/".join(logfile).rsplit(".", 2)[0])[0]
        params = re.findall('WINSIZE_(\d+)_NPEAKS_(\d+)_OVERLAP_(\d+)_SKIPK_(\d+)_SELMETH_(\w+)_ALPHA_(\d+\.?\d?\d?\d?\d?)', file)[0]
        f_scores = summary["F_Score"][~summary["F_Score"].isna()]

        with open(logfile, "r") as f:
            log = f.readlines()
            match_family_time = re.search("([:\d]+)<[:\d]+,[ \d\.it/s\]]+$", log[-1])
            del log

            def pad(array: np.array) -> np.array:
                return np.pad(array, (0, 4937-len(array)))

            if match_family_time is not None:
                res.loc[len(res.index)] = (
                    *params,
                    pad(f_scores).mean(),
                    pad(summary["Precision"]).mean(),
                    pad(summary["Sharpness"]).mean(),
                    (summary["Member_Count"] / summary["Match_Count"]).mean(),
                    pad(summary["Hash_Intersec_Size"]).mean(),
                    match_family_time.group(1)
                )
    except Exception as e:
        if str(e) != "No columns to parse from file":
            raise e

print(res.sort_values("Average_F_Score", ascending=False).to_csv(index=False, float_format="%g"))
