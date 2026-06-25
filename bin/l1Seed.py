import json
import re
from pathlib import Path
from typing import List, Dict

import pandas as pd

def _flatten_trigger_dict(trigger_dict: Dict[str, List[str]]) -> List[str]:
    """Flatten nested dictionary values into a single list of HLT paths."""
    return [p for sub in trigger_dict.values() for p in sub]

def _hlt_path_pattern(hlt_path: str) -> str:
    """Match exactly HLT_Path or its CMSSW versioned HLT_Path_vN form."""
    hlt_path = hlt_path.strip()
    if re.search(r"_v\d+$", hlt_path):
        return rf"^{re.escape(hlt_path)}$"
    return rf"^{re.escape(hlt_path)}(?:_v\d+)?$"

def extract_l1_seeds(trigger_list_json_path: str, csv_path: str) -> Path:
    """Return a JSON file that maps every HLT path → list of L1 seeds."""
    trigger_list_path = Path(trigger_list_json_path).expanduser().resolve()
    csv_path          = Path(csv_path).expanduser().resolve()

    with open(trigger_list_path) as f:
        trigger_dict = json.load(f)

    hlt_paths = _flatten_trigger_dict(trigger_dict)

    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()

    if "path" not in df.columns or "seed" not in df.columns:
        raise ValueError("CSV must contain columns 'path' and 'seed'.")

    # pre-strip every CSV path once up-front
    df["path_clean"] = df["path"].astype(str).str.strip()

    result = {}
    for raw_path in hlt_paths:
        path_stripped = raw_path.strip()

        # Require an exact path match, allowing only the CMSSW _vN suffix.
        # This avoids matching e.g. HLT_Foo and HLT_Foo_DZ as the same path.
        sel = df["path_clean"].str.match(_hlt_path_pattern(path_stripped), na=False)

        if sel.any():
            seed_expr = df.loc[sel, "seed"].iloc[0]
            if isinstance(seed_expr, str) and seed_expr.lower() != "(none)":
                seeds = [s.strip() for s in seed_expr.split("OR") if s.strip()]
                result[raw_path] = seeds
            else:
                result[raw_path] = []
        else:
            result[raw_path] = []

    year = trigger_list_path.parents[1].name if len(trigger_list_path.parents) > 1 else "unknown"
    tag  = trigger_list_path.stem.replace("triggerNames", "seedNames")

    out_dir  = trigger_list_path.parents[2] / "L1SeedLists" / year
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{tag}.json"

    with open(out_file, "w") as f:
        json.dump(result, f, indent=4, sort_keys=True)

    return out_file


