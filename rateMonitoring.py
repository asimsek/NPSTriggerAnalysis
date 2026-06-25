import warnings, re
warnings.filterwarnings("ignore", message=r".*smallest subnormal.*", category=UserWarning)

from bin.utils import *
from bin.plotting import *
from bin.l1Seed import *
from bin.getEraData import *

def flatten_trigger_dict(trigger_dict):
    trigger_names = []
    for triggers in trigger_dict.values():
        trigger_names.extend(triggers)
    return list(dict.fromkeys(trigger_names))

def build_rate_cache(trigger_dict, t, delivered_lumi, selection_mask=None, target_inst_lumi=2.0e34):
    tree_keys = set(t.keys())
    delivered_lumi = np.asarray(delivered_lumi, dtype=float)
    if selection_mask is not None:
        delivered_lumi = delivered_lumi[selection_mask]
    rate_cache = {}
    missing_triggers = []
    trigger_names = flatten_trigger_dict(trigger_dict)

    for trigger_name in trigger_names:
        possible_branches = (f"{trigger_name}_v", trigger_name)
        branch_name = next((b for b in possible_branches if b in tree_keys), None)

        if branch_name is None:
            missing_triggers.append(trigger_name)
            continue

        counts = ak.to_numpy(t[branch_name].array()).astype(float)
        if selection_mask is not None:
            counts = counts[selection_mask]
        rate = np.full_like(counts, np.nan, dtype=float)
        np.divide(counts, delivered_lumi, out=rate, where=delivered_lumi > 0)
        rate_cache[trigger_name] = rate * target_inst_lumi / 1e36

    for trigger_name in missing_triggers:
        print(f" -- Attention {trigger_name} not in tree!")

    return rate_cache

def run_rate_monitoring(args):
    data_path = args.data
    eras_path = args.eras
    outDir = args.outDir
    trigger_dict_path = args.trigger_dict

    t = args.tree

    # Open the trigger dictionary
    with open(trigger_dict_path) as json_file:
        trigger_dict = json.load(json_file)

    # Open the eras path
    #with open(eras_path) as json_file:
    #    eras = json.load(json_file)
    #eras_list = list(eras.keys())

    # Parse eras file into dictionary format: {2024: ['Ev1', 'G', 'H', 'Iv1'], ...}
    input_eras = {}
    with open(eras_path, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip empty or fully commented lines
            if not line or line.startswith("#"):
                continue

            if ":" not in line:
                continue  # skip malformed

            year_str, eras_str = line.split(":", 1)
            year = int(year_str.strip())

            # Remove inline comment (like # comment at end)
            eras_cleaned = eras_str.split("#")[0] if "#" in eras_str and ":" not in eras_str.split("#")[0] else eras_str

            # Extract individual eras and ignore ones prefixed by '#'
            eras = [
                e.strip()
                for e in eras_str.split(",")
                if e.strip() and not e.strip().startswith("#")
            ]

            if eras:
                input_eras[year] = eras

    # Retrieve era run ranges dynamically
    eras = get_run_era_ranges_dict(input_eras)
    eras_list = list(eras.keys())

    print(eras_list)

    metadata_branches = [
        "run",
        "year",
        "month",
        "day",
        "pileup",
        "physics_flag",
        "cms_ready",
        "recorded_lumi_per_lumisection",
        "beams_stable",
        "delivered_lumi_per_lumisection",
    ]
    metadata = {}
    print("\033[91m[INFO] Reading metadata branches ...\033[0m")
    for branch_name in metadata_branches:
        metadata[branch_name] = ak.to_numpy(t[branch_name].array())

    runs = metadata["run"].astype(int)
    year = metadata["year"]
    month = metadata["month"]
    day = metadata["day"]
    pu = metadata["pileup"]
    golden = metadata["physics_flag"]
    cms_ready = metadata["cms_ready"]
    recorded_lumi = metadata["recorded_lumi_per_lumisection"]
    beams_stable = metadata["beams_stable"]
    delivered_lumi = metadata["delivered_lumi_per_lumisection"]
    integrated_lumi = np.cumsum(recorded_lumi) / 1000.0
    dates = pd.to_datetime({'year': year, 'month': month, 'day': day})
    df = pd.DataFrame({
        'run':  runs,
        'date': dates
    })

    df_ls = pd.DataFrame({
        "run": runs,
        "lumi_ls": recorded_lumi,
    })

    int_lumi_per_run_pb   = df_ls.groupby("run")["lumi_ls"].sum()
    int_lumi_per_run_fb   = int_lumi_per_run_pb / 1000.0

    ### Source for LS_seconds = 23.31s : https://cds.cern.ch/record/2890105/files/DP2024_012.pdf?#page=5
    ### 1e36 means: pb-1 --> cm-2 conversation
    avg_inst_lumi_pb_per_s = df_ls.groupby("run")["lumi_ls"].mean() / LS_seconds 
    avg_inst_lumi_cm2_s    = avg_inst_lumi_pb_per_s * 1e36

    targ_inst_lumi = 2.0e34   # cm^-2 s^-1
    rate_correction_factor = targ_inst_lumi / avg_inst_lumi_cm2_s


    #run_query = 393276
    #if run_query in int_lumi_per_run_pb.index:
    #    print(f"Run {run_query}:")
    #    print(f"  Integrated lumi = {int_lumi_per_run_pb.loc[run_query]:.1f} pb^{-1} "
    #          f"({int_lumi_per_run_fb.loc[run_query]:.3f} fb^{-1})")
    #    print(f"  Avg inst lumi   = {avg_inst_lumi_cm2_s.loc[run_query]:.2e} cm^{-2} s{-1}")
    #else:
    #    print(f"Run {run_query} not in the dataset.")

 
    df_by_run = df.sort_values("run")

    era_dates = {}
    table     = PrettyTable()
    table.field_names = ["Era", "Start Run", "End Run", "Start Date", "End Date"]

    for key, (start_run_s, end_run_s) in eras.items():
        start_run = int(start_run_s)
        end_run   = int(end_run_s)

        start_run = max(start_run, df_by_run['run'].iloc[0])
        end_run   = min(end_run,   df_by_run['run'].iloc[-1])

        sub_start     = df_by_run[df_by_run['run'] >= start_run]
        era_start_date = sub_start['date'].iloc[0] if not sub_start.empty else None

        sub_end       = df_by_run[df_by_run['run'] <= end_run]
        era_end_date   = sub_end['date'].iloc[-1] if not sub_end.empty else None

        # swap if reversed
        if era_start_date and era_end_date and era_start_date > era_end_date:
            era_start_date, era_end_date = era_end_date, era_start_date

        era_df = (
            df_by_run[
                (df_by_run['run'] >= start_run) &
                (df_by_run['run'] <= end_run)   #&
                #(df_by_run['date'].dt.year == 2025)
            ]
            .sort_values("date")
        )
        if era_df.empty:
            continue

        era_df = era_df.copy()
        era_df['date_only'] = era_df['date'].dt.date
        unique_days_df      = era_df.drop_duplicates(subset="date_only")

        first = unique_days_df.iloc[0]
        last  = unique_days_df.iloc[-1]

        era_dates[key] = [era_start_date, era_end_date]

        table.add_row([
            key,
            start_run,
            end_run,
            era_start_date.date() if era_start_date else "N/A",
            era_end_date.date()   if era_end_date   else "N/A"
        ])

    if getattr(args, "print_table", True):
        print(table)

    missing_eras = [key for key in eras_list if key not in era_dates]
    if missing_eras:
        print(
            "\033[91m[WARNING] No lumisections found in the loaded ROOT files "
            "for these requested eras:\033[0m " + ", ".join(missing_eras)
        )

    active_eras_list = [key for key in eras_list if key in era_dates]
    if not active_eras_list:
        raise RuntimeError("No requested eras overlap with the loaded ROOT files.")

    mask_pu = pu >= 60
    mask_golden = golden == 1
    mask_cms_ready = cms_ready == 1
    mask_beams_stable = beams_stable == 1
    mask_runs_in_eras = (runs >= int(eras[active_eras_list[0]][0])) & (runs <= int(eras[active_eras_list[-1]][1]))
    mask_delivered_lumi = delivered_lumi > 0.1

    mask = mask_pu & mask_golden & mask_cms_ready & mask_beams_stable & mask_runs_in_eras & mask_delivered_lumi
    runs = runs[mask]
    pu = pu[mask]
    dates = dates[mask]

    print("\033[91m[INFO] Building trigger rate cache ...\033[0m")
    rate_cache = build_rate_cache(trigger_dict, t, delivered_lumi, selection_mask=mask)

    figs_run = []
    figs_date = []
    for group in list(trigger_dict.keys()):
        print (group)
        figs_run.append(plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis='run', x_label="Run Number", runs=runs, dates=dates, pu=pu, eras_list=active_eras_list, print_trig=True, rate_cache=rate_cache, rate_cache_is_masked=True))
        figs_date.append(plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis='date', x_label="Date", runs=runs, dates=dates, pu=pu, eras_list=active_eras_list, print_trig=False, rate_cache=rate_cache, rate_cache_is_masked=True))

    multipage(outDir + "/NPSTriggerMonitoring_run_AllCombined.pdf", figs=figs_run, dpi=50)
    multipage(outDir + "/NPSTriggerMonitoring_date_AllCombined.pdf", figs=figs_date, dpi=50)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CMS NPS Trigger / L1-Seed rate monitoring utility")

    # Positional arguments (existing)
    parser.add_argument("trigger_dict", help="Path to triggerNames_*.json (HLT paths)")
    parser.add_argument("eras", help="Path to eras JSON")
    parser.add_argument("outDir", help="Output directory for plots")
    parser.add_argument("data", help="ROOT file containing rate tree")

    # Optional flags
    parser.add_argument(
        "--l1seed",
        action="store_true",
        help="Generate L1-seed trigger dictionary and run additional monitoring",
    )
    parser.add_argument(
        "--gRun",
        action="store_true",
        help="Regenerate GRun.csv (requires --l1seed) before extracting seeds",
    )

    parser.add_argument(
        "--noHLT",
        action="store_true",
        help="No HLT monitoring",
    )

    args = parser.parse_args()

    tree = open_tree_any(args.data)
    args.tree = tree
    args.print_table = True

    # Ensure base output directory exists
    Path(args.outDir).mkdir(parents=True, exist_ok=True)

    # Run HLT rate monitoring if "--noHLT" is not used
    if not args.noHLT:
        run_rate_monitoring(args)

    # Optional L1-seed flow
    if args.l1seed:
        # If requested, (re)generate the menu CSV for seed extraction
        gRun_csv = Path("jsonFiles/L1SeedLists/GRun.csv")
        if args.gRun:
            cmd = (
                "hltGetConfiguration /dev/CMSSW_15_1_0/GRun | "
                "hltDumpStream --mode csv > jsonFiles/L1SeedLists/GRun.csv"
            )
            print("\033[91m[INFO] Regenerating GRun.csv via hltGetConfiguration ...\033[0m")
            subprocess.run(cmd, shell=True, check=True)
        elif not gRun_csv.exists():
            print(
                "[WARNING] GRun.csv not found and --gRun not specified. "
                "Seed extraction will fail since the file is required."
            )

        # Extract seeds and get new dictionary path
        print("\033[91m[INFO] Extracting L1 seed list ...\033[0m")
        new_dict_path = extract_l1_seeds(args.trigger_dict, str(gRun_csv))

        # Build dedicated output directory
        l1_outdir = Path(args.outDir) / "L1SeedRates"
        l1_outdir.mkdir(parents=True, exist_ok=True)

        # Clone args to reuse run_rate_monitoring
        l1_args = argparse.Namespace(**vars(args))
        l1_args.trigger_dict = str(new_dict_path)
        l1_args.eras = str (args.eras)
        l1_args.outDir = str(l1_outdir)
        l1_args.data = str (args.data)

        l1_args.tree = tree
        l1_args.print_table = False

        print(f"\033[91m[INFO] Running L1-seed rate monitoring ->\033[0m {l1_outdir}")
        run_rate_monitoring(l1_args)
    elif args.gRun:
        # --gRun without --l1seed: warning
        print("[WARNING] --gRun specified without --l1seed. Ignoring --gRun.")




