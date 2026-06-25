import warnings, re, time
warnings.filterwarnings("ignore", message=r".*smallest subnormal.*", category=UserWarning)

from bin.utils import *
from bin.plotting import *
from bin.l1Seed import *
from bin.getEraData import *
from bin.rateChangeMatrix import *
from bin.summaryHeatmaps import *

def format_runtime(seconds):
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.1f} s"
    minutes, rem_seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)} min {rem_seconds:.1f} s"
    hours, rem_minutes = divmod(minutes, 60)
    return f"{int(hours)} h {int(rem_minutes)} min {rem_seconds:.1f} s"

def flatten_trigger_dict(trigger_dict):
    trigger_names = []
    for triggers in trigger_dict.values():
        trigger_names.extend(triggers)
    return list(dict.fromkeys(trigger_names))

def read_tree_arrays(t, branch_names):
    branch_names = list(dict.fromkeys(branch_names))
    if not branch_names:
        return {}
    arrays = t.arrays(branch_names, library="ak")
    return {branch_name: arrays[branch_name] for branch_name in branch_names}

def build_rate_cache(trigger_dict, t, delivered_lumi, selection_mask=None, target_inst_lumi=2.0e34):
    delivered_lumi = np.asarray(delivered_lumi, dtype=float)
    if selection_mask is not None:
        delivered_lumi = delivered_lumi[selection_mask]
    rate_cache = {}
    trigger_names = flatten_trigger_dict(trigger_dict)

    if hasattr(t, "arrays_any"):
        candidates_by_trigger = {
            trigger_name: (f"{trigger_name}_v", trigger_name)
            for trigger_name in trigger_names
        }
        count_arrays, missing_triggers = t.arrays_any(candidates_by_trigger, library="ak")
    else:
        tree_keys = set(t.keys())
        branch_to_trigger = {}
        missing_triggers = []
        for trigger_name in trigger_names:
            possible_branches = (f"{trigger_name}_v", trigger_name)
            branch_name = next((b for b in possible_branches if b in tree_keys), None)

            if branch_name is None:
                missing_triggers.append(trigger_name)
                continue

            branch_to_trigger[branch_name] = trigger_name

        branch_arrays = read_tree_arrays(t, branch_to_trigger.keys())
        count_arrays = {
            trigger_name: branch_arrays[branch_name]
            for branch_name, trigger_name in branch_to_trigger.items()
        }

    for trigger_name, count_array in count_arrays.items():
        if trigger_name in missing_triggers:
            continue
        counts = ak.to_numpy(count_array).astype(float)
        if selection_mask is not None:
            counts = counts[selection_mask]
        rate = np.full_like(counts, np.nan, dtype=float)
        np.divide(counts, delivered_lumi, out=rate, where=delivered_lumi > 0)
        rate_cache[trigger_name] = rate * target_inst_lumi / 1e36

    for trigger_name in missing_triggers:
        print(f" -- Attention {trigger_name} not in tree!")

    return rate_cache

def parse_eras_file(eras_path):
    input_eras = {}
    with open(eras_path, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if ":" not in line:
                continue

            year_str, eras_str = line.split(":", 1)
            year = int(year_str.strip())

            eras = [
                e.strip()
                for e in eras_str.split(",")
                if e.strip() and not e.strip().startswith("#")
            ]

            if eras:
                input_eras[year] = eras
    return input_eras

def prepare_monitoring_context(args):
    t = args.tree
    input_eras = parse_eras_file(args.eras)

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
    print("\033[91m[INFO] Reading metadata branches ...\033[0m")
    metadata_arrays = read_tree_arrays(t, metadata_branches)
    metadata = {
        branch_name: ak.to_numpy(metadata_arrays[branch_name])
        for branch_name in metadata_branches
    }

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
    dates = pd.to_datetime({'year': year, 'month': month, 'day': day})
    df = pd.DataFrame({
        'run':  runs,
        'date': dates
    })

    df_by_run = df.sort_values("run")

    era_dates = {}
    table = PrettyTable()
    table.field_names = ["Era", "Start Run", "End Run", "Start Date", "End Date"]

    for key, (start_run_s, end_run_s) in eras.items():
        start_run = int(start_run_s)
        end_run = int(end_run_s)

        start_run = max(start_run, df_by_run['run'].iloc[0])
        end_run = min(end_run, df_by_run['run'].iloc[-1])

        sub_start = df_by_run[df_by_run['run'] >= start_run]
        era_start_date = sub_start['date'].iloc[0] if not sub_start.empty else None

        sub_end = df_by_run[df_by_run['run'] <= end_run]
        era_end_date = sub_end['date'].iloc[-1] if not sub_end.empty else None

        if era_start_date and era_end_date and era_start_date > era_end_date:
            era_start_date, era_end_date = era_end_date, era_start_date

        era_df = (
            df_by_run[
                (df_by_run['run'] >= start_run) &
                (df_by_run['run'] <= end_run)
            ]
            .sort_values("date")
        )
        if era_df.empty:
            continue

        era_dates[key] = [era_start_date, era_end_date]

        table.add_row([
            key,
            start_run,
            end_run,
            era_start_date.date() if era_start_date else "N/A",
            era_end_date.date() if era_end_date else "N/A"
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

    return {
        "tree": t,
        "eras": eras,
        "era_dates": era_dates,
        "active_eras_list": active_eras_list,
        "mask": mask,
        "runs": runs[mask],
        "pu": pu[mask],
        "dates": dates[mask],
        "delivered_lumi": delivered_lumi,
    }

def run_rate_monitoring(args, monitoring_context=None):
    data_path = args.data
    eras_path = args.eras
    outDir = args.outDir
    trigger_dict_path = args.trigger_dict

    if monitoring_context is None:
        monitoring_context = prepare_monitoring_context(args)

    # Open the trigger dictionary
    with open(trigger_dict_path) as json_file:
        trigger_dict = json.load(json_file)

    t = monitoring_context["tree"]
    eras = monitoring_context["eras"]
    era_dates = monitoring_context["era_dates"]
    active_eras_list = monitoring_context["active_eras_list"]
    mask = monitoring_context["mask"]
    runs = monitoring_context["runs"]
    pu = monitoring_context["pu"]
    dates = monitoring_context["dates"]
    delivered_lumi = monitoring_context["delivered_lumi"]

    print("\033[91m[INFO] Building trigger rate cache ...\033[0m")
    rate_cache = build_rate_cache(trigger_dict, t, delivered_lumi, selection_mask=mask)

    if not getattr(args, "noRateMatrix", False):
        print("\033[91m[INFO] Building adjacent-era rate-change matrix ...\033[0m")
        plot_rate_change_matrix(
            trigger_dict,
            rate_cache,
            runs,
            eras,
            active_eras_list,
            outDir,
            merge_era_versions=not getattr(args, "splitEraVersions", False),
        )

    if not getattr(args, "noSummaryHeatmaps", False):
        print("\033[91m[INFO] Building trigger stability heatmap ...\033[0m")
        plot_trigger_stability_heatmap(
            trigger_dict,
            rate_cache,
            runs,
            eras,
            active_eras_list,
            outDir,
            merge_era_versions=not getattr(args, "splitEraVersions", False),
        )
        print("\033[91m[INFO] Building year-normalized era heatmap ...\033[0m")
        plot_year_normalized_era_heatmap(
            trigger_dict,
            rate_cache,
            runs,
            eras,
            active_eras_list,
            outDir,
            merge_era_versions=not getattr(args, "splitEraVersions", False),
        )

    if getattr(args, "noTrendPlots", False):
        return trigger_dict, rate_cache

    plot_cache = {}
    figs_run = []
    figs_date = []
    make_run_plots = not getattr(args, "noRunPlots", False)
    make_date_plots = not getattr(args, "noDatePlots", False)

    for group in list(trigger_dict.keys()):
        print(group)
        if make_run_plots:
            figs_run.append(plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis='run', x_label="Run Number", runs=runs, dates=dates, pu=pu, eras_list=active_eras_list, print_trig=True, rate_cache=rate_cache, rate_cache_is_masked=True, plot_cache=plot_cache))
        if make_date_plots:
            figs_date.append(plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis='date', x_label="Date", runs=runs, dates=dates, pu=pu, eras_list=active_eras_list, print_trig=not make_run_plots, rate_cache=rate_cache, rate_cache_is_masked=True, plot_cache=plot_cache))

    if figs_run:
        multipage(outDir + "/NPSTriggerMonitoring_run_AllCombined.pdf", figs=figs_run, dpi=50)
    if figs_date:
        multipage(outDir + "/NPSTriggerMonitoring_date_AllCombined.pdf", figs=figs_date, dpi=50)

    return trigger_dict, rate_cache

if __name__ == "__main__":
    total_start_time = time.perf_counter()
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
    parser.add_argument(
        "--noRateMatrix",
        action="store_true",
        help="Skip the adjacent-era rate-change matrix",
    )
    parser.add_argument(
        "--noSummaryHeatmaps",
        action="store_true",
        help="Skip trigger stability, year-normalized, and HLT/dominant-L1 summary heatmaps",
    )
    parser.add_argument(
        "--noTrendPlots",
        action="store_true",
        help="Skip both run and date time-series trend PDFs",
    )
    parser.add_argument(
        "--noRunPlots",
        action="store_true",
        help="Skip run-number time-series trend PDFs",
    )
    parser.add_argument(
        "--noDatePlots",
        action="store_true",
        help="Skip date time-series trend PDFs",
    )
    parser.add_argument(
        "--noL1TrendPlots",
        action="store_true",
        help="With --l1seed, skip only the L1-seed run/date trend PDFs",
    )
    parser.add_argument(
        "--splitEraVersions",
        action="store_true",
        help="Keep v1/v2 era versions separate in the rate-change matrix",
    )

    args = parser.parse_args()

    tree = open_tree_any(args.data)
    args.tree = tree
    args.print_table = True

    # Ensure base output directory exists
    Path(args.outDir).mkdir(parents=True, exist_ok=True)

    monitoring_context = None
    if (not args.noHLT) or args.l1seed:
        monitoring_context = prepare_monitoring_context(args)

    # Run HLT rate monitoring if "--noHLT" is not used
    hlt_result = None
    if not args.noHLT:
        hlt_result = run_rate_monitoring(args, monitoring_context=monitoring_context)

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
        if getattr(args, "noL1TrendPlots", False):
            l1_args.noTrendPlots = True

        print(f"\033[91m[INFO] Running L1-seed rate monitoring ->\033[0m {l1_outdir}")
        l1_result = run_rate_monitoring(l1_args, monitoring_context=monitoring_context)

        if not getattr(args, "noSummaryHeatmaps", False):
            if hlt_result is None:
                print("\033[91m[WARNING] Skipping HLT/dominant-L1 heatmap because HLT monitoring was disabled.\033[0m")
            elif l1_result is None:
                print("\033[91m[WARNING] Skipping HLT/dominant-L1 heatmap because L1-seed monitoring did not return rates.\033[0m")
            else:
                hlt_trigger_dict, hlt_rate_cache = hlt_result
                l1_seed_dict, l1_rate_cache = l1_result
                print("\033[91m[INFO] Building HLT/dominant-L1 seed heatmap ...\033[0m")
                plot_hlt_dominant_l1_seed_heatmap(
                    hlt_trigger_dict,
                    hlt_rate_cache,
                    l1_seed_dict,
                    l1_rate_cache,
                    monitoring_context["runs"],
                    monitoring_context["eras"],
                    monitoring_context["active_eras_list"],
                    args.outDir,
                    merge_era_versions=not getattr(args, "splitEraVersions", False),
                )
    elif args.gRun:
        # --gRun without --l1seed: warning
        print("[WARNING] --gRun specified without --l1seed. Ignoring --gRun.")

    total_runtime = time.perf_counter() - total_start_time
    print(f"\033[91m[INFO] Total run time:\033[0m {format_runtime(total_runtime)}")

