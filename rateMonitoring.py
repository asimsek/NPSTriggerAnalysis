from bin.utils import *
from bin.plotting import *

def run_rate_monitoring(args):
    data_path = args.data
    eras_path = args.eras
    outDir = args.outDir
    trigger_dict_path = args.trigger_dict

    # Open the data file
    f = uproot.open(data_path)
    t = f['tree']

    # Open the trigger dictionary
    with open(trigger_dict_path) as json_file:
        trigger_dict = json.load(json_file)

    # Open the eras path
    with open(eras_path) as json_file:
        eras = json.load(json_file)
    eras_list = list(eras.keys())

    runs = ak.to_numpy(t["run"].array()).astype(int)
    year = ak.to_numpy(t["year"].array())
    month = ak.to_numpy(t["month"].array())
    day = ak.to_numpy(t["day"].array())
    pu = ak.to_numpy(t["pileup"].array())
    golden = ak.to_numpy(t["physics_flag"].array())
    cms_ready = ak.to_numpy(t["cms_ready"].array())
    recorded_lumi = ak.to_numpy(t["recorded_lumi_per_lumisection"].array())
    beams_stable = ak.to_numpy(t["beams_stable"].array())
    delivered_lumi = ak.to_numpy(t["delivered_lumi_per_lumisection"].array())
    integrated_lumi = np.cumsum(recorded_lumi) / 1000.0
    dates = pd.to_datetime({'year': year, 'month': month, 'day': day})

    # Loop to find dates for the given eras - for vertical era separation lines
    era_dates = {}
    table = PrettyTable()
    table.field_names = ["Era", "Start Run", "End Run", "Start Date", "End Date"]
    for key, (start_run, end_run) in eras.items():
        start_run, end_run = int(start_run), int(end_run)  # Ensure both start and end run are integers
        # Find start date
        start_date_idx = np.where(runs >= start_run)[0]
        start_date = dates[start_date_idx[0]] if len(start_date_idx) > 0 else None
        # Find end date
        end_date_idx = np.where(runs <= end_run)[0]
        end_date = dates[end_date_idx[-1]] if len(end_date_idx) > 0 else None
        era_dates[key] = [start_date, end_date]
        table.add_row([
            key,
            start_run,
            end_run,
            start_date.date() if start_date else "N/A",
            end_date.date() if end_date else "N/A"
        ])

    print(table)

    mask_pu = pu >= 62
    mask_golden = golden == 1
    mask_cms_ready = cms_ready == 1
    mask_beams_stable = beams_stable == 1
    mask_runs_in_eras = (runs >= int(eras[eras_list[0]][0])) & (runs <= int(eras[eras_list[-1]][1]))
    mask_delivered_lumi = delivered_lumi > 0.1

    mask = mask_pu & mask_golden & mask_cms_ready & mask_beams_stable & mask_runs_in_eras & mask_delivered_lumi
    runs = runs[mask]
    pu = pu[mask]
    dates = dates[mask]

    figs_run = []
    figs_date = []
    for group in list(trigger_dict.keys()):
        print (group)
        figs_run.append(plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis='run', x_label="Run Number", runs=runs, dates=dates, pu=pu, eras_list=eras_list, print_trig=True))
        figs_date.append(plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis='date', x_label="Date", runs=runs, dates=dates, pu=pu, eras_list=eras_list, print_trig=False))

    multipage(outDir + "/SUSTriggerMonitoring_run_AllCombined.pdf", figs=figs_run, dpi=50)
    multipage(outDir + "/SUSTriggerMonitoring_date_AllCombined.pdf", figs=figs_date, dpi=50)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('trigger_dict')
    parser.add_argument('eras')
    parser.add_argument('outDir')
    parser.add_argument('data')

    args = parser.parse_args()
    print(args.trigger_dict, args.eras, args.data)

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    run_rate_monitoring(args)





