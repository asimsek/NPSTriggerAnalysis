import warnings, re
warnings.filterwarnings("ignore", message=r".*smallest subnormal.*", category=UserWarning)

from bin.utils import *

matplotlib.use('Agg')

def run_rate_check(args):

    data_path = args.data
    eras_path = args.eras
    outDir = args.outDir
    trigger_dict_path = args.trigger_dict

    # Open the data file 
    t = args.tree

    # Open the trigger dictionary
    with open(trigger_dict_path) as json_file:
        trigger_dict = json.load(json_file)

    # Create ungrouped list of trigger names 
    trigger_names = []
    for key, value in trigger_dict.items():
        trigger_names += value
    trigger_names = sorted(set(trigger_names))

    # Open the eras path 
    with open(eras_path) as json_file:
        eras = json.load(json_file)
    eras_list = list(eras.keys())

    avg_rates = {}
    std_rates = {}
    avg_pu = {}
    plotted_triggers = {}
    max_rates = {}

    # Initialize dictionaries for median difference and ratio plots
    median_diff = {group: [] for group in trigger_dict}
    median_ratio = {group: [] for group in trigger_dict}

    # Check that lumi-leveling looks okay 
    fig, ax = plt.subplots(figsize=(16, 10))
    colors = plt.cm.get_cmap('tab10', len(eras))

    for i, era in enumerate(eras):
        avg_rates[era] = []
        std_rates[era] = []
        avg_pu[era] = []
        plotted_triggers[era] = []
        max_rates[era] = []
        print(era)
        runs = t["run"].array()
        lumi = t["lumi"].array()
        pileup = t["pileup"].array()
        golden = t["physics_flag"].array()
        cms_ready = t["cms_ready"].array()
        pileup_smoothed = median_filter(pileup, size=10)
        mask = (runs == eras[era]["run"] ) & (lumi > eras[era]["ls"][0]) & (lumi < eras[era]["ls"][1]) & (golden == 1) & (cms_ready == 1)
        print(runs[mask])
        for trigger_name in trigger_names:
            if trigger_name + "_v" in t.keys():
                trigger = t[trigger_name + "_v"].array()
                trigger_scaled = trigger / LS_seconds
                trigger_scaled_smoothed = median_filter(trigger_scaled, size=50)
                
                # Updated to a scatter plot for a 2D visualization
                sc = ax.scatter(lumi[mask], pileup_smoothed[mask], color=colors(i), alpha=1.0, s=100, label="{0}[{1}]".format(era, eras[era]["run"]), edgecolor='w', rasterized=True)
                #sc = ax.plot(lumi[mask], trigger_scaled_smoothed[mask], color=colors(i), alpha=1.0, label=era)

                avg_rates[era] += [np.mean(trigger_scaled[mask])]
                std_rates[era] += [np.mean(trigger_scaled[mask])]
                avg_pu[era] += [np.mean(pileup[mask])]
                plotted_triggers[era] += [trigger_name]
                max_rates[era].append(np.max(trigger_scaled[mask]))

    # Sort legend labels by max rate
    handles, labels = ax.get_legend_handles_labels()
    #handles = sorted(handles, key=lambda x: max_rates.get(labels[handles.index(x)], 0), reverse=True)
    #labels = sorted(labels, key=lambda x: max_rates.get(x, 0), reverse=True)

    # Zip together any labels that are duplicated
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)

    lgnd = ax.legend(newHandles, newLabels, frameon=False, bbox_to_anchor=(0.45, -0.12), loc='upper center', fontsize=18, ncol=math.ceil(len(eras)/2), handlelength=1, columnspacing=1.5)
    for lh in range(len(lgnd.legend_handles)):
        lgnd.legend_handles[lh]._sizes = [150]
    ax.set_xlabel("Lumi-section")
    ax.set_ylabel("Pileup")
    plt.subplots_adjust(top=0.95, right=0.97, left=0.07, bottom=0.20)
    plt.savefig(outDir + "/lumiLeveling_pileup.pdf")
    plt.close(fig)

    # Plot differences
    plt.figure()
    era_names = list(eras.keys())
    print("---------------------------")
    for era in range(len(eras) - 1):
        era_1 = era_names[era]
        era_2 = era_names[era + 1]
        print(era_1 + "     -       " + era_2)
        difference = np.array(avg_rates[era_2]) - np.array(avg_rates[era_1])
        mean = np.median(difference[~np.isnan(difference) & ~np.isinf(difference)])
        plt.hist(difference, bins=100, range=(-10, 10),
                 alpha=0.3,
                 label=era_2 + " - " + era_1 + "\n Median: " + "{:.2f}".format(mean))
        mask = (difference > 10) | (difference < -10)
        for i in range(np.sum(mask)):  
            print(np.array(plotted_triggers[era_1])[mask][i] + ": " + "{:.0f}".format(difference[mask][i]) + " Hz")
    plt.xlabel("Difference Rate in Era 2 - Rate in Era 1 [Hz]")
    plt.legend(fontsize=16)
    plt.ylabel("Number of trigger paths")
    plt.subplots_adjust(top=0.95, right=0.95, left=0.12, bottom=0.10)
    plt.savefig(outDir + "/rateAnalysis_differences.pdf")
    print("---------------------------")

    # Make ratio plots of rates
    plt.figure()
    for era in range(len(eras) - 1):
        era_1 = era_names[era]
        era_2 = era_names[era + 1]
        #print(era_1 + "-" + era_2 + "---------------------------")
        ratio = np.array(avg_rates[era_2]) / np.array(avg_rates[era_1])
        mean = np.median(ratio[~np.isnan(ratio) & ~np.isinf(ratio)])
        plt.hist(ratio, bins=100, range=(0, 5),
                 alpha=0.3,
                 label=era_2 + " / " + era_1 + "\n Median: " + "{:.2f}".format(mean))
        mask = (ratio > 2) | (ratio < 0.5)
        for i in range(np.sum(mask)):
            print(np.array(plotted_triggers[era_1])[mask][i])
            print("Ratio:" + "{:.2f}".format(ratio[mask][i]) + " Hz")
            print(era_1 + ": " + "{:.2f}".format(np.array(avg_rates[era_1])[mask][i]) + " Hz")
            print(era_2 + ": " + "{:.2f}".format(np.array(avg_rates[era_2])[mask][i]) + " Hz")
    plt.xlabel("Ratio of Rate in Era 2 / Rate in Era 1")
    plt.legend(fontsize=16)
    plt.ylabel("Number of trigger paths")
    plt.subplots_adjust(top=0.95, right=0.95, left=0.12, bottom=0.10)
    plt.savefig(outDir + "/rateAnalysis_ratios.pdf")

    for group, triggers in trigger_dict.items():
        # Create necessary sub-folders if they don't exist
        median_diff_dir = os.path.join(outDir, "rateAnalysis_medianDiff")
        median_ratio_dir = os.path.join(outDir, "rateAnalysis_medianRatio")
        os.makedirs(median_diff_dir, exist_ok=True)
        os.makedirs(median_ratio_dir, exist_ok=True)

        # Plot differences
        plt.figure()
        markers = ['o', 's', '^', 'D', 'v', '>', '<', 'p', '*', 'h']  # List of markers for different lines
        x_labels = [f"{eras_list[era + 1]} - {eras_list[era]}" for era in range(len(eras) - 1)]
        for idx, trig in enumerate(triggers):
            differences = []
            for era in range(len(eras) - 1):
                era_1 = eras_list[era]
                era_2 = eras_list[era + 1]
                
                if trig in trigger_names:
                    group_avg_rate_era1 = avg_rates[era_1][trigger_names.index(trig)]
                    group_avg_rate_era2 = avg_rates[era_2][trigger_names.index(trig)]
                    difference = group_avg_rate_era2 - group_avg_rate_era1
                    differences.append(difference)
            plt.plot(x_labels, differences, label=trig, marker=markers[idx % len(markers)], markersize=8)
        plt.xlabel("Eras", fontsize=18, fontweight='bold')
        plt.ylabel("Difference in Rate [Hz]", fontsize=18, fontweight='bold')
        plt.legend(prop={'size': 14, 'weight': 'bold'})
        plt.ylim(-5, 5)
        plt.title(f"Median Differences Between Eras - Group: {group}", fontsize=18, fontweight='bold')
        plt.xticks(rotation=45, ha='right', fontsize=10, fontweight='bold')
        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.savefig("{}/rateAnalysis_medianDiff_{}.pdf".format(median_diff_dir, group.replace("/", "_")))
        plt.close()

        # Plot ratios
        plt.figure()
        x_labels = [f"{eras_list[era + 1]} / {eras_list[era]}" for era in range(len(eras) - 1)]
        for idx, trig in enumerate(triggers):
            ratios = []
            for era in range(len(eras) - 1):
                era_1 = eras_list[era]
                era_2 = eras_list[era + 1]
                
                if trig in trigger_names:
                    group_avg_rate_era1 = avg_rates[era_1][trigger_names.index(trig)]
                    group_avg_rate_era2 = avg_rates[era_2][trigger_names.index(trig)]
                    ratio = (group_avg_rate_era2 - group_avg_rate_era1) / group_avg_rate_era1 if group_avg_rate_era1 != 0 else np.nan
                    ratios.append(ratio)
            plt.plot(x_labels, ratios, label=trig, marker=markers[idx % len(markers)], markersize=8)
        plt.xlabel("Eras", fontsize=18, fontweight='bold')
        plt.ylabel("Rate Change [%]", fontsize=18, fontweight='bold')
        plt.legend(prop={'size': 14, 'weight': 'bold'})
        plt.ylim(-1, 4)
        plt.title(f"Median Ratios Between Eras - Group: {group}", fontsize=18, fontweight='bold')
        plt.xticks(rotation=45, ha='right', fontsize=10, fontweight='bold')
        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.savefig("{}/rateAnalysis_medianRatio_{}.pdf".format(median_ratio_dir, group.replace("/", "_")))
        plt.close()




    return


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()

    parser.add_argument('trigger_dict')
    parser.add_argument('eras')
    parser.add_argument('outDir')
    parser.add_argument('data') 

    args = parser.parse_args()

    tree = open_tree_any(args.data)
    args.tree = tree

    #print(args.trigger_dict, args.eras, args.outDir, args.data)

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    run_rate_check(args)

