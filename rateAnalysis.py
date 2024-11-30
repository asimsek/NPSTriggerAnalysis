from utils import *
import os

def run_rate_check(args):

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

    # Check that lumi-leveling looks okay 
    plt.figure()
    for era in eras:
        avg_rates[era] = []
        std_rates[era] = []
        avg_pu[era] = []
        plotted_triggers[era] = []
        print("Era: ",era)
        runs = t["run"].array()
        lumi = t["lumi"].array()
        pileup = t["pileup"].array()
        mask = (runs == eras[era]["run"] ) & (lumi > eras[era]["ls"][0]) & (lumi < eras[era]["ls"][1])
        print("Run Numbers: ", runs[mask])
        for trigger_name in trigger_names:
            if trigger_name+"_v" in t.keys():
                trigger = t[trigger_name+"_v"].array()
                trigger_scaled = trigger/LS_seconds
                plt.scatter(lumi[mask],trigger_scaled[mask],marker='.',s=10,label=era)
                avg_rates[era] +=  [np.mean(trigger_scaled[mask])]
                std_rates[era] += [np.mean(trigger_scaled[mask])]
                avg_pu[era] += [np.mean(pileup[mask])]
                plotted_triggers[era] += [trigger_name]
    #plt.legend(frameon=False, bbox_to_anchor=(1, 0.5))
    plt.xlabel("Lumi-section")
    plt.ylabel("Rate")
    #plt.savefig(outDir+"/rates_during_rate_checks.pdf")

    # Plot differences 
    plt.figure()
    era_names = list(eras.keys())
    for era in range(len(eras)-1):
        era_1 = era_names[era]
        era_2 = era_names[era+1]
        print(era_1+"-"+era_2+"---------------------------")
        difference = np.array(avg_rates[era_2])-np.array(avg_rates[era_1])
        #mean = np.mean(difference[~np.isnan(difference) & ~np.isinf(difference)])
        mean = np.median(difference[~np.isnan(difference) & ~np.isinf(difference)])
        plt.hist(difference,bins=100,range=(-10,10),
                 alpha=0.3,
                label=era_2+" - "+era_1+"\n Median: "+"{:.2f}".format(mean))
        mask = (difference > 10) | (difference < -10)
        for i in range(np.sum(mask)):  
            print(np.array(plotted_triggers[era_1])[mask][i]+ ": "+"{:.0f}".format(difference[mask][i])+" Hz")
    #         print(np.array(plotted_triggers[era_1])[mask][i])
    #         print("Difference:"+"{:.2f}".format(difference[mask][i])+" Hz")
    #         print(era_1+": "+"{:.2f}".format(np.array(avg_rates[era_1])[mask][i])+" Hz")
    #         print(era_2+": "+"{:.2f}".format(np.array(avg_rates[era_2])[mask][i])+" Hz")
    plt.xlabel("Difference Rate in Era 2 - Rate in Era 1 [Hz]")
    plt.legend(fontsize=16)
    plt.ylabel("Number of trigger paths")
    #plt.yscale('log')
    plt.savefig(outDir+"/rate_check_differences.pdf")

    # Make ratio plots of rates
    plt.figure()
    for era in range(len(eras)-1):
        era_1 = era_names[era]
        era_2 = era_names[era+1]
        print(era_1+"-"+era_2+"---------------------------")
        ratio = np.array(avg_rates[era_2])/np.array(avg_rates[era_1])
        #mean = np.mean(difference[~np.isnan(difference) & ~np.isinf(difference)])
        mean = np.median(ratio[~np.isnan(ratio) & ~np.isinf(ratio)])
        plt.hist(ratio,bins=100,range=(0,5),
                 alpha=0.3,
                label=era_2+" / "+era_1+"\n Median: "+"{:.2f}".format(mean))
        mask = (ratio > 2) | (ratio < 0.5)
        for i in range(np.sum(mask)):
            print(np.array(plotted_triggers[era_1])[mask][i])
            print("Ratio:"+"{:.2f}".format(ratio[mask][i])+" Hz")
            print(era_1+": "+"{:.2f}".format(np.array(avg_rates[era_1])[mask][i])+" Hz")
            print(era_2+": "+"{:.2f}".format(np.array(avg_rates[era_2])[mask][i])+" Hz")
    plt.xlabel("Ratio of Rate in Era 2 / Rate in Era 1")
    plt.legend(fontsize=16)
    plt.ylabel("Number of trigger paths")
    #plt.yscale('log')
    #outDir = "plots/STEAM_Oct2024"
    plt.savefig(outDir+"/rateAnalysis_ratios.pdf")
    

    return


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()

    parser.add_argument('trigger_dict')           # positional argument
    parser.add_argument('eras')
    parser.add_argument('outDir')
    parser.add_argument('data') 
    

    args = parser.parse_args()
    print(args.trigger_dict, args.eras, args.outDir, args.data)

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    run_rate_check(args) 

