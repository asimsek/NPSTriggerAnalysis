from utils import *
import os

def run_rate_monitoring(args):

	by_run_number = True

	data_path =  args.data
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

	runs = ak.to_numpy(t["run"].array())
	pu = ak.to_numpy(t["pileup"].array())
	golden = ak.to_numpy(t["physics_flag"].array())
	cms_ready = ak.to_numpy(t["cms_ready"].array())
	recorded_lumi = ak.to_numpy(t["recorded_lumi_per_lumisection"].array())
	beams_stable = ak.to_numpy(t["beams_stable"].array())
	delivered_lumi = ak.to_numpy(t["delivered_lumi_per_lumisection"].array())
	integrated_lumi = np.cumsum(recorded_lumi)/1000.0

	mask_pu = pu >= 62
	mask_golden = golden == 1
	mask_cms_ready = cms_ready == 1
	mask_beams_stable = beams_stable == 1
	mask_runs_in_eras = (runs >= int(eras[eras_list[0]][0])) & (runs <= int(eras[eras_list[-1]][1]))
	mask_delivered_lumi = delivered_lumi > 0.1

	mask = mask_pu & mask_golden & mask_cms_ready & mask_beams_stable & mask_runs_in_eras & mask_delivered_lumi # only choose good runs to plot
	runs = runs[mask]
	pu   = pu[mask]

	# ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
	# ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]
	for group in list(trigger_dict.keys()):
		print (group)
		colors =  itertools.cycle(["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"])
		colors2 =  itertools.cycle(["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"])
		max_rates = []
		min_rates = []

		fig, ax = plt.subplots(figsize=(16,10))
		
		#### Plot eras rectangles
		for key, value_ in eras.items():
			value = list(map(int, value_))
			if by_run_number:
				rect = Rectangle((value[0], 0), value[1]-value[0], 10e6, 
								 linewidth=0, facecolor=next(colors2),alpha=0.1,
								label=key)
			else:
				bigger_than_0 = np.where(runs>=value[0])[0]
				bigger_than_1 = np.where(runs>=value[1])[0]
				print("bigger than 0 "+str(bigger_than_0))
				print("bigger than 1"+str(bigger_than_1))
				
				if len(bigger_than_0) == 0:
					value_0 = 0
				else: 
					value_0 = integrated_lumi[np.where(runs>=value[0])[0][0]]
				if len(bigger_than_1) == 0:
					value_1 = np.max(integrated_lumi)
				else:
					value_1 = integrated_lumi[np.where(runs>=value[1])[0][0]]
				print("value 0:"+str(value_0))
				print("value 1:"+str(value_1))
				rect = Rectangle((value_0, 0), value_1-value_0, 10e6, 
					 linewidth=0, facecolor=next(colors2),alpha=0.1,
					label=key)
				
			ax.add_patch(rect)
		
		
		for trigger_name in trigger_dict[group]:
			
			if trigger_name+"_v" not in t.keys():
				print("Attention "+trigger_name+"not in dict -------------------")
				continue 
				
			c = next(colors)
			print("---",trigger_name)
			
			#trigger = t[trigger_name+"_v"].array()/LS_seconds # scale to correct lumi
			trigger = t[trigger_name+"_v"].array()/delivered_lumi*2e34/1e36

			#trigger = np.array(trigger).astype(float)
			#runs = np.array(runs).astype(float)
			#pu = np.array(pu).astype(float)
			#trigger[~mask] = np.nan
			#runs[~mask] = np.nan
			#pu[~mask] = np.nan
			
			trigger = trigger[mask]

			#### Apply a median filter 
			N = 500
			trigger_smoothed = median_filter(trigger, N) # filter out lumi decay

			if by_run_number:
				ax.scatter(runs,
							trigger_smoothed,marker='.',s=1,label=trigger_name,color=c,alpha=1.0,rasterized=True)
				ax.plot(runs,trigger_smoothed,color=c,alpha=1)
			else:
				ax.scatter(integrated_lumi,
							trigger_smoothed,marker='.',s=1,label=trigger_name,color=c,alpha=1.0,rasterized=True)
				ax.plot(integrated_lumi,trigger_smoothed,color=c,alpha=1)
			#print(trigger_smoothed)
			
			#### Store max amd min rates for setting axes later 
			max_rates += [np.nanmax(trigger_smoothed) ]
			if np.sum(trigger_smoothed != 0.0) > 0: 
				
				min_rates += [np.nanmin(trigger_smoothed[trigger_smoothed != 0.0])]
			else:
				min_rates += [0.1] 
		max_rate = np.nanmax(max_rates)
		min_rate = np.nanmin(min_rates)
		
		size = fig.get_size_inches()*fig.dpi
		
		log = False
		if log:
			ymin, ymax = (min_rate,max_rate*10)
			ax.set_ylim(ymin,ymax)
			ax.set_yscale('log')
		else:
			ymin, ymax = (0,max_rate*1.1)
			ax.set_ylim(ymin,ymax)

			
		#### Plot pile-up on second axis 
		ax2 = ax.twinx()
		ax2.scatter(runs,
						median_filter(pu,N),marker='.',s=1,#label="Pileup",
						 color='black',alpha=1,rasterized=True)
		#	ax2.plot(runs[mask],
		#			median_filter(pu[mask],N),color='black',alpha=0.1)
		#		ax2.scatter(runs[mask],
		#			pu[mask],marker='.',s=1,label=key, color='black',alpha=0.5)
		#		ax2.plot(runs[mask],
		#			pu[mask],color='black',alpha=0.5)
		ax2.set_ylim(0,70)

		#### Plot updates 
		# for key,value in updates.items(): 
		#     ax.plot([value,value],[ymin, ymax],'--',label=key)

		#### Sort legend labels by max rate
		handles, labels = ax.get_legend_handles_labels()
		handles = sorted(handles[len(eras.items()):], key=lambda x: max_rates[handles[len(eras.items()):].index(x)],reverse=True)
		labels =  sorted(labels[len(eras.items()):], key=lambda x: max_rates[labels[len(eras.items()):].index(x)],reverse=True)
		
		
		#### Zip together any labels that are duplicated 
		handles2, labels2 = ax2.get_legend_handles_labels()
		newLabels, newHandles = [], []
		for handle, label in zip(handles2, labels2):
			if label not in newLabels:
				newLabels.append(label)
				newHandles.append(handle)
		for handle, label in zip(handles, labels):
			if label not in newLabels:
				newLabels.append(label)
				newHandles.append(handle)

				
		#### Make dot sizes bigger in legend 
		lgnd = plt.legend(newHandles, newLabels,frameon=False, bbox_to_anchor=(-0.1, -0.08),loc='upper left', fontsize=18)
		for lh in range(len(lgnd.legend_handles)): 
			lgnd.legend_handles[lh]._sizes = [150]

		ax.set_xlabel("Run Number")
		ax.set_ylabel(r"Rate at 2.0e34 $cm^{-2} s^{-1}$[Hz]")
		ax2.set_ylabel('Pileup')
		plt.title(group+" Trigger Paths", size=24,fontweight="bold")

		# Create the secondary x-axis for era labels
		ax_upper = ax.secondary_xaxis('top')
		ax_upper.set_xlabel('')

		# Add vertical lines and upper labels using the new function
		x_range = (eras[eras_list[0]][0], eras[eras_list[-1]][1])
		add_era_lines(ax, ax_upper, eras, x_range)


		#plt.savefig("axo_nominal.png")

	multipage(outDir+"/SUSTriggerMonitoring_AllCombined.pdf", figs=plt.get_fignums(), dpi=50)

	return


if __name__ == "__main__": 
	parser = argparse.ArgumentParser()

	parser.add_argument('trigger_dict')           # positional argument
	parser.add_argument('eras')
	parser.add_argument('outDir')
	parser.add_argument('data') 

	args = parser.parse_args()
	print(args.trigger_dict, args.eras, args.data)

	if not os.path.exists(args.outDir):
		os.makedirs(args.outDir)

	run_rate_monitoring(args) 


