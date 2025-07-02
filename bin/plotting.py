from bin.utils import *
from bin.brokenaxes import brokenaxes
from datetime import datetime
from dateutil.relativedelta import relativedelta


matplotlib.use('Agg')

def reset_colors():
    return itertools.cycle(["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"])

def plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis, x_label, runs, dates, pu, eras_list, print_trig):
    fig = plt.figure(figsize=(16, 10))
    min_rates, max_rates = [], []
    colors, colors2 = reset_colors(), reset_colors()
    start, end = runs.min(), runs.max()

    if x_axis == 'run':
        bax = brokenaxes(
            xlims=((start-100, 387300), (392000, end+100)),
            hspace=0.2,
            despine=False,
            d=0.0,
            diag_color='none'
    )

    elif x_axis == 'date':
            start_date, end_date = dates.min(), dates.max()
            start_date_extended, end_date_extended = (start_date - relativedelta(days=3)), (end_date + relativedelta(days=3))
           
            bax = brokenaxes(
                xlims=(
                    ( start_date_extended, datetime(2024, 10, 25) ),
                    ( datetime(2025, 5, 10), end_date_extended )
                ),
            hspace=0.2, despine=False, d=0.0, diag_color='none')

    bax.axs[0].yaxis.set_major_locator(plt.MaxNLocator(nbins=10))
    bax.axs[1].yaxis.set_major_locator(plt.MaxNLocator(nbins=10))



    trigger_to_rate = {}
    for trigger_name in trigger_dict[group]:
        #if trigger_name + "_v" not in t.keys():
        #    print("Attention " + trigger_name + " not in dict -------------------")
        #    continue

        possible_branches = (f"{trigger_name}_v", trigger_name)
        branch_name = next((b for b in possible_branches if b in t.keys()), None)

        if branch_name is None:
            print(f" -- Attention {trigger_name} not in tree!")
            continue

        c = next(colors)
        if print_trig:
            print("---", trigger_name)

        trigger = t[branch_name].array() / delivered_lumi * 2e34 / 1e36
        trigger = trigger[mask]
        if x_axis == 'date':
            trigger_dates = dates[:len(trigger)]
            valid_mask = ~pd.isna(trigger_dates)
            trigger = trigger[valid_mask]
            trigger_dates = trigger_dates[valid_mask]
            
            x_data = mdates.date2num(trigger_dates)
            trigger = trigger[:len(x_data)]  # Ensure equal lengths

        else:
            x_data = runs


        #trigger = pd.Series(trigger).interpolate(limit_direction="both", limit_area="inside").to_numpy()
        # Apply a median filter
        N = 500
        trigger_smoothed = median_filter(trigger, N) # filter out lumi decay
        #trigger_smoothed = savgol_filter(trigger, window_length=500, polyorder=3)

        bax.plot(x_data, trigger_smoothed, color=c, alpha=1.0, label=trigger_name) #, rasterized=True)

        # Store max and min rates for setting axes later
        max_rates += [np.nanmax(trigger_smoothed)]
        if np.sum(trigger_smoothed != 0.0) > 0:
            min_rates += [np.nanmin(trigger_smoothed[trigger_smoothed != 0.0])]
        else:
            min_rates += [0.1]

        trigger_to_rate[trigger_name] = np.nanmax(trigger_smoothed)

    max_rate = np.nanmax(max_rates)
    min_rate = np.nanmin(min_rates)

    log = False
    ymin, ymax = 0, 0
    if log:
        ymin, ymax = (min_rate, max_rate * 10)
        #bax.set_ylim(ymin, ymax)
        bax.set_yscale('log')
    else:
        ymin, ymax = (0, max_rate * 1.3)
    
    bax.set_ylim(ymin, ymax)


    # Sort legend labels by max rate
    handles, labels = bax.axs[0].get_legend_handles_labels()

    trigger_rate_pairs = [(h, l) for h, l in zip(handles, labels) if l in trigger_to_rate]
    trigger_rate_pairs.sort(key=lambda pair: trigger_to_rate[pair[1]], reverse=True)

    newHandles, newLabels = zip(*trigger_rate_pairs) if trigger_rate_pairs else ([], [])

    # Make dot sizes bigger in legend
    lgnd = bax.legend(newHandles, newLabels, 
        title=(group if all(n.startswith("L1_") for n in trigger_dict[group]) else group + " Trigger Paths"), 
        frameon=False, bbox_to_anchor=(-0.1, -0.08), loc='upper left', fontsize=18)
    lgnd.get_title().set_ha('left')
    for lh in range(len(lgnd.legend_handles)):
        lgnd.legend_handles[lh]._sizes = [150]

    bax.axs[1].set_xlabel(x_label)
    bax.axs[0].set_ylabel(r"Rate at 2.0e34 $cm^{-2} s^{-1}$[Hz]")

    # Add annotations and vertical lines for each era
    for idx, key in enumerate(eras_list):
        if x_axis == 'run':
            start_run = runs.min() if key == eras_list[0] else float(eras[key][0])
            end_run = runs.max() if key == eras_list[-1] else float(eras[key][1])
            mid_run = (start_run + end_run) / 2

            for i, subax in enumerate(bax.axs):
                xlim = subax.get_xlim()

                # Annotate only if visible
                if xlim[0] <= mid_run <= xlim[1]:
                    subax.annotate(key, xy=(mid_run, ymax), xycoords='data', textcoords='offset points',
                                   xytext=(0, -10), ha='center', va='top',
                                   fontsize=12, color='black', weight='bold', rotation=90)

                # Vertical line logic
                if idx < len(eras_list) - 1:
                    if i == 0 and xlim[0] <= end_run <= xlim[1]:
                        subax.axvline(end_run, linestyle='--', linewidth=1, color='black', alpha=0.8)
                    elif i == 1:
                        next_start_run = float(eras[eras_list[idx + 1]][0])
                        if xlim[0] <= next_start_run <= xlim[1]:
                            subax.axvline(next_start_run, linestyle='--', linewidth=1, color='black', alpha=0.8)

        elif x_axis == 'date':
            start_date = dates.min() if key == eras_list[0] else era_dates[key][0]
            end_date   = era_dates[key][1]

            start_num = mdates.date2num(start_date)
            end_num   = mdates.date2num(end_date)
            mid_date  = (start_num + end_num) / 2

            for i, subax in enumerate(bax.axs):
                xlim = subax.get_xlim()

                # Annotate only if mid_date is visible in this subaxis
                if xlim[0] <= mid_date <= xlim[1]:
                    subax.annotate(key, xy=(mid_date, ymax), xycoords='data', textcoords='offset points',
                                   xytext=(0, -10), ha='center', va='top',
                                   fontsize=12, color='black', weight='bold', rotation=90)

                # Draw vertical line: end_date in 1st subax, next start_date in 2nd
                if idx < len(eras_list) - 1:
                    if i == 0 and xlim[0] <= end_num <= xlim[1]:
                        subax.axvline(end_num, linestyle='--', linewidth=1, color='black', alpha=0.8)
                    elif i == 1:
                        next_start = era_dates[eras_list[idx + 1]][0]
                        next_start_num = mdates.date2num(next_start)
                        if xlim[0] <= next_start_num <= xlim[1]:
                            subax.axvline(next_start_num, linestyle='--', linewidth=1, color='black', alpha=0.8)


            #print(f"{key} → start={start_date.date()}  end={end_date.date()}  mid={(start_date + (end_date - start_date)/2).date()}")




    #bax.axs[0].yaxis.set_major_locator(plt.MaxNLocator(nbins=10))
    for i, subax in enumerate(bax.axs):
        subax.tick_params(top=True, right=True, direction='in', which='both')
        if i < len(bax.axs) - 1:
            subax.tick_params(labelbottom=True, left=True, bottom=True, right=False, top=True, direction='in', which='both')
            subax.spines['top'].set_visible(True)
            subax.spines['top'].set_linewidth(2)
        if i > 0:
            subax.tick_params(labelbottom=True, labelleft=False, left=False, bottom=True, right=True, top=True, direction='in', which='both')
            subax.spines['top'].set_visible(True)
            subax.spines['top'].set_linewidth(2)
            subax.spines['right'].set_visible(True)
            subax.spines['right'].set_linewidth(2)


    if x_axis == 'date':
        for subax in bax.axs:
            #subax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
            subax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
            subax.xaxis.set_tick_params(rotation=45, labelsize=16)

    plt.close(fig)
    return fig

