from bin.utils import *

# Set up global color iterators for consistent coloring across all plots
colors = itertools.cycle(["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"])
colors2 = itertools.cycle(["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"])

matplotlib.use('Agg')

def plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis, x_label, runs, dates, pu, eras_list, print_trig):
    fig, ax = plt.subplots(figsize=(16, 10))
    plt.subplots_adjust(top=0.85)  # Add a margin on top to visualize era labels properly
    max_rates = []
    min_rates = []

    # Plot eras rectangles
    for key, value_ in eras.items():
        value = list(map(int, value_))
        if x_axis == 'run':
            # Extend the first and last era to the plot boundaries
            if key == eras_list[0]:
                start = runs.min()
            else:
                start = value[0]
            if key == eras_list[-1]:
                end = runs.max()
            else:
                end = value[1]
            rect = Rectangle((start, 0), end - start, 10e6,
                             linewidth=0, facecolor=next(colors2), alpha=0.1, label=key)
            ax.add_patch(rect)
        elif x_axis == 'date':
            if era_dates[key][0] and era_dates[key][1]:
                # Extend the first and last era to the plot boundaries
                if key == eras_list[0]:
                    start_date = dates.min()
                else:
                    start_date = era_dates[key][0]
                if key == eras_list[-1]:
                    end_date = dates.max()
                else:
                    end_date = era_dates[key][1]
                rect_date = Rectangle((mdates.date2num(start_date), 0),
                                      mdates.date2num(end_date) - mdates.date2num(start_date),
                                      10e6, linewidth=0, facecolor=next(colors2), alpha=0.1, label=key)
                ax.add_patch(rect_date)

    for trigger_name in trigger_dict[group]:
        if trigger_name + "_v" not in t.keys():
            print("Attention " + trigger_name + " not in dict -------------------")
            continue

        c = next(colors)
        if print_trig:
            print("---", trigger_name)

        trigger = t[trigger_name + "_v"].array() / delivered_lumi * 2e34 / 1e36
        trigger = trigger[mask]

        # Apply a median filter
        N = 500
        trigger_smoothed = median_filter(trigger, N) # filter out lumi decay
        #trigger_smoothed = savgol_filter(trigger, window_length=500, polyorder=3)

        x_data = runs if x_axis == 'run' else dates[:len(trigger)]
        ax.scatter(x_data, trigger_smoothed, marker='.', s=1, label=trigger_name, color=c, alpha=1.0, rasterized=True)
        ax.plot(x_data, trigger_smoothed, color=c, alpha=1)

        # Store max and min rates for setting axes later
        max_rates += [np.nanmax(trigger_smoothed)]
        if np.sum(trigger_smoothed != 0.0) > 0:
            min_rates += [np.nanmin(trigger_smoothed[trigger_smoothed != 0.0])]
        else:
            min_rates += [0.1]

    max_rate = np.nanmax(max_rates)
    min_rate = np.nanmin(min_rates)

    log = False
    if log:
        ymin, ymax = (min_rate, max_rate * 10)
        ax.set_ylim(ymin, ymax)
        ax.set_yscale('log')
    else:
        ymin, ymax = (0, max_rate * 1.3)
        ax.set_ylim(ymin, ymax)

    # Plot pile-up on second axis
    ax2 = ax.twinx()
    x_data_pu = runs if x_axis == 'run' else dates

    #ax2.scatter(x_data_pu, savgol_filter(pu, window_length=51, polyorder=2), marker='.', s=1, color='black', alpha=1, rasterized=True)
    #ax2.scatter(x_data_pu, median_filter(pu, 50), marker='.', s=1, color='black', alpha=1, rasterized=True)
    ax2.scatter(x_data_pu, pu, marker='.', s=1, color='black', alpha=1, rasterized=True)
    ax2.set_ylim(0, 75)

    # Sort legend labels by max rate
    handles, labels = ax.get_legend_handles_labels()
    handles = sorted(handles[len(eras.items()):], key=lambda x: max_rates[handles[len(eras.items()):].index(x)], reverse=True)
    labels = sorted(labels[len(eras.items()):], key=lambda x: max_rates[labels[len(eras.items()):].index(x)], reverse=True)

    # Zip together any labels that are duplicated
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)

    # Make dot sizes bigger in legend
    lgnd = plt.legend(newHandles, newLabels, title=group + " Trigger Paths", frameon=False, bbox_to_anchor=(-0.1, -0.08), loc='upper left', fontsize=18)
    lgnd.get_title().set_ha('left')  # Set the title alignment to left
    for lh in range(len(lgnd.legend_handles)):
        lgnd.legend_handles[lh]._sizes = [150]

    ax.set_xlabel(x_label)
    ax.set_ylabel(r"Rate at 2.0e34 $cm^{-2} s^{-1}$[Hz]")
    ax2.set_ylabel('Pileup')

    # Create the secondary x-axis for era labels
    ax_upper = ax.secondary_xaxis('top')
    ax_upper.set_xlabel('')
    if x_axis == 'run':
        add_era_lines(ax, ax_upper, eras, (runs.min(), runs.max()))
        ax_upper.xaxis.set_major_locator(plt.MaxNLocator(nbins=10))  # Ensure consistent tick division for run plots
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=10))
    elif x_axis == 'date':
        add_era_lines(ax, ax_upper, era_dates, (dates.min(), dates.max()))
        ax_upper.xaxis.set_major_locator(mdates.AutoDateLocator())  # Ensure consistent tick division for date plots
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())

    # Move era labels just below the top axis line (inside the plot)
    ax_upper.tick_params(axis='x', pad=-20, labelsize=10, labelrotation=45)  # Adjust padding to move labels inside plot

    # Add annotations for each era label manually, centering them
    for key, (start_run, end_run) in eras.items():
        if x_axis == 'run':
            start_run = runs.min() if key == eras_list[0] else float(start_run)
            end_run = runs.max() if key == eras_list[-1] else float(end_run)
            if start_run >= runs.min() and end_run <= runs.max():
                ax.annotate(key, xy=((start_run + end_run) / 2, ymax), xycoords='data', textcoords='offset points',
                            xytext=(0, -10), ha='center', va='top', fontsize=12, color='black', weight='bold', rotation=90)
        elif x_axis == 'date':
            start_date = mdates.date2num(dates.min()) if key == eras_list[0] else mdates.date2num(era_dates[key][0])
            end_date = mdates.date2num(dates.max()) if key == eras_list[-1] else mdates.date2num(era_dates[key][1])
            if start_date >= mdates.date2num(dates.min()) and end_date <= mdates.date2num(dates.max()):
                ax.annotate(key, xy=((start_date + end_date) / 2, ymax), xycoords='data', textcoords='offset points',
                            xytext=(0, -10), ha='center', va='top', fontsize=12, color='black', weight='bold', rotation=90)

    plt.setp(ax_upper.get_xticklabels(), ha='right')
    plt.xticks(rotation=45, ha='right', fontsize=8)

    plt.close(fig)
    return fig
