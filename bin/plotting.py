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
    colors = reset_colors()

    # Determine data range
    if x_axis == 'run':
        start, end = runs.min(), runs.max()
    elif x_axis == 'date':
        start_date, end_date = dates.min(), dates.max()
        start, end = mdates.date2num(start_date), mdates.date2num(end_date)

    # Define potential break areas (modify these as needed)
    break_areas = {
        'run': [
            (387300, 392000),  # Example break area for runs
            (393400, 394500),
            # Add more break areas here, e.g., (400000, 405000)
        ],
        'date': [
            (datetime(2024, 10, 25), datetime(2025, 5, 10)),  # Break between 2024 and 2025
            (datetime(2025, 6, 21), datetime(2025, 7, 10)),
            # Add more break areas here, e.g., (datetime(2024, 6, 1), datetime(2024, 7, 1))
        ]
    }

    # Filter relevant breaks based on data range
    def get_valid_breaks(start, end, breaks, is_date=False):
        valid_breaks = []
        for break_start, break_end in breaks:
            break_start_num = mdates.date2num(break_start) if is_date else break_start
            break_end_num = mdates.date2num(break_end) if is_date else break_end
            # Include break if it overlaps with the data range
            if break_start_num < end and break_end_num > start:
                valid_breaks.append((break_start_num, break_end_num))
        return valid_breaks

    # Get valid breaks for the current axis
    is_date = x_axis == 'date'
    valid_breaks = get_valid_breaks(start, end, break_areas[x_axis], is_date)

    # Create xlims for brokenaxes
    if valid_breaks:
        # Sort breaks by start point
        valid_breaks.sort()
        # Create xlims: from start to first break, between breaks, and from last break to end
        xlims = []
        prev_end = start
        for break_start, break_end in valid_breaks:
            if prev_end < break_start:
                xlims.append((prev_end, break_start))
            prev_end = break_end
        if prev_end < end:
            xlims.append((prev_end, end))
        # Add padding for better visualization
        if is_date:
            xlims = [(mdates.date2num(start_date - relativedelta(days=3)), xlims[0][1]) if i == 0 else x for i, x in enumerate(xlims)]
            xlims[-1] = (xlims[-1][0], mdates.date2num(end_date + relativedelta(days=3)))
        else:
            xlims = [(start - 100, xlims[0][1]) if i == 0 else x for i, x in enumerate(xlims)]
            xlims[-1] = (xlims[-1][0], end + 100)
    else:
        # No valid breaks, use full range
        xlims = [(start - (100 if not is_date else mdates.date2num(start_date - relativedelta(days=3))),
                  end + (100 if not is_date else mdates.date2num(end_date + relativedelta(days=3))))]

    # Initialize brokenaxes or standard axis
    if len(xlims) > 1:
        bax = brokenaxes(
            xlims=xlims,
            hspace=0.2,
            despine=False,
            d=0.0,
            diag_color='none'
        )
    else:
        # Use standard matplotlib axis if no breaks are needed
        bax = type('Bax', (), {'axs': [plt.gca()], 'set_ylim': lambda self, ymin, ymax: plt.gca().set_ylim(ymin, ymax),
                               'set_yscale': lambda self, scale: plt.gca().set_yscale(scale),
                               'legend': lambda self, *args, **kwargs: plt.gca().legend(*args, **kwargs),
                               'axs': [plt.gca()]})()

    # Set y-axis properties
    for ax in bax.axs:
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=10))

    trigger_to_rate = {}
    for trigger_name in trigger_dict[group]:
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
            trigger = trigger[:len(x_data)]
        else:
            x_data = runs

        # Apply median filter
        N = 500
        trigger_smoothed = median_filter(trigger, N)

        # Plot on broken axes
        for ax in bax.axs:
            ax.plot(x_data, trigger_smoothed, color=c, alpha=1.0, label=trigger_name if ax == bax.axs[0] else None)

        max_rates.append(np.nanmax(trigger_smoothed))
        if np.sum(trigger_smoothed != 0.0) > 0:
            min_rates.append(np.nanmin(trigger_smoothed[trigger_smoothed != 0.0]))
        else:
            min_rates.append(0.1)

        trigger_to_rate[trigger_name] = np.nanmax(trigger_smoothed)

    max_rate = np.nanmax(max_rates) if max_rates else 1.0
    min_rate = np.nanmin(min_rates) if min_rates else 0.1

    log = False
    ymin, ymax = (0, max_rate * 1.3) if not log else (min_rate, max_rate * 10)
    bax.set_ylim(ymin, ymax)
    if log:
        bax.set_yscale('log')

    # Sort legend labels by max rate
    handles, labels = bax.axs[0].get_legend_handles_labels()
    trigger_rate_pairs = [(h, l) for h, l in zip(handles, labels) if l in trigger_to_rate]
    trigger_rate_pairs.sort(key=lambda pair: trigger_to_rate[pair[1]], reverse=True)
    newHandles, newLabels = zip(*trigger_rate_pairs) if trigger_rate_pairs else ([], [])

    # Create legend
    lgnd = bax.legend(newHandles, newLabels,
                      title=(group if all(n.startswith("L1_") for n in trigger_dict[group]) else group + " Trigger Paths"),
                      frameon=False, bbox_to_anchor=(-0.1, -0.08), loc='upper left', fontsize=18)
    lgnd.get_title().set_ha('left')
    for lh in lgnd.legend_handles:
        lh._sizes = [150]

    # Set labels
    bax.axs[-1].set_xlabel(x_label)
    bax.axs[0].set_ylabel(r"Rate at 2.0e34 $cm^{-2} s^{-1}$[Hz]")

    # Add annotations and vertical lines for eras
    for idx, key in enumerate(eras_list):
        if x_axis == 'run':
            start_run = runs.min() if key == eras_list[0] else float(eras[key][0])
            end_run = runs.max() if key == eras_list[-1] else float(eras[key][1])
            mid_run = (start_run + end_run) / 2

            for ax in bax.axs:
                xlim = ax.get_xlim()
                if xlim[0] <= mid_run <= xlim[1]:
                    ax.annotate(key, xy=(mid_run, ymax), xycoords='data', textcoords='offset points',
                                xytext=(0, -10), ha='center', va='top',
                                fontsize=12, color='black', weight='bold', rotation=90)
                if idx < len(eras_list) - 1 and end_run >= xlim[0] and end_run <= xlim[1]:
                    ax.axvline(end_run, linestyle='--', linewidth=1, color='black', alpha=0.8)

        elif x_axis == 'date':
            start_date = dates.min() if key == eras_list[0] else era_dates[key][0]
            end_date = era_dates[key][1]
            start_num = mdates.date2num(start_date)
            end_num = mdates.date2num(end_date)
            mid_date = (start_num + end_num) / 2

            for ax in bax.axs:
                xlim = ax.get_xlim()
                if xlim[0] <= mid_date <= xlim[1]:
                    ax.annotate(key, xy=(mid_date, ymax), xycoords='data', textcoords='offset points',
                                xytext=(0, -10), ha='center', va='top',
                                fontsize=12, color='black', weight='bold', rotation=90)
                if idx < len(eras_list) - 1 and end_num >= xlim[0] and end_num <= xlim[1]:
                    ax.axvline(end_num, linestyle='--', linewidth=1, color='black', alpha=0.8)

    # Configure axis ticks and format
    for i, ax in enumerate(bax.axs):
        ax.tick_params(direction='in', which='both')
        if i == 0:
            # Leftmost axis: show left y-tickmarks (and labels), hide right y-tickmarks
            ax.tick_params(left=True, labelleft=True, right=False, labelright=False, top=True, bottom=True, labelbottom=True)
            ax.spines['right'].set_visible(False)
        elif i == len(bax.axs) - 1:
            # Rightmost axis: hide left y-tickmarks, show right y-tickmarks (and labels)
            ax.tick_params(left=False, labelleft=False, right=True, labelright=True, top=True, bottom=True, labelbottom=True)
            ax.spines['left'].set_visible(False)
        else:
            # Middle axes: hide both left and right y-tickmarks
            ax.tick_params(left=False, labelleft=False, right=False, labelright=False, top=True, bottom=True, labelbottom=True)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)

        # Preserve your existing top spine settings if desired (though they may be unnecessary for vertical axes; consider removing if not needed)
        ax.spines['top'].set_visible(True)
        ax.spines['top'].set_linewidth(2)

        if x_axis == 'date':
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
            ax.xaxis.set_tick_params(rotation=45, labelsize=16)

    plt.close(fig)
    return fig


