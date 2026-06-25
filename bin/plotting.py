from bin.utils import *
from bin.brokenaxes import brokenaxes
from datetime import datetime
from dateutil.relativedelta import relativedelta

matplotlib.use('Agg')

def reset_colors():
    return itertools.cycle(["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"])

def get_gap_threshold(x_data, x_axis):
    x_data = np.asarray(x_data, dtype=float)
    if len(x_data) < 2:
        return np.inf

    diffs = np.diff(x_data)
    diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
    if len(diffs) == 0:
        return np.inf

    typical_step = np.nanmedian(diffs)
    if x_axis == 'date':
        return max(3.0, 5.0 * typical_step)
    return max(100.0, 20.0 * typical_step)

def get_continuous_segments(x_data, x_axis):
    x_data = np.asarray(x_data, dtype=float)
    if len(x_data) == 0:
        return []

    gap_threshold = get_gap_threshold(x_data, x_axis)
    finite = np.isfinite(x_data)
    break_mask = (
        (~finite[:-1]) |
        (~finite[1:]) |
        (np.diff(x_data) > gap_threshold)
    )
    starts = np.r_[0, np.where(break_mask)[0] + 1]
    stops = np.r_[starts[1:], len(x_data)]
    return [(start, stop) for start, stop in zip(starts, stops) if stop > start]

def median_filter_segments(x_data, y_data, x_axis, nominal_window=500):
    y_data = np.asarray(y_data, dtype=float)
    y_smoothed = np.full_like(y_data, np.nan, dtype=float)
    segments = get_continuous_segments(x_data, x_axis)

    for start, stop in segments:
        segment = y_data[start:stop]
        segment_len = stop - start
        if segment_len < 3:
            y_smoothed[start:stop] = segment
            continue

        if segment_len >= nominal_window:
            window = nominal_window
        else:
            window = max(3, segment_len // 2)

        if window % 2 == 0:
            window -= 1
        if window < 3:
            y_smoothed[start:stop] = segment
        else:
            y_smoothed[start:stop] = median_filter(segment, size=window)

    return y_smoothed, segments

def plot_rate(group, trigger_dict, t, delivered_lumi, mask, eras, era_dates, x_axis, x_label, runs, dates, pu, eras_list, print_trig, rate_cache=None, rate_cache_is_masked=False, plot_cache=None):
    fig = plt.figure(figsize=(16, 10))
    min_rates, max_rates = [], []
    colors = reset_colors()

    # ------------------------------------------------------------------
    # Remove bad run ranges from runs/dates
    # ------------------------------------------------------------------
    good_mask = get_good_run_mask(runs)

    # Apply to x-axis helper arrays
    runs = runs[good_mask]
    dates = dates[good_mask]
    # pu = pu[good_mask]

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
            (399000, 401500),
            # Add more break areas here, e.g., (400000, 405000)
        ],
        'date': [
            (datetime(2024, 10, 25), datetime(2025, 5, 10)),  # Break between 2024 and 2025
            (datetime(2025, 6, 21), datetime(2025, 7, 10)),
            (datetime(2025, 11, 1), datetime(2026, 3, 1)),
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
        c = next(colors)
        if print_trig:
            print("---", trigger_name)

        cache_key = (x_axis, trigger_name)
        cached = plot_cache.get(cache_key) if plot_cache is not None else None
        if cached is not None:
            x_data = cached["x_data"]
            trigger_smoothed = cached["trigger_smoothed"]
            segments = cached["segments"]
        else:
            if rate_cache is not None:
                if trigger_name not in rate_cache:
                    continue
                trigger = rate_cache[trigger_name]
            else:
                possible_branches = (f"{trigger_name}_v", trigger_name)
                branch_name = next((b for b in possible_branches if b in t.keys()), None)

                if branch_name is None:
                    print(f" -- Attention {trigger_name} not in tree!")
                    continue

                trigger = t[branch_name].array() / delivered_lumi * 2e34 / 1e36
            if not rate_cache_is_masked:
                trigger = trigger[mask]

            # Apply the same bad-run filter as we used for runs/dates
            trigger = trigger[good_mask]

            if x_axis == 'date':
                # dates has already been filtered with good_mask
                trigger_dates = dates[:len(trigger)]
                valid_mask = ~pd.isna(trigger_dates)
                trigger = trigger[valid_mask]
                trigger_dates = trigger_dates[valid_mask]
                x_data = mdates.date2num(trigger_dates)
                trigger = trigger[:len(x_data)]
            else:
                # runs has already been filtered with good_mask
                x_data = runs

            # Apply median filtering only inside continuous data segments.
            # Otherwise long no-data gaps are connected by artificial flat lines.
            trigger_smoothed, segments = median_filter_segments(x_data, trigger, x_axis)
            if plot_cache is not None:
                plot_cache[cache_key] = {
                    "x_data": x_data,
                    "trigger_smoothed": trigger_smoothed,
                    "segments": segments,
                }

        # Plot only the part of each segment that belongs to each broken-axis
        # panel. Passing off-panel data to every axis can create clipped
        # horizontal artifacts.
        label_used = False
        for start_idx, stop_idx in segments:
            if stop_idx <= start_idx:
                continue
            for ax in bax.axs:
                x_min, x_max = sorted(ax.get_xlim())
                x_segment = x_data[start_idx:stop_idx]
                y_segment = trigger_smoothed[start_idx:stop_idx]
                visible_mask = (
                    np.isfinite(x_segment) &
                    np.isfinite(y_segment) &
                    (x_segment >= x_min) &
                    (x_segment <= x_max)
                )
                if np.count_nonzero(visible_mask) < 2:
                    continue
                ax.plot(
                    x_segment[visible_mask],
                    y_segment[visible_mask],
                    color=c,
                    alpha=1.0,
                    label=trigger_name if not label_used else None
                )
                label_used = True

        finite_smoothed = trigger_smoothed[np.isfinite(trigger_smoothed)]
        nonzero_smoothed = finite_smoothed[finite_smoothed != 0.0]
        if len(finite_smoothed) == 0:
            continue

        max_rates.append(np.nanmax(finite_smoothed))
        if len(nonzero_smoothed) > 0:
            min_rates.append(np.nanmin(nonzero_smoothed))
        else:
            min_rates.append(0.1)

        trigger_to_rate[trigger_name] = np.nanmax(finite_smoothed)

    max_rate = np.nanmax(max_rates) if max_rates else 1.0
    min_rate = np.nanmin(min_rates) if min_rates else 0.1

    log = False
    ymin, ymax = (0, max_rate * 1.3) if not log else (min_rate, max_rate * 10)
    bax.set_ylim(ymin, ymax)
    if log:
        bax.set_yscale('log')

    # Sort legend labels by max rate
    handles, labels = [], []
    for ax in bax.axs:
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        handles.extend(ax_handles)
        labels.extend(ax_labels)

    seen_labels = set()
    trigger_rate_pairs = []
    for handle, label in zip(handles, labels):
        if label in trigger_to_rate and label not in seen_labels:
            trigger_rate_pairs.append((handle, label))
            seen_labels.add(label)
    trigger_rate_pairs.sort(key=lambda pair: trigger_to_rate[pair[1]], reverse=True)
    newHandles, newLabels = zip(*trigger_rate_pairs) if trigger_rate_pairs else ([], [])

    # Create legend
    lgnd = bax.legend(newHandles, newLabels,
                      title=(group if all(n.startswith("L1_") for n in trigger_dict[group]) else group + " Trigger Paths"),
                      frameon=False, bbox_to_anchor=(-0.1, -0.08), loc='upper left', fontsize=18)
    lgnd.get_title().set_ha('left')
    legend_handles = getattr(lgnd, "legend_handles", None)
    if legend_handles is None:
        legend_handles = getattr(lgnd, "legendHandles", [])
    for lh in legend_handles:
        if hasattr(lh, "_sizes"):
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
            if key not in era_dates:
                continue
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

