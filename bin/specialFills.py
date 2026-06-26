from pathlib import Path

from matplotlib.colors import TwoSlopeNorm

from bin.brokenaxes import brokenaxes
from bin.rateChangeMatrix import _ordered_unique_triggers
from bin.utils import *


DEFAULT_SPECIAL_FILLS_CONFIG = "jsonFiles/specialFills/specialFills.json"


def _expand_fill_selection(case, fills_key="fills", range_key="fillRange"):
    selected = set()
    for fill in case.get(fills_key, []):
        selected.add(int(fill))

    fill_range = case.get(range_key)
    if fill_range is not None:
        if not isinstance(fill_range, list) or len(fill_range) != 2:
            raise ValueError(f"{range_key} must contain [start, end]")
        start_fill, end_fill = map(int, fill_range)
        if end_fill < start_fill:
            start_fill, end_fill = end_fill, start_fill
        selected.update(range(start_fill, end_fill + 1))

    return selected


def load_special_fill_config(config_path):
    config_path = Path(config_path).expanduser()
    if not config_path.exists():
        print(f"\033[91m[WARNING] Special-fill config not found:\033[0m {config_path}")
        return []

    with open(config_path) as fp:
        config = json.load(fp)

    cases = []
    for idx, case in enumerate(config.get("cases", [])):
        label = str(case.get("label", "")).strip()
        if not label:
            print(f"\033[91m[WARNING] Skipping special-fill case #{idx + 1}: missing label.\033[0m")
            continue

        try:
            special_fills = _expand_fill_selection(case)
            reference_fills = _expand_fill_selection(
                case,
                fills_key="referenceFills",
                range_key="referenceFillRange",
            )
        except (TypeError, ValueError) as exc:
            print(f"\033[91m[WARNING] Skipping special-fill case '{label}': {exc}\033[0m")
            continue

        if not special_fills:
            print(f"\033[91m[WARNING] Skipping special-fill case '{label}': no fills configured.\033[0m")
            continue

        cases.append({
            "label": label,
            "fills": special_fills,
            "reference_fills": reference_fills,
            "reference_mode": case.get(
                "referenceMode",
                "sameEraExcludingSpecialFills",
            ),
            "color": case.get("color", f"C{idx % 10}"),
        })

    return cases


def _fill_statistics(rates, fills, min_lumisections):
    rates = np.asarray(rates, dtype=float)
    fills = np.asarray(fills, dtype=int)
    fill_numbers = []
    medians = []
    p10_values = []
    p90_values = []

    for fill in np.unique(fills):
        values = rates[fills == fill]
        values = values[np.isfinite(values)]
        if len(values) < min_lumisections:
            continue
        p10, median, p90 = np.nanpercentile(values, [10, 50, 90])
        if not np.isfinite(median) or median <= 0.0:
            continue
        fill_numbers.append(int(fill))
        medians.append(median)
        p10_values.append(p10)
        p90_values.append(p90)

    return (
        np.asarray(fill_numbers, dtype=int),
        np.asarray(medians, dtype=float),
        np.asarray(p10_values, dtype=float),
        np.asarray(p90_values, dtype=float),
    )


def _continuous_fill_segments(fill_numbers):
    if len(fill_numbers) == 0:
        return []
    if len(fill_numbers) == 1:
        return [(0, 1)]

    positive_steps = np.diff(fill_numbers)
    positive_steps = positive_steps[positive_steps > 0]
    typical_step = np.nanmedian(positive_steps) if len(positive_steps) else 1.0
    gap_threshold = max(20.0, 5.0 * typical_step)
    breaks = np.where(np.diff(fill_numbers) > gap_threshold)[0] + 1
    starts = np.r_[0, breaks]
    stops = np.r_[breaks, len(fill_numbers)]
    return list(zip(starts, stops))


def _automatic_fill_xlims(fill_numbers):
    fill_numbers = np.asarray(fill_numbers, dtype=float)
    fill_numbers = np.unique(fill_numbers[np.isfinite(fill_numbers)])
    fill_numbers.sort()
    if len(fill_numbers) == 0:
        return []
    if len(fill_numbers) == 1:
        fill = fill_numbers[0]
        return [(fill - 1.5, fill + 1.5)]

    gaps = np.diff(fill_numbers)
    positive_gaps = gaps[gaps > 0]
    typical_gap = np.nanmedian(positive_gaps) if len(positive_gaps) else 1.0
    gap_threshold = max(30.0, 8.0 * typical_gap)
    break_indices = np.where(gaps > gap_threshold)[0] + 1
    starts = np.r_[0, break_indices]
    stops = np.r_[break_indices, len(fill_numbers)]

    xlims = []
    for start, stop in zip(starts, stops):
        segment = fill_numbers[start:stop]
        if len(segment) == 1:
            xlims.append((segment[0] - 1.5, segment[0] + 1.5))
            continue

        segment_gaps = np.diff(segment)
        segment_gaps = segment_gaps[segment_gaps > 0]
        segment_typical_gap = np.nanmedian(segment_gaps) if len(segment_gaps) else typical_gap
        padding = max(1.0, min(5.0, 2.0 * segment_typical_gap))
        xlims.append((segment[0] - padding, segment[-1] + padding))

    merged_xlims = []
    for xmin, xmax in xlims:
        if merged_xlims and xmin <= merged_xlims[-1][1]:
            previous_xmin, previous_xmax = merged_xlims[-1]
            merged_xlims[-1] = (previous_xmin, max(previous_xmax, xmax))
        else:
            merged_xlims.append((xmin, xmax))

    return merged_xlims


def _case_patch(case):
    return Rectangle(
        (0, 0),
        1,
        1,
        facecolor=case["color"],
        edgecolor=case["color"],
        alpha=0.16,
        label=case["label"],
    )


def _draw_special_fill_regions(ax, cases):
    visible_cases = []
    x_min, x_max = ax.get_xlim()

    for case in cases:
        case_visible = False
        sorted_fills = sorted(case["fills"])
        if not sorted_fills:
            continue

        intervals = []
        interval_start = sorted_fills[0]
        previous_fill = sorted_fills[0]
        for fill in sorted_fills[1:]:
            if fill == previous_fill + 1:
                previous_fill = fill
                continue
            intervals.append((interval_start, previous_fill))
            interval_start = fill
            previous_fill = fill
        intervals.append((interval_start, previous_fill))

        for start_fill, end_fill in intervals:
            if end_fill < x_min or start_fill > x_max:
                continue
            ax.axvspan(
                start_fill - 0.5,
                end_fill + 0.5,
                color=case["color"],
                alpha=0.16,
                linewidth=0,
                zorder=0,
            )
            case_visible = True

        if case_visible:
            visible_cases.append(case)

    return visible_cases


def _special_fill_era_mask(cases, fills, runs, eras, eras_list):
    fills = np.asarray(fills, dtype=int)
    runs = np.asarray(runs, dtype=int)
    all_special_fills = set()
    for case in cases:
        all_special_fills.update(case["fills"])

    if not all_special_fills:
        return np.zeros(len(fills), dtype=bool), []

    special_runs = runs[np.isin(fills, list(all_special_fills))]
    if len(special_runs) == 0:
        return np.zeros(len(fills), dtype=bool), []

    selected_eras = []
    for era in eras_list:
        start_run, end_run = map(int, eras[era])
        if np.any((special_runs >= start_run) & (special_runs <= end_run)):
            selected_eras.append(era)

    plot_eras = list(selected_eras)
    if len(selected_eras) == 1:
        era_index = eras_list.index(selected_eras[0])
        neighbor_indices = [era_index - 1, era_index, era_index + 1]
        plot_eras = [
            eras_list[index]
            for index in neighbor_indices
            if 0 <= index < len(eras_list)
        ]

    era_mask = np.zeros(len(fills), dtype=bool)
    for era in plot_eras:
        start_run, end_run = map(int, eras[era])
        era_mask |= (runs >= start_run) & (runs <= end_run)

    return era_mask, plot_eras


def _era_fill_ranges(runs, fills, eras, eras_list):
    era_ranges = []
    for era in eras_list:
        start_run, end_run = map(int, eras[era])
        era_mask = (runs >= start_run) & (runs <= end_run)
        era_fills = fills[era_mask]
        era_fills = era_fills[np.isfinite(era_fills)]
        if len(era_fills) == 0:
            continue
        era_ranges.append((era, float(np.nanmin(era_fills)), float(np.nanmax(era_fills))))
    return era_ranges


def _overview_ymax(values, percentile=98.0, outlier_factor=5.0):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values >= 0.0)]
    if len(values) == 0:
        return 1.0

    max_value = np.nanmax(values)
    robust_value = np.nanpercentile(values, percentile)
    if not np.isfinite(max_value) or max_value <= 0.0:
        return 1.0
    if not np.isfinite(robust_value) or robust_value <= 0.0:
        return max_value * 1.15

    if max_value > outlier_factor * robust_value:
        return robust_value * 1.35
    return max_value * 1.15


def _draw_era_fill_annotations(axes, era_ranges):
    for era_idx, (era, start_fill, end_fill) in enumerate(era_ranges):
        label_candidates = []
        for axis in axes:
            x_min, x_max = axis.get_xlim()
            overlap_start = max(start_fill, x_min)
            overlap_end = min(end_fill, x_max)
            if overlap_start > overlap_end:
                continue
            label_candidates.append((overlap_end - overlap_start, axis, overlap_start, overlap_end))

            if era_idx < len(era_ranges) - 1 and x_min <= end_fill <= x_max:
                axis.axvline(
                    end_fill,
                    linestyle="--",
                    linewidth=1.0,
                    color="black",
                    alpha=0.8,
                    zorder=1,
                )

        if not label_candidates:
            continue

        _, label_axis, overlap_start, overlap_end = max(
            label_candidates,
            key=lambda item: item[0],
        )
        label_axis.text(
            0.5 * (overlap_start + overlap_end),
            0.98,
            era,
            transform=label_axis.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=12,
            color="black",
            weight="bold",
            rotation=90,
            clip_on=False,
        )


def plot_fill_monitoring(
    trigger_dict,
    rate_cache,
    runs,
    fills,
    eras,
    eras_list,
    out_dir,
    cases=None,
    min_lumisections=10,
):
    runs = np.asarray(runs, dtype=int)
    fills = np.asarray(fills, dtype=int)
    cases = cases or []
    era_mask, selected_eras = _special_fill_era_mask(cases, fills, runs, eras, eras_list)
    if not np.any(era_mask):
        print("\033[91m[WARNING] No loaded lumisections found for configured special fills.\033[0m")
        return None

    era_text = ", ".join(selected_eras)
    print(f"\033[91m[INFO] Fill monitoring restricted to eras:\033[0m {era_text}")

    runs = runs[era_mask]
    fills = fills[era_mask]
    era_ranges = _era_fill_ranges(runs, fills, eras, selected_eras)
    output_path = Path(out_dir) / "NPSTriggerMonitoring_fill_AllCombined.pdf"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    colors = plt.get_cmap("tab10").colors
    markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "*")
    pages_written = 0

    with PdfPages(output_path) as pdf:
        for group, trigger_names in trigger_dict.items():
            page_series = []
            all_fill_numbers = []
            y_scale_values = []

            for trigger_idx, trigger_name in enumerate(trigger_names):
                rates = rate_cache.get(trigger_name)
                if rates is None:
                    continue
                rates = np.asarray(rates, dtype=float)[era_mask]

                fill_numbers, medians, p10_values, p90_values = _fill_statistics(
                    rates,
                    fills,
                    min_lumisections,
                )
                if len(fill_numbers) == 0:
                    continue

                color = colors[trigger_idx % len(colors)]
                marker = markers[trigger_idx % len(markers)]
                yerr = np.vstack([
                    medians - p10_values,
                    p90_values - medians,
                ])
                page_series.append({
                    "trigger_name": trigger_name,
                    "fill_numbers": fill_numbers,
                    "medians": medians,
                    "yerr": yerr,
                    "color": color,
                    "marker": marker,
                })
                all_fill_numbers.append(fill_numbers)
                y_scale_values.extend(p90_values)

            if not page_series:
                continue

            all_fill_numbers = np.concatenate(all_fill_numbers)
            xlims = _automatic_fill_xlims(all_fill_numbers)
            if not xlims:
                continue

            trigger_rows = len(page_series) + 1
            case_rows = len(cases) + 1 if cases else 0
            legend_gap = 1.25
            legend_area_height = max(2.4, 0.55 * max(trigger_rows, case_rows) + 0.8)
            plot_height = 8.3
            top_margin = 0.35
            total_height = plot_height + legend_area_height + legend_gap + top_margin
            fig = plt.figure(figsize=(16.0, total_height))
            fig.subplots_adjust(
                left=0.11,
                right=0.98,
                bottom=(legend_area_height + legend_gap) / total_height,
                top=(legend_area_height + legend_gap + plot_height) / total_height,
            )
            bax = brokenaxes(
                xlims=tuple(xlims),
                fig=fig,
                despine=False,
                d=0.007 if len(xlims) > 1 else 0.0,
                diag_color="black",
            )

            visible_case_labels = set()
            for axis in bax.axs:
                for visible_case in _draw_special_fill_regions(axis, cases):
                    visible_case_labels.add(visible_case["label"])
            visible_cases = [
                case for case in cases
                if case["label"] in visible_case_labels
            ]

            legend_handles = []
            legend_labels = []
            for series in page_series:
                errorbar_result = bax.errorbar(
                    series["fill_numbers"],
                    series["medians"],
                    yerr=series["yerr"],
                    fmt=series["marker"],
                    linestyle="none",
                    markersize=5.0,
                    color=series["color"],
                    ecolor=series["color"],
                    elinewidth=0.8,
                    capsize=2.0,
                    alpha=0.9,
                    label=series["trigger_name"],
                    zorder=3,
                )
                legend_handles.append(errorbar_result[0])
                legend_labels.append(series["trigger_name"])

            ymax = _overview_ymax(y_scale_values)
            bax.set_ylim(0.0, ymax)
            bax.set_xlabel("Fill Number", labelpad=50, fontsize=18, fontweight="bold")
            bax.set_ylabel(
                r"Median rate at $2.0\times10^{34}$ [Hz]",
                labelpad=45,
                fontsize=18,
                fontweight="bold",
            )

            for axis in bax.axs:
                is_last_segment = axis.get_subplotspec().is_last_col()
                axis.grid(axis="both", linestyle="--", linewidth=0.6, alpha=0.35)
                axis.tick_params(
                    axis="both",
                    which="both",
                    labelsize=14,
                    top=True,
                    right=is_last_segment,
                    labeltop=False,
                    labelright=False,
                    direction="in",
                )
                axis.spines["top"].set_visible(True)
                axis.spines["right"].set_visible(is_last_segment)
                axis.xaxis.set_major_locator(plt.MaxNLocator(nbins=4, integer=True))
                for tick in axis.get_xticklabels():
                    tick.set_rotation(45)
                    tick.set_ha("right")
                    tick.set_rotation_mode("anchor")
                    tick.set_fontweight("bold")
                for tick in axis.get_yticklabels():
                    tick.set_fontweight("bold")

            _draw_era_fill_annotations(bax.axs, era_ranges)

            legend_top = (legend_area_height - 0.15) / total_height
            extra_artists = []
            trigger_legend = fig.legend(
                handles=legend_handles,
                labels=legend_labels,
                frameon=False,
                title=f"{group} Trigger Paths",
                bbox_to_anchor=(0.11, legend_top),
                bbox_transform=fig.transFigure,
                loc="upper left",
                fontsize=13,
                title_fontsize=16,
            )
            trigger_legend.get_title().set_fontweight("bold")
            extra_artists.append(trigger_legend)

            if visible_cases:
                case_legend = fig.legend(
                    handles=[_case_patch(case) for case in visible_cases],
                    frameon=False,
                    title="Special fills",
                    bbox_to_anchor=(0.98, legend_top),
                    bbox_transform=fig.transFigure,
                    loc="upper right",
                    fontsize=12,
                    title_fontsize=15,
                )
                case_legend.get_title().set_fontweight("bold")
                extra_artists.append(case_legend)

            pdf.savefig(fig, bbox_inches="tight", bbox_extra_artists=extra_artists)
            plt.close(fig)
            pages_written += 1

    if pages_written == 0:
        output_path.unlink(missing_ok=True)
        print("\033[91m[WARNING] No valid pages for fill monitoring.\033[0m")
        return None
    return output_path


def _case_reference_mask(case, cases, fills, runs, eras, eras_list):
    fills = np.asarray(fills, dtype=int)
    runs = np.asarray(runs, dtype=int)
    if case["reference_fills"]:
        return np.isin(fills, list(case["reference_fills"]))

    if case["reference_mode"] != "sameEraExcludingSpecialFills":
        return np.zeros(len(fills), dtype=bool)

    special_mask = np.isin(fills, list(case["fills"]))
    special_runs = runs[special_mask]
    if len(special_runs) == 0:
        return np.zeros(len(fills), dtype=bool)

    relevant_era_mask = np.zeros(len(fills), dtype=bool)
    for era in eras_list:
        start_run, end_run = map(int, eras[era])
        if np.any((special_runs >= start_run) & (special_runs <= end_run)):
            relevant_era_mask |= (runs >= start_run) & (runs <= end_run)

    all_special_fills = set()
    for configured_case in cases:
        all_special_fills.update(configured_case["fills"])
    return relevant_era_mask & ~np.isin(fills, list(all_special_fills))


def plot_special_fill_rate_ratios(
    trigger_dict,
    rate_cache,
    runs,
    fills,
    eras,
    eras_list,
    out_dir,
    cases,
    min_lumisections=10,
):
    rows = _ordered_unique_triggers(trigger_dict)
    if not rows or not cases:
        print("\033[91m[WARNING] Not enough triggers/cases for special-fill ratio heatmap.\033[0m")
        return None

    runs = np.asarray(runs, dtype=int)
    fills = np.asarray(fills, dtype=int)
    values = np.full((len(rows), len(cases)), np.nan, dtype=float)

    for case_idx, case in enumerate(cases):
        special_mask = np.isin(fills, list(case["fills"]))
        reference_mask = _case_reference_mask(
            case,
            cases,
            fills,
            runs,
            eras,
            eras_list,
        )

        for row_idx, (_, trigger_name) in enumerate(rows):
            rates = rate_cache.get(trigger_name)
            if rates is None:
                continue
            rates = np.asarray(rates, dtype=float)
            special_values = rates[special_mask]
            reference_values = rates[reference_mask]
            special_values = special_values[np.isfinite(special_values)]
            reference_values = reference_values[np.isfinite(reference_values)]
            if (
                len(special_values) < min_lumisections or
                len(reference_values) < min_lumisections
            ):
                continue
            reference_median = np.nanmedian(reference_values)
            if not np.isfinite(reference_median) or reference_median <= 0.0:
                continue
            values[row_idx, case_idx] = np.nanmedian(special_values) / reference_median

    log2_values = np.full_like(values, np.nan, dtype=float)
    valid = np.isfinite(values) & (values > 0.0)
    log2_values[valid] = np.log2(values[valid])
    finite_log_values = log2_values[np.isfinite(log2_values)]
    if len(finite_log_values) == 0:
        print("\033[91m[WARNING] No valid entries for special-fill ratio heatmap.\033[0m")
        return None

    color_limit = max(1.0, np.nanpercentile(np.abs(finite_log_values), 95))
    fig_width = max(11.0, 5.0 + 2.2 * len(cases))
    fig_height = max(8.0, 2.5 + 0.35 * len(rows))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#eeeeee")
    image = ax.pcolormesh(
        np.arange(len(cases) + 1) - 0.5,
        np.arange(len(rows) + 1) - 0.5,
        np.ma.masked_invalid(log2_values),
        cmap=cmap,
        norm=TwoSlopeNorm(
            vmin=-color_limit,
            vcenter=0.0,
            vmax=color_limit,
        ),
        shading="flat",
        edgecolors="white",
        linewidth=0.4,
    )
    ax.set_xlim(-0.5, len(cases) - 0.5)
    ax.set_ylim(len(rows) - 0.5, -0.5)

    ax.set_xticks(np.arange(len(cases)))
    ax.set_xticklabels(
        [case["label"] for case in cases],
        rotation=35,
        ha="right",
        fontsize=11,
        fontweight="bold",
    )
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(
        [trigger for _, trigger in rows],
        fontsize=10,
        fontweight="bold",
    )
    ax.tick_params(axis="both", which="both", length=0)

    previous_group = rows[0][0]
    for row_idx, (group, _) in enumerate(rows[1:], start=1):
        if group == previous_group:
            continue
        ax.axhline(row_idx - 0.5, color="black", linewidth=0.8)
        previous_group = group

    for row_idx in range(len(rows)):
        for case_idx in range(len(cases)):
            ratio = values[row_idx, case_idx]
            if not np.isfinite(ratio):
                continue
            log_value = log2_values[row_idx, case_idx]
            text_color = "white" if abs(log_value) > 0.65 * color_limit else "black"
            ax.text(
                case_idx,
                row_idx,
                f"{ratio:.2f}",
                ha="center",
                va="center",
                fontsize=7,
                fontweight="bold",
                color=text_color,
            )

    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label(
        "Median Rate [Special-fill / Reference]",
        labelpad=12,
    )
    ratio_ticks = np.asarray([0.5, 0.75, 1.0, 1.5, 2.0])
    tick_positions = np.log2(ratio_ticks)
    visible_ticks = np.abs(tick_positions) <= color_limit
    colorbar.set_ticks(tick_positions[visible_ticks])
    colorbar.set_ticklabels([
        f"{ratio:g}" for ratio in ratio_ticks[visible_ticks]
    ])
    fig.tight_layout()
    output_path = Path(out_dir) / "specialFillRateRatios.pdf"
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path

