from pathlib import Path

from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm

from bin.rateChangeMatrix import _build_era_groups, _era_mask, _ordered_unique_triggers
from bin.utils import *


def _year_from_era_label(era_label):
    return str(era_label)[:4]


def _era_masks(runs, eras, era_groups):
    runs = np.asarray(runs)
    good_run_mask = get_good_run_mask(runs)
    return [
        _era_mask(runs, eras, era_keys) & good_run_mask
        for _, era_keys in era_groups
    ]


def _finite_values(values):
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def _median_for_mask(rates, mask, min_lumisections):
    values = _finite_values(np.asarray(rates, dtype=float)[mask])
    if len(values) < min_lumisections:
        return np.nan
    return np.nanmedian(values)


def _ratio_text(value):
    if not np.isfinite(value):
        return ""
    if value >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def _small_ratio_text(value):
    if not np.isfinite(value):
        return ""
    if value >= 0.1:
        return f"{value:.2f}"
    if value >= 0.01:
        return f"{value:.3f}"
    if value >= 0.001:
        return f"{value:.4f}"
    return f"{value:.1e}"


def _year_color(label):
    year = _year_from_era_label(label)
    return {
        "2024": "#4C78A8",
        "2025": "#F58518",
        "2026": "#54A24B",
    }.get(year, "#8C8C8C")


def _sum_trigger_rates(trigger_dict, rate_cache, n_entries):
    total_rate = np.zeros(n_entries, dtype=float)
    has_value = np.zeros(n_entries, dtype=bool)

    for _, trigger_name in _ordered_unique_triggers(trigger_dict):
        rates = rate_cache.get(trigger_name)
        if rates is None:
            continue

        rates = np.asarray(rates, dtype=float)
        if len(rates) != n_entries:
            print(
                f"\033[91m[WARNING] Skipping {trigger_name} in NPS rate box plot "
                f"because its rate array length does not match the selected LS.\033[0m"
            )
            continue

        finite = np.isfinite(rates)
        total_rate[finite] += rates[finite]
        has_value |= finite

    total_rate[~has_value] = np.nan
    return total_rate


def _box_data_for_masks(values, masks, min_lumisections, positive_only=False):
    box_data = []
    positions = []
    for idx, mask in enumerate(masks):
        era_values = _finite_values(values[mask])
        if positive_only:
            era_values = era_values[era_values > 0.0]
        if len(era_values) < min_lumisections:
            continue
        box_data.append(era_values)
        positions.append(idx)
    return box_data, positions


def _draw_year_separators(ax, labels):
    for idx in range(1, len(labels)):
        if _year_from_era_label(labels[idx]) != _year_from_era_label(labels[idx - 1]):
            ax.axvline(idx - 0.5, color="black", linestyle="--", linewidth=0.8, alpha=0.7)


def _draw_box_panel(ax, box_data, positions, labels, ylabel, log_scale=False):
    if not box_data:
        return

    box = ax.boxplot(
        box_data,
        positions=positions,
        widths=0.65,
        whis=1.5,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
        medianprops={"color": "black", "linewidth": 2.2},
        whiskerprops={"color": "#333333", "linewidth": 1.2},
        capprops={"color": "#333333", "linewidth": 1.2},
    )

    for patch, pos in zip(box["boxes"], positions):
        color = _year_color(labels[pos])
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        patch.set_edgecolor("#222222")
        patch.set_linewidth(1.0)

    ax.set_ylabel(ylabel, fontsize=15, fontweight="bold")
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.35)
    ax.tick_params(
        axis="both",
        which="both",
        labelsize=11,
        top=True,
        right=True,
        labeltop=False,
        labelright=False,
        direction="in",
    )
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.margins(y=0.12)
    if log_scale:
        ax.set_yscale("log")
    _draw_year_separators(ax, labels)


def _group_color_map(trigger_dict):
    colors = plt.get_cmap("tab10").colors
    return {
        group: colors[idx % len(colors)]
        for idx, group in enumerate(trigger_dict)
    }


def _short_trigger_name(trigger_name):
    return re.sub(r"^HLT_", "", str(trigger_name))


def _trigger_rows(trigger_dict):
    return _ordered_unique_triggers(trigger_dict)


def _draw_group_separators(ax, rows):
    previous_group = rows[0][0] if rows else None
    for row_idx, (group, _) in enumerate(rows[1:], start=1):
        if group == previous_group:
            continue
        ax.axvline(row_idx - 0.5, color="#555555", linewidth=0.9, alpha=0.8)
        previous_group = group


def _binned_median(x_values, y_values, bin_edges, min_lumisections):
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    finite = np.isfinite(x_values) & np.isfinite(y_values)
    x_values = x_values[finite]
    y_values = y_values[finite]

    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    medians = np.full(len(centers), np.nan, dtype=float)
    for idx in range(len(centers)):
        if idx == len(centers) - 1:
            in_bin = (x_values >= bin_edges[idx]) & (x_values <= bin_edges[idx + 1])
        else:
            in_bin = (x_values >= bin_edges[idx]) & (x_values < bin_edges[idx + 1])
        if np.count_nonzero(in_bin) < min_lumisections:
            continue
        medians[idx] = np.nanmedian(y_values[in_bin])
    return centers, medians


def plot_trigger_rate_boxplots_by_era(
    trigger_dict,
    rate_cache,
    runs,
    eras,
    eras_list,
    out_dir,
    merge_era_versions=True,
    min_lumisections=10,
):
    rows = _trigger_rows(trigger_dict)
    era_groups = _build_era_groups(eras_list, merge_era_versions=merge_era_versions)
    if not rows or not era_groups:
        print("\033[91m[WARNING] Not enough eras/triggers for per-era trigger box plots.\033[0m")
        return None

    masks = _era_masks(runs, eras, era_groups)
    group_colors = _group_color_map(trigger_dict)
    output_path = Path(out_dir) / "NPSTriggerRateBoxPlots_byEra_AllCombined.pdf"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pages_written = 0

    with PdfPages(output_path) as pdf:
        for (era_label, _), era_mask in zip(era_groups, masks):
            box_data = []
            positions = []
            box_groups = []
            for row_idx, (group, trigger_name) in enumerate(rows):
                rates = rate_cache.get(trigger_name)
                if rates is None:
                    continue
                values = _finite_values(np.asarray(rates, dtype=float)[era_mask])
                if len(values) < min_lumisections:
                    continue
                box_data.append(values)
                positions.append(row_idx)
                box_groups.append(group)

            if not box_data:
                continue

            fig_width = max(18.0, 0.62 * len(rows))
            group_labels = list(trigger_dict.keys())
            legend_columns = min(4, max(1, len(group_labels)))
            legend_rows = int(np.ceil(len(group_labels) / legend_columns))
            legend_area_height = max(1.6, 0.58 * legend_rows + 0.75)
            tick_label_area_height = 4.0
            plot_height = 8.5
            top_margin = 0.45
            total_height = (
                plot_height +
                tick_label_area_height +
                legend_area_height +
                top_margin
            )
            fig, ax = plt.subplots(figsize=(fig_width, total_height))
            fig.subplots_adjust(
                left=0.08,
                right=0.98,
                bottom=(legend_area_height + tick_label_area_height) / total_height,
                top=(legend_area_height + tick_label_area_height + plot_height) / total_height,
            )
            box = ax.boxplot(
                box_data,
                positions=positions,
                widths=0.62,
                whis=1.5,
                patch_artist=True,
                showfliers=False,
                manage_ticks=False,
                medianprops={"color": "black", "linewidth": 2.0},
                whiskerprops={"color": "#333333", "linewidth": 1.1},
                capprops={"color": "#333333", "linewidth": 1.1},
            )

            for patch, group in zip(box["boxes"], box_groups):
                patch.set_facecolor(group_colors[group])
                patch.set_alpha(0.78)
                patch.set_edgecolor("#222222")
                patch.set_linewidth(0.9)

            ax.set_title(era_label, fontsize=20, fontweight="bold", pad=12)
            ax.set_ylabel(
                r"HLT rate at $2.0\times10^{34}\,\mathrm{cm}^{-2}\mathrm{s}^{-1}$ [Hz]",
                fontsize=15,
                fontweight="bold",
            )
            ax.set_xlim(-0.75, len(rows) - 0.25)
            ax.set_xticks(np.arange(len(rows)))
            ax.set_xticklabels(
                [_short_trigger_name(trigger_name) for _, trigger_name in rows],
                rotation=55,
                ha="right",
                fontsize=9,
                fontweight="bold",
            )
            ax.set_yscale("symlog", linthresh=0.1, linscale=0.6)
            ax.grid(axis="y", which="both", linestyle="--", linewidth=0.55, alpha=0.35)
            ax.tick_params(
                axis="both",
                which="both",
                labelsize=11,
                top=True,
                right=True,
                labeltop=False,
                labelright=False,
                direction="in",
            )
            ax.spines["top"].set_visible(True)
            ax.spines["right"].set_visible(True)
            _draw_group_separators(ax, rows)

            legend_handles = [
                Rectangle((0, 0), 1, 1, facecolor=group_colors[group], alpha=0.78, edgecolor="#222222")
                for group in group_labels
            ]
            legend = fig.legend(
                legend_handles,
                group_labels,
                frameon=False,
                title="Trigger Groups",
                bbox_to_anchor=(0.08, (legend_area_height - 0.15) / total_height),
                bbox_transform=fig.transFigure,
                loc="upper left",
                fontsize=14,
                title_fontsize=16,
                ncol=legend_columns,
                columnspacing=1.9,
                handlelength=1.5,
            )
            legend.get_title().set_fontweight("bold")

            pdf.savefig(fig, bbox_inches="tight", bbox_extra_artists=[legend])
            plt.close(fig)
            pages_written += 1

    if pages_written == 0:
        output_path.unlink(missing_ok=True)
        print("\033[91m[WARNING] No valid pages for per-era trigger box plots.\033[0m")
        return None
    return output_path


def plot_rate_vs_pileup_by_era(
    trigger_dict,
    rate_cache,
    runs,
    pileup,
    eras,
    eras_list,
    out_dir,
    merge_era_versions=True,
    pileup_bin_width=0.5,
    min_lumisections_per_bin=3,
):
    era_groups = _build_era_groups(eras_list, merge_era_versions=merge_era_versions)
    if not trigger_dict or not era_groups:
        print("\033[91m[WARNING] Not enough eras/triggers for rate-vs-pileup plots.\033[0m")
        return None

    pileup = np.asarray(pileup, dtype=float)
    masks = _era_masks(runs, eras, era_groups)
    output_dir = Path(out_dir) / "RateVsPileupByEra"
    output_dir.mkdir(parents=True, exist_ok=True)
    legacy_output = Path(out_dir) / "NPSTriggerRateVsPileup_byEra_AllCombined.pdf"
    if legacy_output.exists():
        legacy_output.unlink()

    output_paths = []
    colors = plt.get_cmap("tab10").colors
    markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "*")

    for (era_label, _), era_mask in zip(era_groups, masks):
        era_pileup = pileup[era_mask]
        finite_pileup = era_pileup[np.isfinite(era_pileup)]
        if len(finite_pileup) < min_lumisections_per_bin:
            continue

        pu_min = 0.0
        pu_max = max(
            75.0,
            np.ceil(np.nanmax(finite_pileup) / pileup_bin_width) * pileup_bin_width,
        )
        if pu_max <= pu_min:
            pu_max = pu_min + pileup_bin_width
        bin_edges = np.arange(pu_min, pu_max + pileup_bin_width, pileup_bin_width)

        safe_era_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", era_label)
        output_path = output_dir / f"NPSTriggerRateVsPileup_{safe_era_label}.pdf"
        pages_written = 0

        with PdfPages(output_path) as pdf:
            for group, trigger_names in trigger_dict.items():
                fig, ax = plt.subplots(figsize=(16.0, 10.0))
                plotted_triggers = []
                for trigger_idx, trigger_name in enumerate(trigger_names):
                    rates = rate_cache.get(trigger_name)
                    if rates is None:
                        continue
                    rates = np.asarray(rates, dtype=float)
                    centers, medians = _binned_median(
                        era_pileup,
                        rates[era_mask],
                        bin_edges,
                        min_lumisections_per_bin,
                    )
                    valid = np.isfinite(medians)
                    if np.count_nonzero(valid) < 2:
                        continue
                    ax.scatter(
                        centers[valid],
                        medians[valid],
                        marker=markers[trigger_idx % len(markers)],
                        s=42,
                        color=colors[trigger_idx % len(colors)],
                        edgecolors="none",
                        label=trigger_name,
                    )
                    plotted_triggers.append(trigger_name)

                if not plotted_triggers:
                    plt.close(fig)
                    continue

                ax.set_title(era_label, fontsize=22, fontweight="bold", pad=14)
                ax.set_xlabel(r"$\langle \mathrm{PU} \rangle$", fontsize=18, fontweight="bold")
                ax.set_ylabel(
                    r"Median rate at $2.0\times10^{34}$ [Hz]",
                    fontsize=18,
                    fontweight="bold",
                )
                ax.set_xlim(0.0, pu_max)
                ax.grid(axis="both", linestyle="--", linewidth=0.6, alpha=0.35)
                ax.tick_params(axis="both", labelsize=14, top=False, right=False)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                for tick in ax.get_xticklabels() + ax.get_yticklabels():
                    tick.set_fontweight("bold")

                legend = ax.legend(
                    frameon=False,
                    title=f"{group} Trigger Paths",
                    bbox_to_anchor=(0.0, -0.16),
                    loc="upper left",
                    fontsize=13,
                    title_fontsize=16,
                    borderaxespad=0.0,
                )
                legend.get_title().set_fontweight("bold")

                fig.subplots_adjust(left=0.11, right=0.98, top=0.91, bottom=0.12)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
                pages_written += 1

        if pages_written == 0:
            output_path.unlink(missing_ok=True)
            continue
        output_paths.append(output_path)

    if not output_paths:
        print("\033[91m[WARNING] No valid per-era rate-vs-pileup PDFs were produced.\033[0m")
        return None
    return output_paths


def _draw_heatmap(
    values,
    row_groups,
    row_labels,
    column_labels,
    out_path,
    colorbar_label,
    cmap,
    norm,
    text_formatter,
    text_color_fn=None,
):
    output_path = Path(out_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig_width = max(14.0, 7.0 + 0.55 * len(column_labels))
    fig_height = max(8.0, 2.5 + 0.35 * len(row_labels))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    cmap = cmap.copy()
    cmap.set_bad("#eeeeee")
    image = ax.imshow(values, aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(np.arange(len(column_labels)))
    ax.set_xticklabels(column_labels, rotation=45, ha="right", fontsize=11, fontweight="bold")
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=10, fontweight="bold")

    for row_idx in range(len(row_labels)):
        for col_idx in range(len(column_labels)):
            value = values[row_idx, col_idx]
            text = text_formatter(value)
            if not text:
                continue
            color = text_color_fn(value) if text_color_fn else "black"
            ax.text(
                col_idx,
                row_idx,
                text,
                ha="center",
                va="center",
                fontsize=6,
                fontweight="bold",
                color=color,
            )

    previous_group = row_groups[0] if row_groups else None
    for row_idx, group in enumerate(row_groups[1:], start=1):
        if group == previous_group:
            continue
        ax.axhline(row_idx - 0.5, color="black", linewidth=0.8)
        previous_group = group

    ax.set_xticks(np.arange(-0.5, len(column_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.4)
    ax.tick_params(axis="both", which="both", length=0)

    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label(colorbar_label, labelpad=14)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_nps_rate_boxplot_by_era(
    trigger_dict,
    rate_cache,
    runs,
    eras,
    eras_list,
    out_dir,
    merge_era_versions=True,
    total_rate=None,
    total_rate_label=None,
    min_lumisections=10,
):
    era_groups = _build_era_groups(eras_list, merge_era_versions=merge_era_versions)
    if len(era_groups) < 1:
        print("\033[91m[WARNING] Not enough eras for NPS rate box plot.\033[0m")
        return None

    runs = np.asarray(runs)
    masks = _era_masks(runs, eras, era_groups)
    labels = [label for label, _ in era_groups]
    total_nps_rate = _sum_trigger_rates(trigger_dict, rate_cache, len(runs))
    rate_box_data, rate_positions = _box_data_for_masks(
        total_nps_rate,
        masks,
        min_lumisections,
    )

    if not rate_box_data:
        print("\033[91m[WARNING] No valid entries for NPS rate box plot.\033[0m")
        return None

    fraction_box_data = []
    fraction_positions = []
    fraction_values = None
    if total_rate is not None:
        total_rate = np.asarray(total_rate, dtype=float)
        if len(total_rate) != len(total_nps_rate):
            print("\033[91m[WARNING] Total-rate denominator length does not match selected LS. Skipping percentage panel.\033[0m")
        else:
            fraction_values = np.full(len(total_nps_rate), np.nan, dtype=float)
            valid = (
                np.isfinite(total_nps_rate) &
                np.isfinite(total_rate) &
                (total_rate > 0.0)
            )
            fraction_values[valid] = 100.0 * total_nps_rate[valid] / total_rate[valid]
            fraction_box_data, fraction_positions = _box_data_for_masks(
                fraction_values,
                masks,
                min_lumisections,
                positive_only=True,
            )
            if not fraction_box_data:
                print("\033[91m[WARNING] No valid NPS/total-rate percentage entries. Skipping percentage panel.\033[0m")

    include_fraction = bool(fraction_box_data)
    n_panels = 2 if include_fraction else 1
    fig_width = max(14.0, 0.72 * len(labels))
    fig_height = 8.5 if include_fraction else 6.2
    fig, axes = plt.subplots(
        n_panels,
        1,
        figsize=(fig_width, fig_height),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]} if include_fraction else None,
    )
    if not include_fraction:
        axes = [axes]

    _draw_box_panel(
        axes[0],
        rate_box_data,
        rate_positions,
        labels,
        r"Total NPS HLT rate at $2.0e34$ [Hz]",
    )

    if include_fraction:
        ylabel = "NPS / total rate [%]"
        if total_rate_label:
            ylabel = f"NPS / {total_rate_label} [%]"
        _draw_box_panel(
            axes[1],
            fraction_box_data,
            fraction_positions,
            labels,
            ylabel,
            log_scale=True,
        )

    axes[-1].set_xlim(-0.5, len(labels) - 0.5)
    axes[-1].set_xticks(np.arange(len(labels)))
    axes[-1].set_xticklabels(labels, rotation=40, ha="right", fontsize=11, fontweight="bold")

    for ax in axes:
        for tick in ax.get_yticklabels():
            tick.set_fontweight("bold")

    fig.tight_layout()
    output_path = Path(out_dir) / "npsRateBoxPlot_byEra.pdf"
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_trigger_stability_heatmap(
    trigger_dict,
    rate_cache,
    runs,
    eras,
    eras_list,
    out_dir,
    merge_era_versions=True,
    min_lumisections=10,
):
    rows = _ordered_unique_triggers(trigger_dict)
    era_groups = _build_era_groups(eras_list, merge_era_versions=merge_era_versions)
    if len(era_groups) < 1 or not rows:
        print("\033[91m[WARNING] Not enough eras/triggers for trigger-stability heatmap.\033[0m")
        return None

    masks = _era_masks(runs, eras, era_groups)
    values = np.full((len(rows), len(era_groups)), np.nan, dtype=float)

    for row_idx, (_, trigger_name) in enumerate(rows):
        rates = rate_cache.get(trigger_name)
        if rates is None:
            continue
        rates = np.asarray(rates, dtype=float)
        for era_idx, mask in enumerate(masks):
            finite = _finite_values(rates[mask])
            if len(finite) < min_lumisections:
                continue
            median = np.nanmedian(finite)
            if not np.isfinite(median) or median <= 0.0:
                continue
            p10, p90 = np.nanpercentile(finite, [10, 90])
            values[row_idx, era_idx] = (p90 - p10) / median

    finite_values = values[np.isfinite(values)]
    if len(finite_values) == 0:
        print("\033[91m[WARNING] No valid entries for trigger-stability heatmap.\033[0m")
        return None

    vmax = max(0.5, np.nanpercentile(finite_values, 95))
    return _draw_heatmap(
        values,
        [group for group, _ in rows],
        [trigger for _, trigger in rows],
        [label for label, _ in era_groups],
        Path(out_dir) / "triggerStabilityHeatmap.pdf",
        r"$(P90 - P10) / median$",
        plt.get_cmap("YlOrRd"),
        Normalize(vmin=0.0, vmax=vmax),
        _ratio_text,
        text_color_fn=lambda value: "white" if value > 0.7 * vmax else "black",
    )


def plot_year_normalized_era_heatmap(
    trigger_dict,
    rate_cache,
    runs,
    eras,
    eras_list,
    out_dir,
    merge_era_versions=True,
    min_lumisections=10,
):
    rows = _ordered_unique_triggers(trigger_dict)
    era_groups = _build_era_groups(eras_list, merge_era_versions=merge_era_versions)
    if len(era_groups) < 1 or not rows:
        print("\033[91m[WARNING] Not enough eras/triggers for year-normalized heatmap.\033[0m")
        return None

    masks = _era_masks(runs, eras, era_groups)
    era_labels = [label for label, _ in era_groups]
    era_years = [_year_from_era_label(label) for label in era_labels]
    years = list(dict.fromkeys(era_years))
    year_masks = {
        year: np.logical_or.reduce([
            mask for mask, era_year in zip(masks, era_years) if era_year == year
        ])
        for year in years
    }

    values = np.full((len(rows), len(era_groups)), np.nan, dtype=float)
    for row_idx, (_, trigger_name) in enumerate(rows):
        rates = rate_cache.get(trigger_name)
        if rates is None:
            continue
        rates = np.asarray(rates, dtype=float)
        year_medians = {
            year: _median_for_mask(rates, year_mask, min_lumisections)
            for year, year_mask in year_masks.items()
        }
        for era_idx, mask in enumerate(masks):
            year_median = year_medians[era_years[era_idx]]
            era_median = _median_for_mask(rates, mask, min_lumisections)
            if not np.isfinite(year_median) or year_median <= 0.0:
                continue
            if not np.isfinite(era_median):
                continue
            values[row_idx, era_idx] = era_median / year_median

    log2_values = np.full_like(values, np.nan, dtype=float)
    valid = np.isfinite(values) & (values > 0.0)
    log2_values[valid] = np.log2(values[valid])

    return _draw_heatmap(
        log2_values,
        [group for group, _ in rows],
        [trigger for _, trigger in rows],
        era_labels,
        Path(out_dir) / "yearNormalizedEraHeatmap.pdf",
        "Era median / same-year median",
        plt.get_cmap("RdBu_r"),
        TwoSlopeNorm(vmin=-2.0, vcenter=0.0, vmax=2.0),
        lambda value: _ratio_text(2.0 ** value) if np.isfinite(value) else "",
        text_color_fn=lambda value: "white" if np.isfinite(value) and abs(value) > 1.25 else "black",
    )


def plot_hlt_dominant_l1_seed_heatmap(
    hlt_trigger_dict,
    hlt_rate_cache,
    l1_seed_dict,
    l1_rate_cache,
    runs,
    eras,
    eras_list,
    out_dir,
    merge_era_versions=True,
    min_lumisections=10,
):
    rows = _ordered_unique_triggers(hlt_trigger_dict)
    era_groups = _build_era_groups(eras_list, merge_era_versions=merge_era_versions)
    if len(era_groups) < 1 or not rows:
        print("\033[91m[WARNING] Not enough eras/triggers for HLT/L1 heatmap.\033[0m")
        return None

    masks = _era_masks(runs, eras, era_groups)
    values = np.full((len(rows), len(era_groups)), np.nan, dtype=float)

    for row_idx, (_, hlt_path) in enumerate(rows):
        hlt_rates = hlt_rate_cache.get(hlt_path)
        if hlt_rates is None:
            continue
        hlt_rates = np.asarray(hlt_rates, dtype=float)
        l1_seeds = l1_seed_dict.get(hlt_path, [])
        if not l1_seeds:
            continue

        for era_idx, mask in enumerate(masks):
            hlt_median = _median_for_mask(hlt_rates, mask, min_lumisections)
            if not np.isfinite(hlt_median):
                continue

            l1_medians = []
            for l1_seed in l1_seeds:
                l1_rates = l1_rate_cache.get(l1_seed)
                if l1_rates is None:
                    continue
                l1_median = _median_for_mask(l1_rates, mask, min_lumisections)
                if np.isfinite(l1_median):
                    l1_medians.append(l1_median)

            if not l1_medians:
                continue
            dominant_l1_median = np.nanmax(l1_medians)
            if dominant_l1_median <= 0.0:
                continue
            values[row_idx, era_idx] = hlt_median / dominant_l1_median

    finite_values = values[np.isfinite(values) & (values > 0.0)]
    if len(finite_values) == 0:
        print("\033[91m[WARNING] No valid entries for HLT/dominant-L1 heatmap.\033[0m")
        return None

    log_min = np.floor(np.log10(np.nanmin(finite_values)))
    log_max = np.ceil(np.log10(np.nanmax(finite_values)))
    if log_min == log_max:
        log_min -= 1
        log_max += 1
    norm = LogNorm(vmin=10 ** log_min, vmax=10 ** log_max)

    def text_color(value):
        if not np.isfinite(value) or value <= 0:
            return "black"
        return "white" if norm(value) < 0.35 else "black"

    return _draw_heatmap(
        values,
        [group for group, _ in rows],
        [trigger for _, trigger in rows],
        [label for label, _ in era_groups],
        Path(out_dir) / "hltDominantL1SeedHeatmap.pdf",
        "Median HLT rate / median dominant L1 seed rate",
        plt.get_cmap("viridis"),
        norm,
        _small_ratio_text,
        text_color_fn=text_color,
    )

