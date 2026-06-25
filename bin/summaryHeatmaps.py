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

