from bin.utils import *
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path


def _base_era_name(era_name):
    return re.sub(r"v\d+$", "", era_name)


def _build_era_groups(eras_list, merge_era_versions=True):
    if not merge_era_versions:
        return [(era, [era]) for era in eras_list]

    era_groups = []
    for era in eras_list:
        label = _base_era_name(era)
        if era_groups and era_groups[-1][0] == label:
            era_groups[-1][1].append(era)
        else:
            era_groups.append((label, [era]))
    return era_groups


def _ordered_unique_triggers(trigger_dict):
    rows = []
    seen = set()
    for group, triggers in trigger_dict.items():
        for trigger_name in triggers:
            if trigger_name in seen:
                continue
            rows.append((group, trigger_name))
            seen.add(trigger_name)
    return rows


def _era_mask(runs, eras, era_keys):
    mask = np.zeros(len(runs), dtype=bool)
    for era_key in era_keys:
        start_run, end_run = int(eras[era_key][0]), int(eras[era_key][1])
        mask |= (runs >= start_run) & (runs <= end_run)
    return mask


def _ratio_text(value):
    if not np.isfinite(value):
        return ""
    if value >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def plot_rate_change_matrix(
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

    if len(era_groups) < 2 or not rows:
        print("\033[91m[WARNING] Not enough eras/triggers for rate-change matrix.\033[0m")
        return None

    runs = np.asarray(runs)
    good_run_mask = get_good_run_mask(runs)
    group_masks = []
    group_lumisections = []

    for _, era_keys in era_groups:
        mask = _era_mask(runs, eras, era_keys) & good_run_mask
        group_masks.append(mask)
        group_lumisections.append(int(np.count_nonzero(mask)))

    era_labels = [label for label, _ in era_groups]
    transition_labels = [
        f"{era_labels[idx]} / {era_labels[idx - 1]}"
        for idx in range(1, len(era_labels))
    ]

    median_rates = np.full((len(rows), len(era_groups)), np.nan, dtype=float)
    finite_lumisections = np.zeros((len(rows), len(era_groups)), dtype=int)

    for row_idx, (_, trigger_name) in enumerate(rows):
        rates = rate_cache.get(trigger_name)
        if rates is None:
            continue
        rates = np.asarray(rates, dtype=float)

        for era_idx, mask in enumerate(group_masks):
            values = rates[mask]
            finite = values[np.isfinite(values)]
            finite_lumisections[row_idx, era_idx] = len(finite)
            if len(finite) < min_lumisections:
                continue
            median_rates[row_idx, era_idx] = np.nanmedian(finite)

    ratios = np.full((len(rows), len(transition_labels)), np.nan, dtype=float)
    for era_idx in range(1, len(era_groups)):
        previous = median_rates[:, era_idx - 1]
        current = median_rates[:, era_idx]
        valid = np.isfinite(previous) & np.isfinite(current) & (previous > 0.0)
        ratios[valid, era_idx - 1] = current[valid] / previous[valid]

    log2_ratios = np.full_like(ratios, np.nan, dtype=float)
    valid_ratios = np.isfinite(ratios) & (ratios > 0.0)
    log2_ratios[valid_ratios] = np.log2(ratios[valid_ratios])

    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    row_groups = [group for group, _ in rows]
    row_triggers = [trigger for _, trigger in rows]

    fig_width = max(14.0, 7.0 + 0.55 * len(transition_labels))
    fig_height = max(8.0, 2.5 + 0.35 * len(rows))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#eeeeee")
    norm = TwoSlopeNorm(vmin=-2.0, vcenter=0.0, vmax=2.0)
    image = ax.imshow(log2_ratios, aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(np.arange(len(transition_labels)))
    ax.set_xticklabels(transition_labels, rotation=45, ha="right", fontsize=11, fontweight="bold")
    ax.set_yticks(np.arange(len(row_triggers)))
    ax.set_yticklabels(row_triggers, fontsize=10, fontweight="bold")

    for row_idx in range(len(rows)):
        for col_idx in range(len(transition_labels)):
            text = _ratio_text(ratios[row_idx, col_idx])
            if not text:
                continue
            color = "white" if abs(log2_ratios[row_idx, col_idx]) > 1.25 else "black"
            ax.text(col_idx, row_idx, text, ha="center", va="center", fontsize=6, fontweight="bold", color=color)

    previous_group = row_groups[0]
    for row_idx, group in enumerate(row_groups[1:], start=1):
        if group == previous_group:
            continue
        ax.axhline(row_idx - 0.5, color="black", linewidth=0.8)
        previous_group = group

    ax.set_xticks(np.arange(-0.5, len(transition_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.4)
    ax.tick_params(axis="both", which="both", length=0)

    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label("Median rate ratio")
    colorbar.set_ticks([-2, -1, 0, 1, 2])
    colorbar.set_ticklabels(["0.25", "0.5", "1", "2", "4"])

    fig.tight_layout()
    output_path = output_dir / "rateChangeMatrix_adjacentEraRatios.pdf"
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

    return output_path

