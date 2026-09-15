"""
Interactive Wright-Fisher Model (v5)
-------------------------------------
Uses matplotlib's built-in Slider / TextBox / Button widgets only
(no ipywidgets, no Jupyter widget extension required).

Run as a standalone script:
    python wright_fisher_interactive.py
(requires a GUI backend, e.g. TkAgg, which ships with most Python installs)

Or inside Jupyter Lab:
    %matplotlib widget      # requires: pip install ipympl
    (then paste/run this code in a cell)
Do NOT use the default inline backend - sliders/textboxes need a live canvas.

Layout:
- Top-left    : trajectory plot, width:height = 2:1 (wide), ALL paths shown
- Bottom-left : two SEPARATE square histograms side by side
                (generation of fixation / generation of loss)
- Top-right   : parameter labels (own line) + slider + numeric textbox
                (only the textbox shows the number - slider's own value
                display is hidden to avoid duplication)
- Bottom-right: results panel with bold section headers
                (Path counts / Probability of being Fixed or Lost /
                 Theoretical Fixation Probability)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------------
# Core simulation
# ---------------------------------------------------------------------------
def wright_fisher_simulate(n_alleles, p0, generations, num_paths=200):
    """
    Simulate the neutral Wright-Fisher model.

    Each generation, the number of A alleles is drawn from a
    Binomial(n_alleles, current_frequency) distribution
    (no selection, no mutation). A fresh random seed is used every call.
    """
    rng = np.random.default_rng()  # no fixed seed - fresh randomness each run
    p0 = float(np.clip(p0, 0.0, 1.0))

    freqs = np.full(num_paths, p0)
    trajectory = np.zeros((generations + 1, num_paths))
    trajectory[0] = freqs

    for t in range(1, generations + 1):
        counts = rng.binomial(n_alleles, freqs)
        freqs = counts / n_alleles
        trajectory[t] = freqs

    return trajectory


def analyze_absorption(trajectory):
    """
    For every path, find the first generation at which it becomes fixed
    (frequency = 1) or lost (frequency = 0).
    """
    boundary = (trajectory == 0.0) | (trajectory == 1.0)  # (gens+1, num_paths)
    any_hit = boundary.any(axis=0)

    first_hit = np.where(any_hit, boundary.argmax(axis=0), -1)

    safe_idx = first_hit.copy()
    safe_idx[safe_idx == -1] = 0
    hit_values = trajectory[safe_idx, np.arange(trajectory.shape[1])]

    is_fixed = any_hit & (hit_values == 1.0)
    is_lost = any_hit & (hit_values == 0.0)

    return first_hit, is_fixed, is_lost


# ---------------------------------------------------------------------------
# Initial parameter values
# ---------------------------------------------------------------------------
init_n_alleles = 100
init_generations = 100
init_num_paths = 300
init_p0 = 0.5

PARAM_BOUNDS = {
    "n_alleles": (2, 5000),
    "generations": (1, 3000),
    "num_paths": (10, 5000),
    "p0": (0.0, 1.0),
}

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(14, 10))

# --- left column ---
ax_traj = fig.add_axes([0.05, 0.56, 0.57, 0.40])
ax_traj.set_box_aspect(0.5)  # width:height = 2:1

ax_hist_fixed = fig.add_axes([0.05, 0.10, 0.27, 0.38])
ax_hist_lost = fig.add_axes([0.35, 0.10, 0.27, 0.38])
ax_hist_fixed.set_box_aspect(1)
ax_hist_lost.set_box_aspect(1)

ax_traj.set_xlabel("Generation")
ax_traj.set_ylabel("Frequency of allele A")
ax_traj.set_ylim(-0.02, 1.02)
ax_traj.grid(alpha=0.3)

ax_hist_fixed.set_xlabel("Generation of fixation")
ax_hist_fixed.set_ylabel("Count")
ax_hist_fixed.grid(alpha=0.3)

ax_hist_lost.set_xlabel("Generation of loss")
ax_hist_lost.set_ylabel("Count")
ax_hist_lost.grid(alpha=0.3)

# legend handles for fixed/lost markers (colors match the dots on the plot)
legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor="green",
           markeredgecolor="black", markersize=7, label="Fixed"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="red",
           markeredgecolor="black", markersize=7, label="Lost"),
]

# ---------------------------------------------------------------------------
# Right column, top half: parameter label (own line) + slider + textbox
# only the textbox displays the number (slider's own value text is hidden)
# ---------------------------------------------------------------------------
RX = 0.66  # right-column x anchor

param_rows = [
    ("Number of alleles (2N)", "n_alleles", init_n_alleles, True, 0.93, 0.885),
    ("Number of generations", "generations", init_generations, True, 0.83, 0.785),
    ("Number of paths", "num_paths", init_num_paths, True, 0.73, 0.685),
    ("Initial p(A)", "p0", init_p0, False, 0.63, 0.585),
]

sliders = {}
textboxes = {}

for label_txt, key, init_val, is_int, label_y, slider_y in param_rows:
    fig.text(RX, label_y, label_txt, fontsize=12, fontname="Arial", va="bottom")
    ax_s = fig.add_axes([RX, slider_y, 0.24, 0.03])
    ax_t = fig.add_axes([RX + 0.26, slider_y, 0.07, 0.03])
    lo, hi = PARAM_BOUNDS[key]
    step = 1 if is_int else 0.01
    slider = Slider(ax_s, "", lo, hi, valinit=init_val, valstep=step)
    slider.valtext.set_visible(False)  # avoid duplicating the textbox's number
    textbox = TextBox(ax_t, "", initial=str(init_val))
    sliders[key] = slider
    textboxes[key] = textbox

ax_button = fig.add_axes([RX, 0.52, 0.20, 0.045])
button = Button(ax_button, "Run Simulation")

# ---------------------------------------------------------------------------
# Right column, bottom half: results panel (bold section headers, no dashes)
# ---------------------------------------------------------------------------
title_counts = fig.text(RX, 0.46, "Path counts", fontsize=14,
                         fontweight="bold", fontname="Arial", va="top")
body_counts = fig.text(RX, 0.415, "", fontsize=11, fontname="Arial", va="top")

title_prob = fig.text(RX, 0.28, "Probability of being Fixed or Lost", fontsize=14,
                       fontweight="bold", fontname="Arial", va="top")
body_prob = fig.text(RX, 0.225, "", fontsize=11, fontname="Arial", va="top")

title_theory = fig.text(RX, 0.13, "Theoretical Fixation Probability",
                         fontsize=14, fontweight="bold", fontname="Arial", va="top")
body_theory = fig.text(RX, 0.085, "", fontsize=11, fontname="Arial", va="top")


# ---------------------------------------------------------------------------
# Slider <-> TextBox syncing (only syncs values; simulation runs on button
# click so the UI stays responsive with large numbers of paths)
# ---------------------------------------------------------------------------
def make_slider_to_text(slider, textbox, is_int):
    def _on_slider_change(val):
        val = int(val) if is_int else round(val, 4)
        textbox.set_val(str(val))
    return _on_slider_change


def make_text_to_slider(slider, textbox, bounds, is_int):
    lo, hi = bounds

    def _on_text_submit(text):
        try:
            val = float(text)
        except ValueError:
            textbox.set_val(str(slider.val))
            return
        val = max(lo, min(hi, val))
        if is_int:
            val = int(round(val))
        slider.set_val(val)  # triggers slider->text sync above, keeps in sync

    return _on_text_submit


for key in PARAM_BOUNDS:
    is_int = key != "p0"
    sliders[key].on_changed(make_slider_to_text(sliders[key], textboxes[key], is_int))
    textboxes[key].on_submit(
        make_text_to_slider(sliders[key], textboxes[key], PARAM_BOUNDS[key], is_int)
    )


# ---------------------------------------------------------------------------
# Main run / redraw
# ---------------------------------------------------------------------------
def run(event=None):
    n_alleles = int(sliders["n_alleles"].val)
    generations = int(sliders["generations"].val)
    num_paths = int(sliders["num_paths"].val)
    p0 = float(sliders["p0"].val)

    trajectory = wright_fisher_simulate(n_alleles, p0, generations, num_paths)
    first_hit, is_fixed, is_lost = analyze_absorption(trajectory)

    n_fixed = int(is_fixed.sum())
    n_lost = int(is_lost.sum())
    n_segregating = num_paths - n_fixed - n_lost

    p_fixed = n_fixed / num_paths
    p_lost = n_lost / num_paths

    # ---------------- trajectory plot ----------------
    # For performance, only draw a random subset of ~100 trajectory LINES.
    # However, the fixation/loss dots are still plotted for EVERY path
    # (drawing a single scatter point is cheap even for thousands of paths).
    ax_traj.cla()
    ax_traj.set_box_aspect(0.5)
    ax_traj.set_xlabel("Generation")
    ax_traj.set_ylabel("Frequency of allele A")
    ax_traj.grid(alpha=0.3)

    n_lines_to_show = min(num_paths, 100)
    if num_paths > n_lines_to_show:
        display_rng = np.random.default_rng()
        line_indices = display_rng.choice(num_paths, size=n_lines_to_show, replace=False)
    else:
        line_indices = np.arange(num_paths)

    for i in line_indices:
        ax_traj.plot(trajectory[:, i], alpha=0.3, linewidth=0.7, color="steelblue")

    # scatter ALL fixation/loss events, regardless of which lines are drawn
    hit_idx = np.where(first_hit != -1)[0]
    if hit_idx.size > 0:
        x_dots = first_hit[hit_idx]
        y_dots = trajectory[x_dots, hit_idx]
        dot_colors = np.where(is_fixed[hit_idx], "green", "red")
        ax_traj.scatter(x_dots, y_dots, c=dot_colors, s=28,
                         edgecolors="black", linewidths=0.4, zorder=5)

    ax_traj.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax_traj.axhline(0.0, color="gray", linestyle="--", linewidth=0.8)

    ax_traj.set_xlim(0, generations)
    ax_traj.set_ylim(-0.02, 1.02)
    ax_traj.set_title(
        f"2N={n_alleles}, initial p(A)={p0:.2f}, {num_paths} paths total\n"
        f"(randomly showing {n_lines_to_show} paths as lines; "
        f"all fixation/loss points are marked)",
        fontsize=10,
    )
    ax_traj.legend(handles=legend_handles, loc="center left",
                   bbox_to_anchor=(1.01, 0.5), fontsize=9, frameon=True)

    # ---------------- two separate histograms ----------------
    fixed_gens = first_hit[is_fixed]
    lost_gens = first_hit[is_lost]
    max_bins = max(min(30, generations), 5)

    ax_hist_fixed.cla()
    ax_hist_fixed.set_box_aspect(1)
    ax_hist_fixed.set_xlabel("Generation of fixation")
    ax_hist_fixed.set_ylabel("Count")
    ax_hist_fixed.grid(alpha=0.3)
    ax_hist_fixed.set_title(f"Fixed (n={n_fixed})", color="darkgreen", fontsize=11)

    if fixed_gens.size > 0:
        bins_f = np.linspace(0, max(fixed_gens.max(), 1), max_bins)
        ax_hist_fixed.hist(fixed_gens, bins=bins_f, color="green", alpha=0.65)
        mean_f, std_f = fixed_gens.mean(), fixed_gens.std()
        ax_hist_fixed.axvline(mean_f, color="darkgreen", linestyle="--", linewidth=1.2)
        ax_hist_fixed.text(
            0.02, 0.97, f"mean={mean_f:.1f}\nstd={std_f:.1f}",
            transform=ax_hist_fixed.transAxes, fontsize=8, va="top", ha="left",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="gray"),
        )
    else:
        ax_hist_fixed.text(0.5, 0.5, "no fixed paths", transform=ax_hist_fixed.transAxes,
                            fontsize=9, ha="center", va="center", color="gray")

    ax_hist_lost.cla()
    ax_hist_lost.set_box_aspect(1)
    ax_hist_lost.set_xlabel("Generation of loss")
    ax_hist_lost.set_ylabel("Count")
    ax_hist_lost.grid(alpha=0.3)
    ax_hist_lost.set_title(f"Lost (n={n_lost})", color="darkred", fontsize=11)

    if lost_gens.size > 0:
        bins_l = np.linspace(0, max(lost_gens.max(), 1), max_bins)
        ax_hist_lost.hist(lost_gens, bins=bins_l, color="red", alpha=0.65)
        mean_l, std_l = lost_gens.mean(), lost_gens.std()
        ax_hist_lost.axvline(mean_l, color="darkred", linestyle="--", linewidth=1.2)
        ax_hist_lost.text(
            0.02, 0.97, f"mean={mean_l:.1f}\nstd={std_l:.1f}",
            transform=ax_hist_lost.transAxes, fontsize=8, va="top", ha="left",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="gray"),
        )
    else:
        ax_hist_lost.text(0.5, 0.5, "no lost paths", transform=ax_hist_lost.transAxes,
                           fontsize=9, ha="center", va="center", color="gray")

    # ---------------- bottom-right results panel ----------------
    body_counts.set_text(
        f"Total number of paths    = {num_paths}\n"
        f"Paths fixed              = {n_fixed}\n"
        f"Paths lost               = {n_lost}\n"
        f"Paths still segregating  = {n_segregating}"
    )

    body_prob.set_text(
        f"Probability of being fixed at generation {generations} = {p_fixed:.4f}\n"
        f"Probability of being lost at generation {generations} = {p_lost:.4f}"
    )

    body_theory.set_text(
        f"P(eventually fixed) = initial p(A) = {p0:.4f}"
    )

    fig.canvas.draw_idle()


button.on_clicked(run)

# initial draw so the figure isn't empty before the first click
run(None)

plt.show()