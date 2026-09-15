"""
Interactive diploid random-mating model with viability selection
and two-way mutation (v3).

Life cycle in each generation
-----------------------------
1. Viability selection acts on the current zygote genotype frequencies.
2. Selected adults produce gametes.
3. Mutation occurs in the gametes:
       A -> a at rate u
       a -> A at rate v
4. Gametes unite at random, producing Hardy-Weinberg zygote frequencies
   in the next generation.

Selection coefficients are defined by
       w_AA = 1 - s_AA
       w_Aa = 1 - s_Aa
       w_aa = 1 - s_aa

Only f_AA and f_Aa are editable - f_aa is filled in automatically as
1 - f_AA - f_Aa (displayed, not an input box).

All parameters use plain number input boxes (no sliders) - type a value
and press Enter.

Run as a standalone script:
    python random_mating_selection_mutation_v3.py

Or in JupyterLab:
    %matplotlib widget
    # requires: pip install ipympl
    # then paste and run this script in a cell

Do not use the default inline backend if you want the controls to remain
interactive.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox

plt.rcParams["font.family"] = "Arial"


# ---------------------------------------------------------------------------
# Core model
# ---------------------------------------------------------------------------
def normalize_genotype_frequencies(f_AA, f_Aa, f_aa):
    """Normalize three nonnegative genotype frequencies so they sum to 1."""
    frequencies = np.array([f_AA, f_Aa, f_aa], dtype=float)
    frequencies = np.clip(frequencies, 0.0, None)
    total = frequencies.sum()
    if total <= 0:
        raise ValueError("At least one initial genotype frequency must be positive.")
    return frequencies / total


def simulate_random_mating(f_AA0, f_Aa0, f_aa0, s_AA, s_Aa, s_aa, u, v, generations):
    """
    Simulate genotype and allele-frequency dynamics.

    Returns
    -------
    genotype_history : ndarray, shape (generations + 1, 3)   columns AA, Aa, aa
    allele_history   : ndarray, shape (generations + 1, 2)   columns A, a
    mean_fitness_history : ndarray, shape (generations,)
    """
    genotype_history = np.zeros((generations + 1, 3), dtype=float)
    allele_history = np.zeros((generations + 1, 2), dtype=float)
    mean_fitness_history = np.zeros(generations, dtype=float)

    genotype_history[0] = normalize_genotype_frequencies(f_AA0, f_Aa0, f_aa0)

    w = np.array([max(0.0, 1.0 - s_AA), max(0.0, 1.0 - s_Aa), max(0.0, 1.0 - s_aa)])

    for t in range(generations + 1):
        f_AA, f_Aa, f_aa = genotype_history[t]
        p_A = f_AA + 0.5 * f_Aa
        p_a = f_aa + 0.5 * f_Aa
        allele_history[t] = [p_A, p_a]

        if t == generations:
            break

        weighted = genotype_history[t] * w
        mean_fitness = weighted.sum()
        mean_fitness_history[t] = mean_fitness

        if mean_fitness <= 0:
            raise ValueError(
                "All genotypes have zero fitness. "
                "At least one selection coefficient must be less than 1."
            )

        selected = weighted / mean_fitness
        selected_AA, selected_Aa, selected_aa = selected

        p_selected = selected_AA + 0.5 * selected_Aa
        p_mutated = p_selected * (1.0 - u) + (1.0 - p_selected) * v
        p_mutated = np.clip(p_mutated, 0.0, 1.0)
        q_mutated = 1.0 - p_mutated

        genotype_history[t + 1] = [p_mutated**2, 2.0 * p_mutated * q_mutated, q_mutated**2]

    return genotype_history, allele_history, mean_fitness_history


# ---------------------------------------------------------------------------
# Initial parameters
# ---------------------------------------------------------------------------
init_f_AA = 0.25
init_f_Aa = 0.50
# f_aa is auto-filled = 1 - f_AA - f_Aa = 0.25

init_s_AA = 0.00
init_s_Aa = 0.10
init_s_aa = 0.20

init_u = 0.001
init_v = 0.001
init_generations = 100

BOUNDS = {
    "f_AA": (0.0, 1.0),
    "f_Aa": (0.0, 1.0),
    "s_AA": (0.0, 1.0),
    "s_Aa": (0.0, 1.0),
    "s_aa": (0.0, 1.0),
    "u": (0.0, 0.5),
    "v": (0.0, 0.5),
    "generations": (1, 5000),
}

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(15, 12))

fig.text(0.27, 0.975, "Random Mating with Selection and Two-Way Mutation",
          fontsize=15, fontweight="bold", ha="center", va="bottom")
fig.text(
    0.475, 0.975,
    "Designed by Huarui Zhou",
    fontsize=9,
    fontstyle="italic",
    ha="left",
    va="bottom",
    color="dimgray"
)

# left column: two 4:3 subplots
ax_genotypes = fig.add_axes([0.06, 0.55, 0.42, 0.30])
ax_alleles = fig.add_axes([0.06, 0.10, 0.42, 0.30])
ax_genotypes.set_box_aspect(3 / 4)
ax_alleles.set_box_aspect(3 / 4)

fig.text(0.27, 0.875, "Genotype Frequencies (AA, Aa, aa)",
          fontsize=13, fontweight="bold", ha="center", va="bottom")
fig.text(0.27, 0.425, "Allele Frequencies (A, a)",
          fontsize=13, fontweight="bold", ha="center", va="bottom")

for ax in (ax_genotypes, ax_alleles):
    ax.set_xlabel("Generation")
    ax.set_ylabel("Frequency")
    ax.set_ylim(-0.02, 1.02)

# ---------------------------------------------------------------------------
# Parameter controls (text boxes only, no sliders)
# ---------------------------------------------------------------------------
RX = 0.54
LABEL_W = 0.19

textbox_rows = [
    ("Initial frequency AA (f_AA)", "f_AA", init_f_AA, 0.93),
    ("Initial frequency Aa (f_Aa)", "f_Aa", init_f_Aa, 0.887),
    # f_aa is auto-filled - handled separately below (not a real textbox)
    ("Selection coefficient s_AA", "s_AA", init_s_AA, 0.801),
    ("Selection coefficient s_Aa", "s_Aa", init_s_Aa, 0.758),
    ("Selection coefficient s_aa", "s_aa", init_s_aa, 0.715),
    ("Mutation rate A \u2192 a (u)", "u", init_u, 0.672),
    ("Mutation rate a \u2192 A (v)", "v", init_v, 0.629),
    ("Number of generations", "generations", init_generations, 0.586),
]

textboxes = {}
for label_txt, key, init_val, y in textbox_rows:
    fig.text(RX, y, label_txt, fontsize=11, va="center")
    ax_t = fig.add_axes([RX + LABEL_W, y - 0.015, 0.10, 0.03])
    textboxes[key] = TextBox(ax_t, "", initial=str(init_val))

# auto-filled f_aa display (not an input)
fig.text(RX, 0.844, "Initial frequency aa (auto)", fontsize=11, va="center")
f_aa_display = fig.text(RX + LABEL_W, 0.844, "", fontsize=11, va="center",
                         fontweight="bold", color="#1a5276")


def clamp(value, lo, hi):
    return max(lo, min(hi, value))


def refresh_f_aa_display(event=None):
    try:
        f_AA_val = clamp(float(textboxes["f_AA"].text), *BOUNDS["f_AA"])
    except ValueError:
        f_AA_val = init_f_AA
    try:
        f_Aa_val = clamp(float(textboxes["f_Aa"].text), *BOUNDS["f_Aa"])
    except ValueError:
        f_Aa_val = init_f_Aa

    f_aa_val = max(0.0, 1.0 - f_AA_val - f_Aa_val)
    if f_AA_val + f_Aa_val > 1.0:
        f_aa_display.set_text(f"{f_aa_val:.4f}  (sum>1, auto-renorm.)")
    else:
        f_aa_display.set_text(f"{f_aa_val:.4f}")
    fig.canvas.draw_idle()


textboxes["f_AA"].on_submit(refresh_f_aa_display)
textboxes["f_Aa"].on_submit(refresh_f_aa_display)
refresh_f_aa_display()

ax_button = fig.add_axes([RX, 0.49, 0.20, 0.045])
button = Button(ax_button, "Run Simulation")

# ---------------------------------------------------------------------------
# Results panel: simulation parameters (text) + two comparison tables
# ---------------------------------------------------------------------------
title_params = fig.text(RX, 0.435, "Simulation Parameters", fontsize=14,
                         fontweight="bold", va="top")
body_params = fig.text(RX, 0.40, "", fontsize=10.5, va="top")

title_geno = fig.text(RX, 0.308, "Genotype Frequencies", fontsize=14,
                       fontweight="bold", va="top")
ax_table_geno = fig.add_axes([RX, 0.18, 0.22, 0.115])
ax_table_geno.axis("off")

title_allele = fig.text(RX, 0.18, "Allele Frequencies", fontsize=14,
                         fontweight="bold", va="top")
ax_table_allele = fig.add_axes([RX, 0.07, 0.22, 0.075])
ax_table_allele.axis("off")


# ---------------------------------------------------------------------------
# Main redraw
# ---------------------------------------------------------------------------
def run(event=None):
    try:
        f_AA0 = clamp(float(textboxes["f_AA"].text), *BOUNDS["f_AA"])
        f_Aa0 = clamp(float(textboxes["f_Aa"].text), *BOUNDS["f_Aa"])
        s_AA = clamp(float(textboxes["s_AA"].text), *BOUNDS["s_AA"])
        s_Aa = clamp(float(textboxes["s_Aa"].text), *BOUNDS["s_Aa"])
        s_aa = clamp(float(textboxes["s_aa"].text), *BOUNDS["s_aa"])
        u = clamp(float(textboxes["u"].text), *BOUNDS["u"])
        v = clamp(float(textboxes["v"].text), *BOUNDS["v"])
        generations = int(clamp(float(textboxes["generations"].text), *BOUNDS["generations"]))
    except ValueError:
        body_params.set_text("Error: please enter valid numbers in all fields.")
        fig.canvas.draw_idle()
        return

    f_aa0 = max(0.0, 1.0 - f_AA0 - f_Aa0)
    refresh_f_aa_display()

    try:
        genotype_history, allele_history, mean_fitness_history = simulate_random_mating(
            f_AA0=f_AA0, f_Aa0=f_Aa0, f_aa0=f_aa0,
            s_AA=s_AA, s_Aa=s_Aa, s_aa=s_aa,
            u=u, v=v, generations=generations,
        )
    except ValueError as error:
        body_params.set_text(f"Error:\n{error}")
        fig.canvas.draw_idle()
        return

    x = np.arange(generations + 1)

    # ---------------- genotype frequency plot ----------------
    ax_genotypes.cla()
    ax_genotypes.set_box_aspect(3 / 4)
    ax_genotypes.plot(x, genotype_history[:, 0], linewidth=2.0, color="steelblue", label="AA")
    ax_genotypes.plot(x, genotype_history[:, 1], linewidth=2.0, color="seagreen", label="Aa")
    ax_genotypes.plot(x, genotype_history[:, 2], linewidth=2.0, color="firebrick", label="aa")
    ax_genotypes.set_xlabel("Generation")
    ax_genotypes.set_ylabel("Genotype frequency")
    ax_genotypes.set_xlim(0, generations)
    ax_genotypes.set_ylim(-0.02, 1.02)
    ax_genotypes.legend(loc="upper right", fontsize=10)

    # ---------------- allele frequency plot ----------------
    ax_alleles.cla()
    ax_alleles.set_box_aspect(3 / 4)
    ax_alleles.plot(x, allele_history[:, 0], linewidth=2.2, color="darkorange", label="A")
    ax_alleles.plot(x, allele_history[:, 1], linewidth=2.2, color="purple", label="a")
    ax_alleles.set_xlabel("Generation")
    ax_alleles.set_ylabel("Allele frequency")
    ax_alleles.set_xlim(0, generations)
    ax_alleles.set_ylim(-0.02, 1.02)
    ax_alleles.legend(loc="upper right", fontsize=10)

    # ---------------- results panel: simulation parameters ----------------

    body_params.set_text(
        f"Selection coefficients\n"
        f"  s_AA = {s_AA:.4f},  s_Aa = {s_Aa:.4f},  s_aa = {s_aa:.4f} \n"
        f"Mutation rates\n"
        f"  A \u2192 a (u) = {u:.4f}, a \u2192 A (v) = {v:.4f}\n"
        f"Generations simulated = {generations}"
    )

    # ---------------- results panel: genotype frequency table ----------------
    ax_table_geno.cla()
    ax_table_geno.axis("off")
    geno_init = genotype_history[0]
    geno_final = genotype_history[-1]
    tbl_geno = ax_table_geno.table(
        cellText=[
            [f"{geno_init[0]:.4f}", f"{geno_final[0]:.4f}"],
            [f"{geno_init[1]:.4f}", f"{geno_final[1]:.4f}"],
            [f"{geno_init[2]:.4f}", f"{geno_final[2]:.4f}"],
        ],
        rowLabels=["AA", "Aa", "aa"],
        colLabels=["Initial", "Final"],
        cellLoc="center",
        loc="center",
    )
    tbl_geno.auto_set_font_size(False)
    tbl_geno.set_fontsize(10)
    tbl_geno.scale(1.0, 1.6)

    # ---------------- results panel: allele frequency table ----------------
    ax_table_allele.cla()
    ax_table_allele.axis("off")
    allele_init = allele_history[0]
    allele_final = allele_history[-1]
    tbl_allele = ax_table_allele.table(
        cellText=[
            [f"{allele_init[0]:.4f}", f"{allele_final[0]:.4f}"],
            [f"{allele_init[1]:.4f}", f"{allele_final[1]:.4f}"],
        ],
        rowLabels=["A", "a"],
        colLabels=["Initial", "Final"],
        cellLoc="center",
        loc="center",
    )
    tbl_allele.auto_set_font_size(False)
    tbl_allele.set_fontsize(10)
    tbl_allele.scale(1.0, 1.6)

    fig.canvas.draw_idle()


button.on_clicked(run)
run(None)

plt.show()