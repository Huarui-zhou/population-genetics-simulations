# Population Genetics Simulations

Interactive Python simulations for exploring fundamental concepts in population genetics.

This repository contains two visual and interactive population-genetics models:

1. A deterministic model of random mating with viability selection and two-way mutation.
2. A stochastic Wright–Fisher model for genetic drift in a finite population.

The simulations are intended primarily for teaching, visualization, and exploratory use.

---

## Repository Structure

```text
population-genetics-simulations/
├── README.md
└── codes/
    ├── Random_mating.py
    └── Wright_Fisher_drift.py
```

---

## 1. Random Mating with Selection and Two-Way Mutation

**File:** `codes/Random_mating.py`

This program is an interactive deterministic diploid population-genetics model incorporating:

- Random mating
- Viability selection
- Forward mutation: A → a
- Reverse mutation: a → A

### Model life cycle

In each generation:

1. Viability selection acts on the current zygote genotype frequencies.
2. Selected adults produce gametes.
3. Mutation occurs in the gamete pool.
4. Gametes unite at random.
5. The next generation is formed under Hardy–Weinberg proportions.

The genotype fitnesses are defined as:

$$
w_{AA} = 1 - s_{AA}
$$

$$
w_{Aa} = 1 - s_{Aa}
$$

$$
w_{aa} = 1 - s_{aa}
$$

### Adjustable parameters

Users can modify:

- Initial frequency of genotype `AA`
- Initial frequency of genotype `Aa`
- Selection coefficient $s_{AA}$
- Selection coefficient $s_{Aa}$
- Selection coefficient $s_{aa}$
- Forward mutation rate $u$, for A → a
- Reverse mutation rate $v$, for a → A
- Number of generations

The initial frequency of genotype `aa` is calculated automatically from:

$$
f_{aa} = 1 - f_{AA} - f_{Aa}
$$

### Output

The simulation displays:

- Genotype-frequency dynamics for `AA`, `Aa`, and `aa`
- Allele-frequency dynamics for `A` and `a`
- Initial and final genotype frequencies
- Initial and final allele frequencies
- Simulation parameters

This model is deterministic: the same parameter values produce the same population trajectory.

---

## 2. Genetic Drift Model — Wright–Fisher Model

**File:** `codes/Wright_Fisher_drift.py`

This program is an interactive stochastic simulation of genetic drift based on the neutral Wright–Fisher model.

In a finite population, allele frequencies fluctuate from generation to generation because of random sampling. In each generation, the number of copies of allele `A` is sampled from a binomial distribution determined by the allele frequency in the previous generation.

For a population containing $2N$ allele copies,

$$
X_{t+1} \sim \mathrm{Binomial}(2N, p_t)
$$

and

$$
p_{t+1} = \frac{X_{t+1}}{2N}
$$

where $p_t$ is the frequency of allele `A` in generation $t$.

### Adjustable parameters

Users can modify:

- Population size, represented as the number of allele copies (`2N`)
- Number of generations
- Number of independent simulation paths
- Initial allele frequency, $p(A)$

### Output

The simulation displays:

- Multiple allele-frequency trajectories
- Fixation events
- Loss events
- Distribution of fixation times
- Distribution of loss times
- Number of populations fixed
- Number of populations lost
- Number of populations still segregating
- Simulated probability of fixation
- Simulated probability of loss
- Theoretical fixation probability

For the neutral Wright–Fisher model,

$$
P(\text{eventual fixation of } A) = p_0
$$

where $p_0$ is the initial frequency of allele `A`.

Because this model is stochastic, repeated simulations with the same parameters can produce different trajectories.

---

## Comparison of the Two Models

| Model | Type | Main processes |
|---|---|---|
| `Random_mating.py` | Deterministic | Random mating, viability selection, mutation |
| `Wright_Fisher_drift.py` | Stochastic | Genetic drift caused by finite-population sampling |

The random-mating model illustrates how deterministic evolutionary forces change genotype and allele frequencies.

The Wright–Fisher model illustrates how random sampling alone can cause an allele to fluctuate in frequency and eventually become fixed or lost, even in the absence of selection.

Together, the two simulations provide a simple comparison between deterministic and stochastic processes in population genetics.

---

## Requirements

The simulations require Python 3 and the following Python packages:

```bash
pip install numpy matplotlib
```

For interactive use in JupyterLab, install:

```bash
pip install ipympl
```

Then enable the interactive Matplotlib backend:

```python
%matplotlib widget
```

Do not use the default inline Matplotlib backend if you want the interactive controls to remain functional.

---

## Running the Simulations

Clone or download this repository, then run the programs from the repository root.

### Random mating, selection, and mutation model

```bash
python codes/Random_mating.py
```

### Wright–Fisher genetic drift model

```bash
python codes/Wright_Fisher_drift.py
```

Both programs use Matplotlib's built-in interactive widgets.

---

## Educational Purpose

These simulations are intended to help visualize and explore concepts such as:

- Hardy–Weinberg random mating
- Genotype and allele frequencies
- Viability selection
- Mutation
- Genetic drift
- Finite-population sampling
- Allele fixation
- Allele loss
- Deterministic versus stochastic evolutionary dynamics

The models are deliberately simplified so that individual population-genetic processes can be explored interactively.

---

## Author

**Huarui Zhou**

Simulation concepts, model design, interface design, and scientific specifications by Huarui Zhou.

Generative AI was used as a coding assistant during implementation, debugging, and code refinement.

---

## Disclaimer

These programs are intended for educational and exploratory use.

They are simplified population-genetic models and are not intended to represent all biological processes acting in natural populations.
