# Population Genetics Models: Selection, Mutation, and Genetic Drift

Interactive Python simulations for exploring fundamental concepts in population genetics.

These programs were designed as visual and interactive tools for teaching and exploring how genotype and allele frequencies change across generations under selection, mutation, random mating, and genetic drift.

## Simulations

### 1. Random Mating with Selection and Two-Way Mutation

`HW.py`

An interactive deterministic diploid population-genetics model incorporating:

* Viability selection
* Random mating
* Forward mutation \(A \rightarrow a\)
* Reverse mutation \(a \rightarrow A\)

The mating cycle in each generation is:

1. Viability selection acts on the current genotype frequencies.
2. Selected adults produce gametes.
3. Mutation changes allele frequencies in the gamete pool.
4. Gametes unite randomly.
5. The next generation is formed under Hardy–Weinberg proportions.

Genotype fitnesses are defined as

$$
w_{AA}=1-s_{AA}
$$

$$
w_{Aa}=1-s_{Aa}
$$

$$
w_{aa}=1-s_{aa}
$$

Users can modify:

* Initial genotype frequencies
* Selection coefficients
* Forward mutation rate \(u\)
* Reverse mutation rate \(v\)
* Number of generations

The program displays genotype-frequency and allele-frequency dynamics through time.

---

### 2. Genetic Drift Model (Wright–Fisher Model)

`WF.py`

An interactive stochastic simulation of genetic drift based on the neutral Wright–Fisher model.

The model represents a finite population in which allele frequencies change from generation to generation because of random sampling. For each generation, the number of copies of allele **A** is sampled from a binomial distribution determined by its frequency in the previous generation.

Users can modify:

* Population size (`2N`)
* Number of generations
* Number of simulation paths
* Initial allele frequency, `p(A)`

The simulation displays:

* Allele-frequency trajectories
* Fixation and loss events
* Distribution of fixation times
* Distribution of loss times
* Numbers of populations fixed, lost, or still segregating
* Simulated fixation and loss probabilities
* The theoretical fixation probability

Under the neutral Wright–Fisher model,

$$
P(\text{eventual fixation of A}) = p_0
$$

where \(p_0\) is the initial frequency of allele **A**.

---

## Requirements

The simulations require Python 3 and the following packages:

```bash
pip install numpy matplotlib
```

For interactive use inside JupyterLab:

```bash
pip install ipympl
```

Then enable the interactive Matplotlib backend:

```python
%matplotlib widget
```

The default inline backend should not be used because the interactive controls require a live Matplotlib canvas.

---

## Running the Simulations

Run the random-mating selection and mutation model with:

```bash
python HW.py
```

Run the genetic drift model with:

```bash
python WF.py
```

The programs use Matplotlib's built-in interactive widgets, including sliders, text boxes, and buttons.

---

## Purpose

These simulations are intended primarily for:

* Teaching population genetics
* Visualizing evolutionary dynamics
* Exploring selection, mutation, random mating, and genetic drift
* Comparing deterministic and stochastic population-genetic models
* Developing intuition about allele-frequency change, fixation, and loss

`HW.py` is a deterministic model of genotype and allele-frequency change under viability selection, mutation, and random mating.

`WF.py` is a stochastic model of genetic drift in a finite population.

Together, the two simulations illustrate the contrast between deterministic evolutionary forces and stochastic changes caused by finite-population sampling.

---

## Author

**Huarui Zhou**

Simulation concepts, model design, interface design, and scientific specifications by Huarui Zhou.

Generative AI was used as a coding assistant during implementation, debugging, and code refinement.

---

## Disclaimer

These programs are intended for educational and exploratory use. They are simplified population-genetic models and are not intended to represent every biological process affecting natural populations.

---

## License

See the `LICENSE` file for licensing information.
