# Validation of Variance for Sobol' and CvM rank-based estimators

Companion code to the paper

> **A Martingale Approach To Fluctuations of Rank Estimators in Sensitivity Analysis**
> Reda Chhaibi, Fabrice Gamboa and Clément Pellegrini (2026).

The repository contains the numerical validation of Appendix C of the paper: for several models
we compare our closed-form asymptotic variances with crude Monte Carlo estimates.

## The paper in a nutshell

Given an i.i.d. sample (X<sub>i</sub>, Y<sub>i</sub>), i = 1..n, sensitivity indices measure how
much an input X influences an output Y. The paper studies two of them:

- the **Sobol' index**: Var(E[Y | X]) / Var(Y);
- the **Cramér–von Mises (CvM) index**, which compares the distribution of Y with its
  conditional distribution given X.

Both can be estimated from a single sample with the rank-based estimators of Chatterjee: sort the
sample by X, and pair each observation with its right neighbour. Writing Y = f(X, ε) with ε
independent of X, the paper splits these estimators into three parts: a part that depends on
the X<sub>i</sub>'s only, a martingale part, and a remainder. This gives consistency and central
limit theorems under minimal regularity assumptions:

- **Theorem 2.1** (Sobol'): the estimator is asymptotically Gaussian at rate √n, with an explicit,
  structured asymptotic variance σ²<sub>Sobol'</sub>;
- **Theorem 2.2**: a multivariate extension to several functions of the same scalar input;
- **Theorem 2.3** (CvM): the same for the CvM estimator, with what is, to our knowledge, the first
  explicit formula for its asymptotic variance σ²<sub>CvM</sub>.

Earlier variance formulas in the literature turned out to contain errors, so we check ours
numerically. That check is what this repository does.

## Content

```bash
./
|-- README.md
|-- ipynb/
    |-- joint_code_AOS_CGP_v1.0.ipynb   # The notebook: reproduces every figure of Appendix C
    |-- sobol_ishigami_results.npz      # Saved Monte Carlo results, Section 1 (Fig. C.1)
    |-- gaussian_results.npz            # Saved Monte Carlo results, Section 2 (Figs. C.2, C.3)
    |-- random_support_results.npz      # Saved Monte Carlo results, Section 3 (Fig. C.4)
```

### What the notebook does

| Notebook section | Model | What is computed | Paper |
|---|---|---|---|
| 1. Ishigami | Ishigami function, inputs uniform on (−π, π) | Sobol' index of X<sub>1</sub> and its asymptotic variance, derived **symbolically** with SymPy from the matrices of Eqs. (2.7)–(2.15), then compared with Monte Carlo for a = 10 and b from 0 to 2 | Fig. C.1 |
| 2.1 Gaussian | Y = ρX + √(1−ρ²) ε | Sobol' index and its asymptotic variance in closed form, compared with Monte Carlo for ρ from −1 to 1 | Fig. C.2 |
| 2.2 Gaussian | same | Sobol' and CvM indices side by side | Fig. C.3 |
| 3. Random support | Y = Xε, with X ~ Beta(α, 1) and ε uniform on (0, 1) | Covariance kernels of Appendix C.3.2, integrated numerically with `scipy.integrate.dblquad` to get the CvM asymptotic variance, compared with Monte Carlo for α from 0.5 to 7 | Fig. C.4 |

In each Monte Carlo experiment, the estimator is computed `n_repetitions` times on independent
samples of size `n_samples`. The empirical variance, multiplied by `n_samples`, is then plotted
against the theoretical asymptotic variance:

| Experiment | `n_samples` | `n_repetitions` | Grid points |
|---|---|---|---|
| Ishigami (Fig. C.1) | 10 000 | 10 000 | 201 values of b |
| Gaussian (Fig. C.2) | 10 000 | 10 000 | 201 values of ρ |
| Random support (Fig. C.4) | 20 000 | 20 000 | 93 values of α |

The notebook also checks its own implementation. The vectorised, multithreaded estimator is
compared with a literal implementation based on `scipy.stats.pearsonr`, and the floating-point
kernels are compared with their SymPy transcriptions.

## Installation

1. Create and activate a virtual environment (on Debian/Ubuntu, run `sudo apt install python3-venv` first if needed)

```bash
$ python3 -m venv .venv
$ source .venv/bin/activate
```

2. Install the dependencies

```bash
$ pip install --upgrade pip
$ pip install numpy scipy sympy matplotlib jupyter
```

The code was last run with numpy 1.26, scipy 1.11, sympy 1.12 and matplotlib 3.6.

3. (Optional) Register the environment as a Jupyter kernel

```bash
$ pip install ipykernel
$ python -m ipykernel install --user --name=.venv
```

(see https://janakiev.com/blog/jupyter-virtual-envs/ for details)

## Reproducing the figures

Launch Jupyter **from the `ipynb/` folder**, because the notebook reads and writes the `.npz` files
by relative path:

```bash
$ cd ipynb
$ jupyter notebook joint_code_AOS_CGP_v1.0.ipynb
```

Each figure is drawn from its `.npz` file, not from variables in memory. There are two ways to use
the notebook:

- **Fast path (minutes): redraw the figures from the saved results.** Run all cells *except* the
  three main Monte Carlo loops: the cells that start with `# Main loop` in Sections 1 and 3, and with
  `# Monte Carlo experiments` in Section 2. The symbolic computations of Section 1.1 still take a
  little while.
- **Full path: redo the Monte Carlo experiments.** Run every cell. Each experiment
  draws and sorts between 2×10¹⁰ and 4×10¹⁰ sample points, so it runs for a
  long time even when multithreaded. Each main loop **overwrites** its `.npz` file. To try things out, reduce
  `n_repetitions` / `n_samples`, or use the coarser grids left commented in the code.

Running time depends on the `hardware_options` cell at the top of the notebook:

```python
hardware_options = {
    "n_threads": 16,    # threads the parameter grid is spread over; 1 for a serial run
    "batch_size": 25,   # repetitions per block; small so that a block stays in cache
}
```

These settings affect speed only, never the results. Their defaults were tuned for a
14-core / 20-thread laptop with 30 GB of memory.

The figures are saved as `sobol_ishigami.png`, `sobol_gaussian.png`, `cvm_vs_sobol_gaussian.png` and
`cvm_random_support.png`. Image files are excluded by `.gitignore`.

## Citation

If you use this code, please cite the paper:

```bibtex
@article{ChhaibiGamboaPellegrini2026,
  title   = {A Martingale Approach To Fluctuations of Rank Estimators in Sensitivity Analysis},
  author  = {Chhaibi, Reda and Gamboa, Fabrice and Pellegrini, Cl{\'e}ment},
  year    = {2026},
  note    = {Preprint}
}
```

## Contact

- Reda Chhaibi, Université Côte d'Azur, LJAD — reda.chhaibi@univ-cotedazur.fr
- Fabrice Gamboa, Institut de Mathématiques de Toulouse and ANITI — fabrice.gamboa@math.univ-toulouse.fr
- Clément Pellegrini, Institut de Mathématiques de Toulouse and ANITI — clement.pellegrini@math.univ-toulouse.fr
