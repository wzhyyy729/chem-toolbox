# Chem Toolbox

**English** | [简体中文](README_zh.md)

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Dependencies](https://img.shields.io/badge/dependencies-none-3FB950)](#installation)
[![Tests](https://img.shields.io/badge/tests-78%20passed-3FB950)](#self-test)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

Four kinds of calculations a chemistry student does every week, bundled into one command-line tool: **formula analysis, solution preparation and pH, lab-data statistics and curve fitting, and spectral peak finding**.

Built on the **Python standard library only** — no third-party dependencies, nothing to `pip install`. Copy the folder onto any machine with Python and it runs.

> Started by a second-year chemistry undergraduate. The code was written with AI assistance, but **every formula was verified by hand and is covered by 78 automated assertions**. It is not fancy; the numbers are checked.

---

## Contents

- [Features](#features)
- [Why this exists](#why-this-exists)
- [Installation](#installation)
- [Usage](#usage)
- [The chemistry and maths behind it](#the-chemistry-and-maths-behind-it)
- [Self-test](#self-test)
- [Known limitations](#known-limitations)
- [Project layout](#project-layout)

---

## Features

| Command | What it does |
| :--- | :--- |
| `formula` | Parse a chemical formula → molar mass, element mass fractions, degree of unsaturation (DBE) |
| `solution prepare` | Solution preparation: mass ↔ concentration (accepts a formula and looks up molar mass for you) |
| `solution dilute` | Dilution law c₁V₁ = c₂V₂, leave exactly one term unknown |
| `solution ph` | pH of strong acids, strong bases, monoprotic weak acids and weak bases (exact solution, not the approximation) |
| `solution beer` | Beer–Lambert law A = εbc |
| `solution unknown` | Back-calculate an unknown concentration from a calibration curve |
| `stats summary` | Mean, sample/population standard deviation, standard error, range, median |
| `stats fit` | Least-squares linear fit: slope, intercept, R², parameter standard errors |
| `peak` | Peak finding: position, height, signal-to-noise ratio, full width at half maximum |

**Also handled**: nested brackets in formulas (`K3[Fe(CN)6]`), hydrates (`CuSO4.5H2O`), instrument exports with `#` comment lines (skipped automatically), and comma / tab / whitespace delimiters (auto-detected).

---

## Why this exists

Lab reports kept demanding the same calculations: how much solid to weigh out, the standard deviation of replicate titrations, the slope and R² of a calibration curve, the peak positions in a UV-Vis spectrum. Typing these into a calculator is slow and easy to get wrong.

More importantly, **every one of these formulas has conditions attached**. The strong-acid approximation fails below about 10⁻⁶ mol/L because water's self-ionisation can no longer be ignored. The weak-acid shortcut √(Kₐc) drifts noticeably once the degree of ionisation passes 5%. Plenty of scripts and online calculators hand you a number without mentioning any of that.

So this tool gives you **the result, the reasoning, and a warning when the assumptions do not hold**.

---

## Installation

Requires Python 3.10 or newer (the code uses `X | None` type annotations).

```bash
git clone https://github.com/wzhyyy729/chem-toolbox.git
cd chem-toolbox
python -m chem_toolbox --help
```

There is nothing to install. No `requirements.txt`.

---

## Usage

### Formula analysis

```console
$ python -m chem_toolbox formula "CuSO4.5H2O"
Formula       : CuSO4.5H2O
Composition   : H10  O9  Cu1  S1
Molar mass    : 249.6770 g/mol
Degree of unsaturation: not applicable - contains Cu, S
```

Nested brackets work too:

```console
$ python -m chem_toolbox formula "K3[Fe(CN)6]"
Molar mass    : 329.2470 g/mol
```

And organics get a real DBE:

```console
$ python -m chem_toolbox formula "C6H6"
Degree of unsaturation: 4.0000     # benzene: 3 double bonds + 1 ring
```

### pH

```console
$ python -m chem_toolbox solution ph --acid-kind weak-acid --pka 4.76 --conc 0.1
Concentration : 0.1000 mol/L
Ka            : 1.7378e-05
pH            : 2.883
Ionisation α  : 1.310%
Basis         : exact solution of x² + Ka·x - Ka·c = 0
```

Warnings fire when the assumptions break:

```console
$ python -m chem_toolbox solution ph --acid-kind strong-acid --conc 1e-8
pH            : 8.000
Warning       : below 1e-6 mol/L, water ionisation is no longer negligible
```

### Preparing a solution

```console
$ python -m chem_toolbox solution prepare --formula NaCl --conc 0.1 --volume 0.5
Molar mass    : 58.4400 g/mol (from formula NaCl)
Mass required : 2.9220 g
```

### Calibration curve

```console
$ python -m chem_toolbox stats fit --file examples/standard_curve.csv --x-col 0 --y-col 1
n             = 8
slope k       = 1.9957 ± 0.0096
intercept b   = 0.0020 ± 0.0024
R²            = 0.999861
```

### Peak finding

```console
$ python -m chem_toolbox peak --file examples/uvvis_spectrum.csv --x-col 0 --y-col 1
Noise level   : 0.0199
Threshold     : 3.0000 × noise = 0.0598

2 peaks found
  #        position        height     prominence     SNR        FWHM
  1         455.000        0.4420        0.3080      15.5      37.635
  2         605.000        0.2850        0.1510       7.6      44.359
```

### As a library

```python
from chem_toolbox import analyse_formula, ph_weak_acid, linear_fit

report = analyse_formula("C6H12O6")
print(report.molar_mass)            # 180.156

ph, alpha = ph_weak_acid(ka=1.75e-5, concentration=0.1)
print(round(ph, 4))                 # 2.8813

fit = linear_fit([0.1, 0.2, 0.3], [0.2, 0.4, 0.6])
print(fit.slope, fit.r_squared)     # 2.0 1.0
```

Add `--json` to any subcommand for machine-readable output.

---

## The chemistry and maths behind it

### Degree of unsaturation

$$\mathrm{DBE} = \frac{2C + 2 + N - H - X}{2}$$

Oxygen and sulfur do not appear. The formula only holds for **C, H, N, O and halogen** organic compounds — for anything containing a metal (say `K3[Fe(CN)6]`) the tool refuses to compute rather than returning a meaningless number. A non-integer result means the formula is mistyped; a negative one means the composition is impossible.

### Monoprotic weak acid, exact solution

With $[\mathrm{H^+}] = x$ from mass and charge balance:

$$K_a = \frac{x^2}{c - x} \quad\Longrightarrow\quad x = \frac{-K_a + \sqrt{K_a^2 + 4K_a c}}{2}$$

The common shortcut $x \approx \sqrt{K_a c}$ is only valid while the degree of ionisation stays small; above roughly 5% the error becomes visible, so this tool solves the quadratic and warns you.

### Least squares

$$k = \frac{\sum (x_i - \bar{x})(y_i - \bar{y})}{\sum (x_i - \bar{x})^2}, \qquad b = \bar{y} - k\bar{x}, \qquad R^2 = 1 - \frac{SS_{\text{res}}}{SS_{\text{tot}}}$$

With residual standard error $s = \sqrt{SS_{\text{res}}/(n-2)}$:

$$s_k = \frac{s}{\sqrt{\sum (x_i-\bar{x})^2}}, \qquad s_b = s\sqrt{\frac{1}{n} + \frac{\bar{x}^2}{\sum (x_i-\bar{x})^2}}$$

### Peak finding

The **baseline** is the median, not the mean, because peaks drag a mean upward while the median barely moves. **Noise** is estimated robustly from successive differences: $\hat{\sigma} \approx 1.4826 \times \mathrm{median}(|\Delta y|) / \sqrt{2}$, where 1.4826 converts a median absolute deviation into a standard deviation and the $\sqrt{2}$ accounts for differencing doubling the variance. A peak must be a local maximum whose prominence exceeds 3σ. **FWHM** is found by linear interpolation at half prominence; for a Gaussian this equals $2\sqrt{2\ln 2}\,\sigma \approx 2.3548\sigma$.

---

## Self-test

```console
$ python -m chem_toolbox.tests
==============================================================
chem_toolbox self-test
==============================================================
  [pass] Formula parsing
  [pass] Solutions and pH
  [pass] Statistics and fitting
  [pass] Peak finding
--------------------------------------------------------------
78 assertions passed, 0 failed
```

Every expected value can be checked against a textbook or by hand — this is not a "it ran without raising" test suite. Examples:

- Glucose 180.156 g/mol, copper sulfate pentahydrate 249.68 g/mol, `K3[Fe(CN)6]` 329.24 g/mol
- Benzene DBE = 4, cyclohexane = 1, ethane = 0; metals rejected
- 0.1 mol/L acetic acid (pKₐ 4.76) → pH 2.8813 ± 0.0005; 0.1 mol/L ammonia (K_b 1.8×10⁻⁵) → pH 11.1248
- A synthetic Gaussian peak must yield FWHM within 1.5 of 2.3548σ
- Invalid input (empty formula, lowercase element symbols, unbalanced brackets, negative concentrations, constant x) must raise, not silently return nonsense

Works without pytest.

---

## Known limitations

Documenting the edges is more useful than pretending there are none:

- The atomic mass table covers roughly 40 common elements (H through Pb). Add a line to `ATOMIC_MASS` in `chem_toolbox/core.py` for anything else.
- **Monoprotic acids and bases only.** Polyprotic systems such as H₃PO₄ or H₂CO₃ need simultaneous stepwise equilibria.
- **No activity coefficients** — concentrations stand in for activities, so results drift at high ionic strength (roughly ≥ 0.1 mol/L strong electrolyte).
- **No peak-shape fitting**, just local maxima with a threshold. The x-axis must be increasing; decreasing series are reversed, but no interpolation is done.
- **No plotting**, deliberately: matplotlib would break the zero-dependency property, which is the main reason this tool is easy to run anywhere. Export the CSV to Excel or Origin instead.

---

## Project layout

```
chem-toolbox/
├── chem_toolbox/
│   ├── __init__.py     public API (from chem_toolbox import ...)
│   ├── __main__.py     CLI: argument parsing and output formatting
│   ├── core.py         pure computation - no printing, no file I/O
│   └── tests.py        78 assertions, no pytest required
├── examples/
│   ├── standard_curve.csv    spectrophotometric calibration data
│   ├── titration.csv         replicate titration volumes
│   └── uvvis_spectrum.csv    UV-Vis spectrum with two bands
├── LICENSE
└── README.md
```

`core.py` and `__main__.py` are kept apart on purpose: the algorithms do no input or output, so they can be driven from the CLI, imported into a script, or asserted against directly in tests.

---

## License

[MIT](LICENSE)
