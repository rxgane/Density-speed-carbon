# Cost, Carbon, and the Speed of Replenishment in Distribution Network Design

Replication materials for the analytical paper *"Cost, Carbon, and the Speed of
Replenishment in Distribution Network Design"* (Ram Ganeshan, Raymond A. Mason
School of Business, William & Mary).

A continuous-approximation location–inventory model in which a firm chooses **network density**
(the number of distribution centers, `n`) and **replenishment speed** (`v`) to
minimize cost subject to a carbon cap and a service requirement. Carbon splits
into a transport component and a holding component; pricing carbon selects which
operational substitute — stock vs. speed for availability, facilities vs.
last-mile for proximity — is lower-emission. Key results: a closed-form speed
threshold, a pooling-adjusted density threshold, their complementarity, and a
chance-constrained (high-confidence) carbon cap.

## Repository contents

| File | Description |
|------|-------------|
| `model.tex`, `references.bib` | The paper source and bibliography |
| `model.pdf` | Compiled paper |
| `numerics.py` | Core model: cost/carbon functions and optimizers |
| `computational_experiments.py` | Calibration, archetypes, value-of-joint-redesign, robustness |
| `chance_experiment.py` | Chance-constrained carbon-cap illustrations |
| `figures.py`, `gen_chance_fig.py`, `gen_regimes_fig.py` | Figure generation |
| `fig_regimes.pdf`, `fig_frontier.pdf`, `fig_chance.pdf` | Generated figures |
| `verify_*.wl`, `check_*.{wl,py}` | Symbolic/numeric checks (Mathematica + Python) |
| `results.txt` | Saved numerical output |

## Reproducing

Numerics require **NumPy** (`pip install numpy`); everything else is standard
library (`statistics.NormalDist`, no SciPy). With Python 3:

```bash
python3 numerics.py
python3 computational_experiments.py
python3 chance_experiment.py
```

Figures:

```bash
python3 gen_regimes_fig.py
python3 gen_chance_fig.py
python3 figures.py
```

Paper:

```bash
pdflatex model && bibtex model && pdflatex model && pdflatex model
```

Symbolic checks use `wolframscript` (Mathematica), e.g. `wolframscript -file verify_hessian.wl`.

> Note: copyrighted third-party journal PDFs used as references during the
> project are intentionally excluded from this repository.
