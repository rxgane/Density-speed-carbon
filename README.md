# Cost, Carbon, and the Speed of Replenishment in Distribution Network Design

Paper source and computations for *"Cost, Carbon, and the Speed of Replenishment in Distribution
Network Design"* (Ram Ganeshan, Raymond A. Mason School of Business, William & Mary), version 2,
prepared for *Transportation Research Part E*.

A continuous-approximation location–inventory model in which a firm chooses **network density**
(the number of distribution centers, `n`) and **replenishment speed** (`v`) to minimize cost under a
carbon price or cap and a service requirement. Carbon arises from facilities, holding inventory,
line-haul freight, and last-mile delivery. The joint design has a single best solution at every carbon
price, and the carbon-minimizing direction of both decisions follows from four carbon-intensity-per-dollar
ratios and the share of a distribution center's marginal cost that is lost pooling. A calibration to
public US data shows how the direction depends on how speed is bought and how the last mile is driven.

## Repository contents

| File | Description |
|------|-------------|
| `network_design_v2.tex`, `references.bib` | Manuscript source (anonymized for review) and bibliography |
| `network_design_v2.pdf` | Compiled manuscript |
| `title_page.tex`, `cover_letter.tex`, `highlights.tex` | Title page with declarations, cover letter, highlights |
| `revision/calibration_v2.py` | Carbon intensity per dollar of each activity from public sources; product classification; carbon share of priced cost |
| `revision/application_v2.py` | Calibrated national network in physical units (designs, carbon-price path, integer rounding) |
| `revision/v2_checks.py` | Checks of the analytical results against brute-force optimization |
| `revision/*_output.txt` | Saved output of the three scripts above |
| `revision/v2_analysis.md` | Statements and proofs of the version 2 results |
| `numerics.py` | Model functions (cost, carbon, safety stock) used by the scripts above |
| `computational_experiments.py`, `chance_experiment.py`, figure scripts | Version 1 illustrations, kept for reference |

## Reproducing the numbers in the paper

Python 3 with NumPy and SciPy:

```bash
python3 revision/calibration_v2.py
python3 revision/application_v2.py
python3 revision/v2_checks.py
```

## Building the paper

```bash
pdflatex network_design_v2 && bibtex network_design_v2 && pdflatex network_design_v2 && pdflatex network_design_v2
```

> Note: copyrighted third-party journal PDFs used as references during the project are excluded from
> this repository.
