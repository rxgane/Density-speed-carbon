# Cost, Carbon, and the Speed of Replenishment in Distribution Network Design

Paper source and computations for *"Cost, Carbon, and the Speed of Replenishment in Distribution
Network Design"* (Ram Ganeshan, Raymond A. Mason School of Business, William & Mary), version 2,
prepared for *Transportation Research Part E*.

A continuous-approximation location–inventory model in which a firm chooses **network density**
(the number of distribution centers, `n`) and **replenishment speed** (`v`) to minimize cost under a
carbon price or cap and a service requirement. Carbon arises from facilities, holding inventory,
line-haul freight, and last-mile delivery; cost also includes the financing of stock in transit. The joint
design has a unique optimum at every carbon price, and the carbon-minimizing direction of both decisions
follows from four carbon-intensity-per-dollar ratios, the share of a distribution center's marginal cost that
is lost pooling, and the share of the inventory saving from speed that is safety stock. A calibration to
public US data traces how a carbon price moves a reference network that finances its stock in transit.

## Repository contents

| File | Description |
|------|-------------|
| `network_design_v2.tex`, `references.bib` | Manuscript source (anonymized for review) and bibliography |
| `network_design_v2.pdf` | Compiled manuscript |
| `title_page.tex`, `cover_letter.tex`, `highlights.tex` | Title page with declarations, cover letter, highlights |
| `revision/application_v21.py` | Benchmark carbon intensities per dollar (Table 3), allocation shares, carbon shares of priced cost |
| `revision/application_v22.py` | Reference network that finances stock in transit: response to a carbon price, line-haul link, robustness, rounding |
| `revision/fig_price.py` | Figure of the reference network's response to a carbon price (`fig_price.pdf`) |
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
cd revision
python3 application_v21.py
python3 application_v22.py
python3 fig_price.py
```

## Building the paper

```bash
pdflatex network_design_v2 && bibtex network_design_v2 && pdflatex network_design_v2 && pdflatex network_design_v2
```

> Note: copyrighted third-party journal PDFs used as references during the project are excluded from
> this repository.
