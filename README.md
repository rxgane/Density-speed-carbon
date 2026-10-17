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
| `revision/application_v21.py` | Benchmark carbon intensities per dollar (Table 2), with warehouse carbon at the EPA input–output rate; carbon shares of priced cost |
| `revision/application_v22.py` | Reference network that finances stock in transit: inputs, response to a carbon price, rounding |
| `revision/targets_v22.py` | Least-cost design for each required cut in emissions, three action sets; carbon-price path |
| `revision/levers_data.py`, `revision/fig_levers.py` | Cost and carbon as each decision moves on its own (Figure 2, `fig_levers.pdf`) |
| `revision/fig_target.py` | Least-cost design for each required cut (Figure 3, `fig_target.pdf`) |
| `revision/inbound_data.py`, `revision/inbound_fine.py`, `revision/fig_inbound.py` | Carbon intensity of faster inbound service and the number of distribution centers (Figure 4, `fig_inbound.pdf`) |
| `revision/robust_v22.py` | Robustness of the least-cost design for a 30% cut (Table 4) |
| `revision/sensitivity_v24.py` | Sensitivity of that design to the carbon inputs and the storage environment (Table 5) |
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
python3 calibration_v2.py
python3 targets_v22.py && python3 levers_data.py && python3 fig_levers.py && python3 fig_target.py
python3 inbound_data.py && python3 inbound_fine.py && python3 fig_inbound.py
python3 robust_v22.py
python3 sensitivity_v24.py > sensitivity_v24_output.txt
```

## Building the paper

```bash
pdflatex network_design_v2 && bibtex network_design_v2 && pdflatex network_design_v2 && pdflatex network_design_v2
```

> Note: copyrighted third-party journal PDFs used as references during the project are excluded from
> this repository.
