# Predicting VRC01 Neutralization of HIV-1 with Gaussian Process Regression

Exploratory machine-learning work I did as an ML Researcher at the Broad Institute of MIT and Harvard.

VRC01 is a broadly neutralizing antibody against HIV-1. Its potency against an HIV-1 viral variant is measured as an IC50 (lower IC50 means more potent). Some viral variants are unable to be neutralized at all, which is recorded as undetectable (`UD`). This project predicts:

1. whether VRC01 will neutralize it at all (UD vs. not UD), and
2. if so, its log IC50, along with an uncertainty estimate for each prediction.

from a virus's aligned Env protein sequence.

## Data

`data/vrc01-catnap-data.csv` contains VRC01 neutralization measurements paired with aligned HIV-1 Env sequences (1,061 alignment positions), drawn from the Los Alamos CATNAP database.

| | Sequences |
|---|---|
| Rows in the file | 2,096 |
| After dropping rows with no sequence | 1,926 |
| With a numeric IC50 (used for regression) | 1,663 |
| Undetectable (`UD`) | 263 (13.7%) |

IC50 values span 0.001 to 100 µg/mL, so they are modeled on a log10 scale (−3 to 2).

## Why a Gaussian process

About 1,900 sequences, each with roughly 1,000 positions, is a high-dimensional, low-sample problem. The IC50 values are also noisy, and the UD values are effectively infinite outliers. A Gaussian process suits this setting because it:

- has only a few parameters to fit, instead of thousands of network weights;
- predicts by similarity to known sequences, with similarity defined by a kernel I control; and
- returns an uncertainty with every prediction, so it can flag a virus that looks unlike anything in the training data.

## Method

**Sequence kernel.** For two aligned sequences, the similarity is the fraction of matching residues, counted only over positions where at least one sequence is not a gap. The kernel is vectorized in NumPy, so the full 1,330 × 1,330 training matrix computes in about 8 seconds.

**Step 1: UD classifier.** A class-balanced logistic regression on the kernel matrix predicts whether a sequence is neutralizable. Because only 13.7% of sequences are UD, I swept the decision threshold to trade UD recall against false positives and settled on 0.4.

**Step 2: GP regression.** For sequences predicted to be neutralizable, a Gaussian process using the same kernel (with fixed noise of 0.01 added to the diagonal) predicts log IC50 and its standard deviation.

**Pipeline.** `predict_ic50()` runs both steps on new sequences and returns the classification, the UD probability, the predicted IC50 and the uncertainty.

## Results

Held-out 20% test sets, random split (`random_state=42`).

**Regression** (333 test sequences)

| Model | Pearson r |
|---|---|
| scikit-learn `GaussianProcessRegressor`, precomputed kernel, no noise term | 0.43 |
| GP with noise term added to the kernel diagonal | **0.56** (p ≈ 1e-28) |

The uncertainty estimates are informative. Test predictions with below-median uncertainty had a mean absolute error of 0.32 log10 units, and those with above-median uncertainty had 0.59, nearly twice as large.

**UD classification** (386 test sequences, 53 UD)

| Model | AUC-ROC | UD recall | UD precision |
|---|---|---|---|
| Logistic regression, unweighted | 0.77 | 0.04 | 1.00 |
| Logistic regression, class-balanced, threshold 0.5 | 0.79 | 0.64 | 0.28 |
| Logistic regression, class-balanced, threshold 0.4 | 0.79 | 0.81 | 0.23 |

Without class balancing, the classifier labels almost everything as neutralizable. Balancing and lowering the threshold recover most UD viruses, at the cost of more false positives.

## Limitations and next steps

- **One antibody.** All results are for VRC01. Other antibodies would need their own models or a multi-task setup.
- **Random split.** Closely related viruses can appear in both training and test sets, which likely makes these scores optimistic. About 20 sequences in the file are exact duplicates. Deduplicating and splitting by clade or sequence-similarity cluster would give a harder, more realistic test.
- **Fixed hyperparameters.** The noise level and kernel scale were set by hand rather than fit by maximum likelihood.
- **Simple kernel.** Identity matching treats every substitution as equally different. Planned alternatives: weighting substitutions by amino-acid similarity (e.g., BLOSUM), or using cosine similarity of ESM-2 embeddings.
- **Escape mutations.** Two sequences can be highly similar overall yet differ at a single residue that lets the virus escape VRC01. Combining the kernel with known escape positions is a natural extension.

## Repository structure

```
gpr-ic50-prediction/
  data/vrc01-catnap-data.csv          IC50 values + aligned Env sequences
  notebooks/gpr-ic50-prediction.ipynb full analysis and pipeline
```

## Running it

```bash
pip install numpy pandas scipy scikit-learn matplotlib jupyter
jupyter notebook gpr-ic50-prediction/notebooks/gpr-ic50-prediction.ipynb
```

The notebook loads `vrc01-catnap-data.csv` from its working directory. Either run it with the data file alongside it or change the path in the "Prepare Data" cell to `../data/vrc01-catnap-data.csv`.
