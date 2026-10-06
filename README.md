# ParSlet: Prevalence-Aware Feature Selection Improves Biomarker Identification in Microbiome Studies

This repository contains R code associated with our published work on prevalence-aware feature selection for microbiome biomarker discovery.

## Publication

**Yang R., Li Y., Sankaran K., Mace T.A., Hart P.A., Ma Q., Wang X.-W., Ke S.**  
*Prevalence aware feature selection improves biomarker identification in microbiome studies.*  
**Bioinformatics**, 2026, 42(7), btag371.

**DOI:** https://doi.org/10.1093/bioinformatics/btag371

## Overview

Machine-learning approaches are widely used for microbiome biomarker discovery, but selected biomarkers can vary substantially across datasets and sample sizes.

**ParSlet** is a prevalence-aware feature selection framework that integrates microbial prevalence with Random Forest feature importance to improve the stability and reproducibility of biomarker identification.

The framework was evaluated using simulated microbiome data and multiple real-world microbiome datasets.

## Core Idea

We observed a relationship between:

- taxon prevalence $p$
- frequency of being selected as a biomarker $f$

This relationship approximately follows a power-law scaling:

$$
f(p) \propto p^{\alpha}
$$

where $\alpha \approx 2.5$ in both simulated and real microbiome datasets.

Based on this observation, the original Random Forest feature importance is adjusted using taxon prevalence:

$$
I_j^{\mathrm{adj}} = I_j \times p_j^{\alpha}
$$

where:

- $I_j$: original feature importance (Gini or MDA)
- $p_j$: prevalence of taxon $j$
- $\alpha$: prevalence exponent (default = 2.5)

The resulting **prevalence-integrated importance score** is used to rank candidate microbial biomarkers.

## Method Workflow

The general analysis workflow is:

1. Prepare the microbiome abundance table.
2. Split data into training and testing sets.
3. Compute taxon prevalence using the training data.
4. Train a Random Forest model.
5. Extract feature importance scores.
6. Adjust feature importance using taxon prevalence.
7. Rank features based on prevalence-integrated importance.
8. Select top-ranked candidate biomarkers.
9. Evaluate feature-selection stability and predictive performance.

## Installation

Clone this repository:

```bash
git clone https://github.com/yruoxi00/ParSlet-Feature-Selection.git
cd ParSlet-Feature-Selection
```

## Required R Packages

The analyses were conducted in R. Major packages used in this project include:

```r
install.packages(c(
  "randomForest",
  "dplyr",
  "ggplot2"
))
```

Simulation analyses additionally use **SparseDOSSA2**.

## Repository Contents

This repository contains R scripts for the main analyses presented in the study, including:

- Microbiome data preprocessing
- Random Forest modeling
- Taxon prevalence calculation
- Prevalence-aware feature selection
- Simulation experiments
- Feature-selection stability evaluation
- Statistical analysis and visualization

## Data Availability

The analyses use simulated data and publicly available microbiome datasets.

Real microbiome datasets were obtained from resources including:

- **MDAD (Microbiome Differential Abundance Datasets)**
- **curatedMetagenomicData**

Simulation data were generated using **SparseDOSSA2**.

## Citation

If you use this method or code, please cite:

> Yang R, Li Y, Sankaran K, Mace TA, Hart PA, Ma Q, Wang X-W, Ke S.  
> **Prevalence aware feature selection improves biomarker identification in microbiome studies.**  
> *Bioinformatics*. 2026;42(7):btag371.  
> https://doi.org/10.1093/bioinformatics/btag371