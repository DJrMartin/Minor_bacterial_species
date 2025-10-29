# A conceptual framework for revealing minor bacterial species in microbiome data through guided data transformation

This repository provides data and simulated data, as well as the R scripts used to realize the figures and the experiments of the guided transformation of microbiome data.

## Installation

git clone <https://github.com/DJrMartin/Minor_bacterial_species.git>

Developed in **R**. Required packages:

``` r
renv::init()
renv::install(readLines("requirements.txt"))

### Import packages.
packages <- readLines("requirements.txt")
invisible(lapply(packages, function(pkg){library(pkg, character.only = TRUE)}))
```

## Repository structure

-   `data/`: real-world data-sets and simulated data from different algorithms.

-   `scr/`:

    -   `Numerical_experiments/` provides the scripts to reproduce the numerical experiments (Figure 3 of the paper).
    -   `Real_world_testing/` provides the scripts to reproduce the analyses of the real-world data-sets (Figures 1, 4, 5, S4).
    -   `Validation_simulation/` provides the script to reproduce Figures S6 (i.e., description of the simulated data) and S7 (i.e., augmented data in machine learning context).

<!-- -->

-   `functions/`: provides the script for the out_of_bag procedures for machine learning algorithms.

## Utilisation

### Load the data

The dataset provided in this repository can be directly used for demonstration or analysis.
For example, Martin et al. (2025) described a metagenomic shotgun dataset composed of 50 samples.

```r
# Load example data
X <- normalised_data   # Feature matrix (e.g., taxonomic counts or abundances)
Y <- target_variables  # Response variable (e.g., fat oxidation levels)
```

### Enterotype definition
The first step of the guided transformation consists of defining enterotypes. It means clustering the samples based on their microbial composition so that each sample belongs to exactly one enterotype.

```r
# Compute Bray–Curtis dissimilarity
bray_curtis <- vegan::vegdist(X, method = "bray")

# Hierarchical clustering and enterotype assignment
Z <- cutree(hclust(bray_curtis, method = "ward.D"), k = 2)

table(Z)  # Check the number of samples per enterotype
```

### Guided transformation

The guided transformation aims to remove the linear effect of the enterotype on each variable,
thus focusing on variation not explained by this grouping factor.

```r
# Compute residuals after adjusting for enterotype
model <- lm(Z ~ X)
X_transf <- residuals(model)
```

### Machine learning for binary classification

You can directly assess the performance of your machine learning models, with or without guided transformation, using the provided out-of-bag testing function.

```r
# Option 1: guided transformation computed on raw data
source("functions/out_of_bag.R")

res <- out_of_bag_prediction(
  as.matrix(X),
  Y,
  cv = 10
)

# Option 2: guided transformation applied on transformed data (e.g., CLR)
res_clr <- out_of_bag_prediction(
  X,
  Y,
  transformed = TRUE,
  X_transformed = compositions::clr(X),
  cv = 10
)
```

Each call to out_of_bag_prediction() returns model performance metrics (e.g., Accuracy, Spe., Sens., or F1) under cross-validation.

This allows for a direct comparison between raw data, guided-transformed data, and other normalization strategies.

## License

This repository is made available under the [MIT License](LICENSE).

## Citation

If you use this resource in your research, please cite it as:

> [Martin D. et al.] *A conceptual framework for revealing minor bacterial species in microbiome data through guided data transformation*. BioRxiv, [2025].

## Contact

For questions, feedback, or collaboration opportunities, please contact:\
**David Martin (PhD)** -- *University of Rennes*\
[[david.martin.2\@univ-rennes.fr](mailto:david.martin.2@univ-rennes.fr)]
