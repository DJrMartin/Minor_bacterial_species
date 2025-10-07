# A conceptual framework for revealing minor bacterial species in microbiome data through guided data transformation

The microbiome is a rich source of biological data that offers promising insights into personalized medicine. However, inferring host health from gut bacterial composition using statistical analytic methods remains a challenge. Here, we show that groups of bacterial species with high abundance and variance (referred to as dominant bacterial species and often associated with enterotypes) exert a disproportionately large influence on microbiome analyses, hiding the contribution of less expressed species (referred to as minor bacterial species). To address this limitation, we propose a guided data transformation highlighting minor bacterial species while minimizing the impact of dominant bacterial species on microbiome statistical analyses. This transformation (i) leads to alternative clustering more closely associated with host health and (ii) helps to improve the performance of supervised machine learning algorithms in high-dimensional settings (𝒏 ≪ 𝒑). Applying to a real dataset, our results suggest that enterotypes may act as a confounding variable to predict host health.

This repository provides data and simulated data, as well as the R scripts used to realize the figures and the experiments.

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

    -   `experiments/` provides the scripts to reproduce the numerical experiments (Figure 3 of the paper).
    -   `real_world/` provides the scripts to reproduce the analyses of the real-world data-sets (Figures 4, 5, S4 and S5).
    -   `simulation_validation/` provides the script to reproduce Figures S6 (i.e., description of the simulated data) and S7 (i.e., augmented data in machine learning context).

<!-- -->

-   `functions/`: provides the script for the out_of_bag procedures for machine learning algorithms. This function is essential to avoid over fitting.

## License

This repository is made available under the [MIT License](LICENSE).

## Citation

If you use this resource in your research, please cite it as:

> [Martin D. et al.] *A conceptual framework for revealing minor bacterial species in microbiome data through guided data transformation*. BioRxiv, [2025].

## Contact

For questions, feedback, or collaboration opportunities, please contact:\
**David Martin (PhD)** -- *University of Rennes*\
[[david.martin.2\@univ-rennes.fr](mailto:david.martin.2@univ-rennes.fr)]
