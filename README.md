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
