# DEEPOL-Marine-Biogeography

**Biogeography and Distribution Drivers of Three Novel Marine Microbe Lineages**

This repository contains a modular analytical pipeline for the biogeographic study of three newly identified marine microbial lineages (**DEEPOL1, DEEPOL2, and DEEPOL3**) using large-scale metabarcoding datasets.

## Project Overview
This study explores the ecological niches, community assembly, and phylogenetic specialization of these lineages across global ocean basins. The workflow is designed to be **modular**, allowing for independent configuration and execution of specific analyses for each lineage.

## Key Analytical Modules
The pipeline (`DEEPOLs_modular.Rmd`) is organized into the following components:

* **Data Integration:** Merging ASV presence/absence matrices with global environmental metadata and size fraction data.
* **Community Structure:** * **UpSet Plots:** Visualizing ASV intersections across size fractions.
    * **PCoA (Jaccard):** Visualizing community similarity with environmental vector fitting (`envfit`).
* **Drivers of Distribution:** * **dbRDA:** Distance-based Redundancy Analysis to identify specific environmental correlates.
    * **Variation Partitioning:** Decoupling the relative influence of Environment vs. Geography (spatial coordinates).
* **Statistical Significance:** * **PERMANOVA:** Testing for significant differences between ocean basins and layers.
    * **PERMDISP:** Assessing multivariate dispersion and centroid distances to validate community shifts.
* **Phylogenetic Niche Traits:** * Mapping habitat specialization (Polarity, Ocean Layer) and niche breadth (Occupancy) onto phylogenetic trees using `ggtree`.

## Repository Structure
* `DEEPOLs_modular.Rmd`: The primary R Markdown analysis notebook.
* `data/`: (Expected) Contains lineage-specific filtered abundance matrices and metadata.
* `results/`: Automatically generated output directories containing PDFs and CSVs for each lineage.
* `scripts/`: Contains helper functions for niche trait calculations and reporting.

## Requirements
The analysis is performed in **R** and utilizes the following library suite:
* **Data Wrangling:** `tidyverse`, `tibble`, `reshape2`, `here`, `conflicted`
* **Statistics:** `vegan`, `ape`
* **Visualization:** `ggplot2`, `patchwork`, `ggforce`, `ggtree`, `ggnewscale`

