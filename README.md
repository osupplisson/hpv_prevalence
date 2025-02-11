# General description
This repository contains the codes for the study 'High-risk human papillomavirus cervical infection prevalence in France, 2020-2023: a nationwide, large-scale, and spatially resolved retrospective study comparing opportunistic and organised screening'. 

# History

## Peer-review/Preprint

* 21/10/2024: first public version posted on MedRxiv (https://www.medrxiv.org/content/10.1101/2024.10.20.24315479v1, doi: https://doi.org/10.1101/2024.10.20.24315479).
* 21/10/2024: submitted to Eurosurveillance
* 29/10/2024: reviewers invited
* 08/11/2024: under review
* 14/01/2025: reviews received, revised version due to 07/02/2025
* 09/02/2025: revision submitted, preprint updated (https://www.medrxiv.org/content/10.1101/2024.10.20.24315479v2), .RMD uploaded on the github page.


## Major updates

* 19/11/2025: Additional neighbourhood matrix (2nd order Queen) considered. Replacement of the single-level prior of each of the three additional intercepts by a hierarchical prior (PC-prior). Results remain the same.
* 01/12/2025: Fixing an issues with the reported ETI95 (before 01/12/2025 there were reported as $[q_{a},q_{1-a}]$ instead of the correct way, $[q_{a/2},q_{1-a/2}]$). Adding more draws from the joint posterior to compute the quantities (2,000 for PP and 1,500 for other quantities). Percentiles reported in supplementary files are now the 50th, 90th, 97.5th, and 99th quantiles (and the complement to 1 of these percentiles).
* 21/01/2025: Run with the latest R-INLA version (24.12.11 version, see https://github.com/hrue/r-inla/blob/devel/rinla/NEWS.md) and update of the build used on the HPC (Rocky Linux-8.10 (Green Obsidian) instead of Centos 7 in previous versions).
* 09/02/2025: Run with the 25.01.23 R-INLA version. All quantities but the Average Marginal Effect ('MDEP') are now computed based on 3,000 draws from the joint posterior distribution (1,550 for the MDEP). Also, *all* quantities for a given model are now computed from the *same* draws from the joint distribution.

# Description of codes

The `sources.bib` contains the `.bib` information used by the `.Rmd` files for citations (Nb: not all entries are used).

In the `code` folder you will find the following R code:
1. `00_package.R`: Import needed packages and define some functions and objects used in subsequent codes;
2. `01_mask.R`: Create the boundary for the spatial domain. Notably removing small islands around France mainland;
3. `02_spatial.R`: Import various datasets containing spatial information;
4. `03_import.R`: Import the data, subpopulation selection (see flowchart);
5. `04_descriptive_statistics.R`: Descriptive statistics figures for the Descriptive sample;
6. `05_inlabru.R`: Keeping analytic sample, fit all models, computed LGO-CV log-score, perform PP-check, perform post-fit analyses, select model, fit sensitivity analyses for fitted model;
7. `99_post_fit_function.R`: Functions used for PP check and post-fit analyses;
8. `paper.Rmd`: contains the .RMD files that produced the complete .PDF manuscript;


Any suggestions, comments, or error reporting are welcomed! Please, reach out to `osupplis@gmail.com`.


## Funding

* OS (1st author) is funded by a Ph.D. grant from the ANRS|Maladies infectieuses émergentes (grant number 22485).
* NT (2nd author) is funded by a research grant from the ANRS|Maladies infectieuses émergentes (grant number 21290). He is the recipient of a Fellowship from the ExposUM Institute. The grant is a combined contribution from the Université de Montpellier, Région d'Occitanie, and France 2030.
* This project was also funded by the ANRS|Maladies infectieuses émergentes through the research project ANRS-0681-MIST.

## Data sharing

Important: Data cannot be shared due to legal constraints associated with using the French reference methodology MR-004. 


## Acknowledgements

We are grateful to the genotoul bioinformatics platform Toulouse Occitanie (Bioinfo Genotoul, https://doi.org/10.15454/1.5572369328961167E12) for providing help and computing and storage resources.

The authors acknowledge the  ISO 9001 certified IRD i-Trop HPC (South Green Platform) at IRD montpellier for providing HPC resources that have contributed to the research results reported within this paper. https://bioinfo.ird.fr/-http://www.southgreen.fr

We acknowledge Bioclust, the computing cluster of PB-IBENS (LABEX MEMOLIFE members) for providing storage resources. https://www.ibens.bio.ens.psl.eu/spip.php?rubrique55#bioclust

