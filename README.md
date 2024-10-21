# Description
This repository contains the codes for the study 'High-risk human papillomavirus cervical infection prevalence in France, 2020-2023: a nationwide, large-scale, and spatially-resolved study comparing opportunistic and organised screening'

Important: Data cannot be shared due to legal constraints associated with using the French reference methodology MR-004. 

The `paper.Rmd` contains the `.Rmd` file that *exactly* produced the manuscript (main text and supplementary files) posted on MedRxiv (DOI to be updated once accepted).
The `sources.bib` contains the `.bib` informations used by the `.Rmd` files for citations (Nb: not all entries are used).

In the `code` folder you will find the following R code:
1. `00_package.R`: Import needed packages and define some functions and objects used in subsequent codes;
2. `01_mask.R`: Create the boundary for the spatial domain. Notably removing small islands around France mainland;
3. `02_spatial.R`: Import various datasets containing spatial information;
4. `03_import.R`: Import the data, subpopulation selection (see flowchart);
5. `04_descriptive_statistics.R`: Descriptive statistics figures for the Descriptive sample;
6. `05_inlabru.R`: Keeping analytic sample, fit all models, computed LGO-CV log-score, perform PP-check, perform post-fit analyses, select model, fit sensitivity analyses for fitted model;
7. `99_post_fit_function.R`: Functions used for PP check and post-fit analyses


Any suggestions, comments, or error reporting are welcomed! Please, reach out to `osupplis@gmail.com`.

## Funding

* OS (1st author) is funded by a Ph.D. grant from the ANRS|Maladies infectieuses émergentes (grant number 22485).
* NT (2nd author) is funded by a research grant from the ANRS|Maladies infectieuses émergentes (grant number 21290).
* This project was also funded by the ANRS|Maladies infectieuses émergentes through the research project ANRS-0681-MIST.

## Acknowledgements

We are grateful to the genotoul bioinformatics platform Toulouse Occitanie (Bioinfo Genotoul, https://doi.org/10.15454/1.5572369328961167E12) for providing help and computing and storage resources".

The authors acknowledge the  ISO 9001 certified IRD i-Trop HPC (South Green Platform) at IRD montpellier for providing HPC resources that have contributed to the research results reported within this paper. https://bioinfo.ird.fr/-http://www.southgreen.fr

We acknowledge Bioclust, the computing cluster of PB-IBENS (LABEX MEMOLIFE members) for providing storage resources. https://www.ibens.bio.ens.psl.eu/spip.php?rubrique55#bioclust
