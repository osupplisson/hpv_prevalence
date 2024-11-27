# General description
This repository contains the codes for the study 'High-risk human papillomavirus cervical infection prevalence in France, 2020-2023: a nationwide, large-scale, and spatially resolved retrospective study comparing opportunistic and organised screening'. 

# History

## Peer-review/Preprint

* 21/10/2024: first public version posted on MedRxiv (https://www.medrxiv.org/content/10.1101/2024.10.20.24315479v1, doi: https://doi.org/10.1101/2024.10.20.24315479).
* 21/10/2024: submitted to Eurosurveillance
* 29/10/2024: reviewers invited
* 08/11/2024: under review


## Typo correction

Several typos were identified after the document was released on MedRxiv. The next medrxiv update will be released following the first peer review round. In the meantime, a typo-corrected version can be found on my personal Google drive, following this url: https://drive.google.com/drive/folders/1on-osm5h_M3CVgz2mXA0YuDsN1kFgKJC?usp=drive_link

## Major updates

* 19/11/2025: Additional neighbourhood matrix (2nd order Queen) considered. Replacement of the single-level prior of each of the three additional intercepts by a hierarchical prior (PC-prior). Results remain the same.
* 27/11/2025: Fixing an issues with the reported ETI95 (there were defined as $[q_{a},q_{1-a}]$ instead of $[q_{a/2},q_{1-a/2}]$). Adding more draws from the joint posterior to compute the quantities (2,000 for PP and 1,200 for other quantities).


# Description of codes

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

Note: The code `paper.Rmd` will be posted once the manuscript is accepted for publication. It contains the `.Rmd` file that produced the whole manuscript (main text and supplementary files).

## Funding

* OS (1st author) is funded by a Ph.D. grant from the ANRS|Maladies infectieuses émergentes (grant number 22485).
* NT (2nd author) is funded by a research grant from the ANRS|Maladies infectieuses émergentes (grant number 21290).
* This project was also funded by the ANRS|Maladies infectieuses émergentes through the research project ANRS-0681-MIST.

## Data sharing

Important: Data cannot be shared due to legal constraints associated with using the French reference methodology MR-004. 


## Acknowledgements

We are grateful to the genotoul bioinformatics platform Toulouse Occitanie (Bioinfo Genotoul, https://doi.org/10.15454/1.5572369328961167E12) for providing help and computing and storage resources.

The authors acknowledge the  ISO 9001 certified IRD i-Trop HPC (South Green Platform) at IRD montpellier for providing HPC resources that have contributed to the research results reported within this paper. https://bioinfo.ird.fr/-http://www.southgreen.fr

We acknowledge Bioclust, the computing cluster of PB-IBENS (LABEX MEMOLIFE members) for providing storage resources. https://www.ibens.bio.ens.psl.eu/spip.php?rubrique55#bioclust

