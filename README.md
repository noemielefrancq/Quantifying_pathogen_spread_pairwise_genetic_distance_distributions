# Quantifying pathogen spread across scales from pairwise genetic distance distributions

Code of the paper "Quantifying pathogen spread across scales from pairwise genetic distance distributions".

The repo contains the code of the method (pairate) and the scripts to reproduce the analyses of the paper.

## Structure

- `1_code`: pairate, the R code of the method (functions and Stan model). Load the script with `source('1_code/load_pairate.R', chdir = T)`.
- `2_theory`: theoretical figures (approximations and likelihood profiles).
- `3_simulations_ReMaster`: ReMaster XML templates used to simulate phylogenies, and example simulated trees.
- `4_SC2_geography`: spread of SARS-CoV-2 between countries.
- `5_age_groups`: spread of SARS-CoV-2, influenza A/H1N1pdm, A/H3N2, RSV-A and RSV-B across age groups.
- `6_MERS`: spread of MERS-CoV between camels and humans.

## Running the analyses

SARS-CoV-2 between countries, in this order: `Construct_matrix_mobility.R` (optional), `SC2_geo_run_model.R`, `SC2_geo_sensitivity.R` and `SC2_geo_figures.R`.

Age groups: `Age_group_run_model.R`, then `Age_group_figures.R`. This analysis requiers GISAID data (see below).

MERS-CoV: `MERS_analysis.R`. The scripts in `2_theory` can be run independently.

## Data

SARS-CoV-2 between countries: sequences from NCBI. The list of isolates (Table S1), the phylogeny and the distance matrix are included. The distance matrix only contains the pairs used in the analysis (up to 40 SNPs apart, sampled within 0.08 years of each other). Mobility data are from GTMD2 and can be downloaded from https://zenodo.org/records/18496028 (GTMD2_trips.csv). Population sizes are from the World Bank (2022) and can be downloaded from https://data.worldbank.org/indicator/SP.POP.TOTL.

MERS-CoV: sequences from Dudas et al. (2018) and BEAST2 analysis from https://github.com/tgvaughan/MultiTypeTrajectoryAnalyses/tree/master/mers (included).

Age groups: contact matrices from SocialCov (can be downloaded from https://zenodo.org/records/13834596) and COMES-F (downloaded with socialmixr), and population sizes from INSEE (can be downloaded from https://catalogue-donnees.insee.fr/fr/catalogue/recherche/DS_RP_TD_POPULATION_AGESEX_PRINC). Sequences and metadata are from GISAID and cannot be shared. They can be retrieved with the following EPI_SET IDs: EPI_SET_261110ge (SARS-CoV-2), EPI_SET_261110qd (H1N1), EPI_SET_261110xe (H3N2), EPI_SET_261110mh (RSV-A) and EPI_SET_261110rs (RSV-B).
