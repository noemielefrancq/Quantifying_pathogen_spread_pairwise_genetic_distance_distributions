#########################################################################################################################################
## Spread across countries worldwide - SC2 NCBI data - sensitivity - subsampling
#########################################################################################################################################

#########################################################################################################################################
## Libraries
#########################################################################################################################################
library(ape)
library(stringr)
library(lubridate)
library(Matrix)
library(MetBrewer)

source('1_code/load_pairate.R', chdir = T)
#########################################################################################################################################

#########################################################################################################################################
## Parameters SC2
#########################################################################################################################################
genome_length_SC2 = 29903
substitution_rate_SC2 = 5*10^-4 ## Nextstrain
#########################################################################################################################################

#########################################################################################################################################
## Load full data, distances from tree
#########################################################################################################################################
## Read file precomputed on cluster
dist_mat = readRDS('4_SC2_geography/Data_countries/Seqs/mat_sparse_7_countries.rds')
seq_names = colnames(dist_mat)
n = length(seq_names)
colnames(dist_mat) = rownames(dist_mat) = 1:ncol(dist_mat)
dist_mat@x = dist_mat@x*genome_length_SC2 ## Rescale branch lengths, as iqtree default is N subs per site
dist_mat@x[dist_mat@x > 40] <- 0
dist_mat <- drop0(dist_mat)

# Metadata from NCBI
metadata_raw = read.csv('4_SC2_geography/Data_countries/Seqs/Table_S1_List_isolates_SC2_geo.csv')
metadata_raw = metadata_raw[match(sapply(seq_names,function(x)stringr::str_split(x, '\\.1')[[1]][1]), 
                                  metadata_raw$Accession),]
metadata_raw$id = 1:n

#########################################################################################################################################

#########################################################################################################################################
## Subset data for specific countries
#########################################################################################################################################
countries = c("Denmark", "Japan",  "Slovakia",
              "Switzerland", "Thailand",  "United Kingdom", "USA")
countries_iso2 = c('DK', 'JP', 'SK', 'CH', 'TH', 'UK', 'USA')
countries_colors = met.brewer(name = 'Java', n = length(countries_iso2))

keep = which(!is.na(match(metadata_raw$Country, countries)))
metadata_raw = metadata_raw[keep,]
dist_mat = dist_mat[keep,]
dist_mat = dist_mat[,keep]

## Subset dataset to 200 sequences max per country
max_seqs = 200
idx_remove = NULL
for (ctry in countries) {
  idx_ctry = which(metadata_raw$Country == ctry)
  set.seed(42)
  if(length(idx_ctry) > max_seqs){
    idx_remove = c(idx_remove, sample(idx_ctry, length(idx_ctry) - max_seqs))
  }
}
metadata_raw = metadata_raw[-idx_remove,]
dist_mat = dist_mat[-idx_remove,]
dist_mat = dist_mat[,-idx_remove]
saveRDS(dist_mat, '4_SC2_geography/Data_countries/Seqs/Dist_mat_7_countries_relative_sensitivity.rds')
saveRDS(metadata_raw, '4_SC2_geography/Data_countries/Seqs/Metadata_raw_7_countries_relative_sensitivity.rds')
#########################################################################################################################################

#########################################################################################################################################
## Compute histograms
#########################################################################################################################################
# Named vectors of metadata, matching dist_mat row/colnames                                                                                                                         
traits = setNames(metadata_raw$Country, metadata_raw$id)
times  = setNames(metadata_raw$Collection_Date, metadata_raw$id)                                                                                                                               

bootstrap = 200

# Compute histograms                                                                                                                     
hists = compute_histograms_from_dist(                                                                                                                                               
  dist_mat = dist_mat,                                                                                                                                                              
  traits   = traits,                                                                                                                                                                
  times    = times,                                                                                                                                                                 
  demes    = countries, 
  
  hist_dist_min_max = c(0,15),
  hist_break_method = c('Manual'),
  first_bin_npairs = seq(0,15,1),
  
  collection_time_diff_max = 0.08,
  
  bootstrap = bootstrap
)
countries_pairs_colors = met.brewer(name = 'Signac', n = length(hists$pair_labels))
saveRDS(hists, '4_SC2_geography/Data_countries/Seqs/Hists_7_countries_relative_sensitivity.rds')
######################################

######################################
## Fit migration rates
######################################
fit = fit_hist_model(
  model = 'SYM',
  parallel_chains = 3, chains = 3,
  iter_warmup = 2000, iter_sampling = 2000, 
  refresh = 50,
  
  precomputed_hists = hists, 
  n_bins_to_fit = rep(15, 42),
  stan_data_override = list(rate_dispersion_prior = 1),
  
  bootstrap = 0
)

saveRDS(fit, '4_SC2_geography/Data_countries/Seqs/Fit_7_countries_0_15_2000_2000_relrate_sensitivity.rds')
######################################







