#########################################################################################################################################
## Age group fit, for different pathogens in ARA
#########################################################################################################################################
library(ape)
library(stringr)
library(lubridate)

source('1_code/load_pairate.R', chdir = T)

########################################
## This file needs GISAID DATA TO RUN ##
########################################


#########################################################################################################################################
## Load data
#########################################################################################################################################
Months = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec')

######################################
## SC2 (EPI_SET_261110ge)
######################################
metadata_SC2 = read.csv('5_age_groups/Data/SC2/EPI_SET_261110ge_SC2.tsv', sep = '\t')
metadata_SC2 = metadata_SC2[which(unlist(lapply(metadata_SC2$Location, function(x)str_split(x, ' / ')[[1]][3])) == "Auvergne-Rhone-Alpes"),]
metadata_SC2$Collection.date = decimal_date(as_date(metadata_SC2$Collection.date))
metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == '1M')] = 1/12
tmp = which(unlist(lapply(metadata_SC2$Patient.age, function(x)str_split(x, ' ')[[1]][2])) == "years")
for(i in 1:length(tmp)){
  metadata_SC2$Patient.age[tmp[i]] = str_split(metadata_SC2$Patient.age[tmp[i]], ' years')[[1]][1]
}
tmp = which(unlist(lapply(metadata_SC2$Patient.age, function(x)str_split(x, ' ')[[1]][2])) == "years,")
for(i in 1:length(tmp)){
  metadata_SC2$Patient.age[tmp[i]] = str_split(metadata_SC2$Patient.age[tmp[i]], ' years,')[[1]][1]
}
tmp = which(unlist(lapply(metadata_SC2$Patient.age, function(x)str_split(x, ' ')[[1]][2])) == "months")
for(i in 1:length(tmp)){
  metadata_SC2$Patient.age[tmp[i]] = as.character(as.numeric(str_split(metadata_SC2$Patient.age[tmp[i]], ' months')[[1]][1])/12)
}
metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == 'unknown')] = NA
metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == '1 month')] = 1/12
metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == '<2')] = 1
metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == '2-17')] = 10

metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == '18-64')] = NA
metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == '>64')] = '>65'
metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == '>=65')] = '>65'
metadata_SC2$Patient.age[which(metadata_SC2$Patient.age == '>65')] = NA
metadata_SC2$Patient.age = as.numeric(metadata_SC2$Patient.age)

metadata_SC2 = metadata_SC2[which(metadata_SC2$Collection.date <= 2025.5),]
metadata_SC2$surveillance_system = NA
idx = (which(!is.na(sapply(metadata_SC2$Virus.name, function(x)str_split(x, 'HCL')[[1]][2]))))
metadata_SC2$surveillance_system[idx] = 'HCL'
idx = (which(!is.na(sapply(metadata_SC2$Virus.name, function(x)str_split(x, 'RELAB')[[1]][2]))))
metadata_SC2$surveillance_system[idx] = 'RELAB'
tree_SC2 = read.tree('5_age_groups/Data/SC2/EPI_SET_261110ge_SC2.treefile')
tree_SC2 = root(phy = tree_SC2, outgroup = 'hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31', resolve.root = T) 
genome_length_SC2 = 29903
substitution_rate_SC2 = 5*10^-4 ## Nextstrain
tree_SC2$edge.length = tree_SC2$edge.length*genome_length_SC2
idx = match(unlist(lapply(tree_SC2$tip.label, function(x)str_split(x, '\\|')[[1]][2])), metadata_SC2$Accession.ID)
tree_SC2 = drop.tip(tree_SC2, tree_SC2$tip.label[which(is.na(idx))])
######################################
## H1N1 (EPI_SET_261110qd)
######################################
metadata_h1n1 = readxl::read_excel('5_age_groups/Data/H1N1/EPI_SET_261110qd_H1N1.xls')
metadata_h1n1 = as.data.frame(metadata_h1n1)
metadata_h1n1 = metadata_h1n1[which(metadata_h1n1$Location == "Europe / France / Auvergne-Rhone-Alpes" | 
                                      metadata_h1n1$Location =="Europe / France / Rhone-Alpes" |
                                      metadata_h1n1$Location =="Europe / France / Auvergne"),]
metadata_h1n1$Collection_Date = decimal_date(as_date(metadata_h1n1$Collection_Date))
metadata_h1n1 = metadata_h1n1[which(metadata_h1n1$Collection_Date <= 2025.5),]
tree_h1n1 = read.tree('5_age_groups/Data/H1N1/EPI_SET_261110qd_H1N1.treefile')
tree_h1n1 = root(phy = tree_h1n1, outgroup = 'EPI_ISL_404527|A/Wisconsin/588/2019|2019-12-19|HA|A_/_H1N1|EPI1661758', resolve.root = T) 
genome_length_h1n1 = 1752
substitution_rate_h1n1 = 3.3*10^-3 ## https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0093632
tree_h1n1$edge.length = tree_h1n1$edge.length*genome_length_h1n1
idx = match(unlist(lapply(tree_h1n1$tip.label, function(x)str_split(x, '\\|')[[1]][1])), metadata_h1n1$Isolate_Id)
tree_h1n1 = drop.tip(tree_h1n1, tree_h1n1$tip.label[which(is.na(idx))])

######################################
## H3N2 (EPI_SET_261110xe)
######################################
metadata_h3n2 = readxl::read_excel('5_age_groups/Data/H3N2/EPI_SET_261110xe_H3N2.xls')
metadata_h3n2 = as.data.frame(metadata_h3n2)
metadata_h3n2 = metadata_h3n2[which(metadata_h3n2$Location == "Europe / France / Auvergne-Rhone-Alpes" | 
                                      metadata_h3n2$Location =="Europe / France / Rhone-Alpes" |
                                      metadata_h3n2$Location =="Europe / France / Auvergne"),]
metadata_h3n2$Collection_Date = decimal_date(as_date(metadata_h3n2$Collection_Date))
metadata_h3n2 = metadata_h3n2[which(metadata_h3n2$Collection_Date <= 2025.5),]
tree_h3n2 = read.tree('5_age_groups/Data/H3N2/EPI_SET_261110xe_H3N2.treefile')
tree_h3n2 = root(phy = tree_h3n2, outgroup = 'EPI_ISL_1563628|A/Darwin/6/2021|2021-03-16|HA|A_/_H3N2|EPI1857216', resolve.root = T) 
genome_length_h3n2 = 1718
substitution_rate_h3n2 = 4*10^-3 #https://www.pnas.org/doi/pdf/10.1073/pnas.94.15.7712
tree_h3n2$edge.length = tree_h3n2$edge.length*genome_length_h3n2
idx = match(unlist(lapply(tree_h3n2$tip.label, function(x)str_split(x, '\\|')[[1]][1])), metadata_h3n2$Isolate_Id)
tree_h3n2 = drop.tip(tree_h3n2, tree_h3n2$tip.label[which(is.na(idx))])

######################################
## RSVA (EPI_SET_261110mh)
######################################
metadata_RSVA = read.csv('5_age_groups/Data/RSVA/EPI_SET_261110mh_RSVA.tsv', sep = '\t')
metadata_RSVA = metadata_RSVA[which(metadata_RSVA$Location == "Europe / France / Auvergne-Rhone-Alpes" | 
                                      metadata_RSVA$Location =="Europe / France / Rhone-Alpes" |
                                      metadata_RSVA$Location =="Europe / France / Auvergne"),]
metadata_RSVA$Collection.date = decimal_date(as_date(metadata_RSVA$Collection.date))
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '1M')] = 1/12
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '1 month')] = 1/12
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '2 months')] = 2/12
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '3 months')] = 3/12
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '4 months')] = 4/12
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '5 months')] = 5/12
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '6 months')] = 6/12
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '2-17')] = 10
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '18-64')] = 40
metadata_RSVA$Patient.age[which(metadata_RSVA$Patient.age == '>64')] = 65
metadata_RSVA$Patient.age = as.numeric(metadata_RSVA$Patient.age)
metadata_RSVA = metadata_RSVA[which(metadata_RSVA$Collection.date <= 2025.5),]
tree_RSVA = read.tree('5_age_groups/Data/RSVA/EPI_SET_261110mh_RSVA.treefile')
tree_RSVA = root(phy = tree_RSVA, outgroup = 'hRSV/A/England/397/2017|EPI_ISL_412866|2017|2017', resolve.root = T) 
genome_length_RSVA = 15225
substitution_rate_RSVA = 6.7*10^-4 # https://pmc.ncbi.nlm.nih.gov/articles/PMC8217232/
tree_RSVA$edge.length = tree_RSVA$edge.length*genome_length_RSVA
idx = match(unlist(lapply(tree_RSVA$tip.label, function(x)str_split(x, '\\|')[[1]][2])), metadata_RSVA$Accession.ID)
tree_RSVA = drop.tip(tree_RSVA, tree_RSVA$tip.label[which(is.na(idx))])

######################################
## RSVB (EPI_SET_261110rs)
######################################
metadata_RSVB = read.csv('5_age_groups/Data/RSVB/EPI_SET_261110rs_RSVB.tsv', sep = '\t')
metadata_RSVB = metadata_RSVB[which(metadata_RSVB$Location == "Europe / France / Auvergne-Rhone-Alpes" | 
                                      metadata_RSVB$Location =="Europe / France / Rhone-Alpes" |
                                      metadata_RSVB$Location =="Europe / France / Auvergne"),]
metadata_RSVB$Collection.date = decimal_date(as_date(metadata_RSVB$Collection.date))
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '1M')] = 1/12
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '1 month')] = 1/12
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '2 months')] = 2/12
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '3 months')] = 3/12
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '4 months')] = 4/12
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '5 months')] = 5/12
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '6 months')] = 6/12
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '2-17')] = 10
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '18-64')] = 40
metadata_RSVB$Patient.age[which(metadata_RSVB$Patient.age == '>64')] = 65
metadata_RSVB$Patient.age = as.numeric(metadata_RSVB$Patient.age)
metadata_RSVB = metadata_RSVB[which(metadata_RSVB$Collection.date <= 2025.5),]
tree_RSVB = read.tree('5_age_groups/Data/RSVB/EPI_SET_261110rs_RSVB.treefile')
tree_RSVB = root(phy = tree_RSVB, outgroup = 'hRSV/B/Australia/VIC-RCH056/2019|EPI_ISL_1653999|2019|2019', resolve.root = T) 
genome_length_RSVB = 15222
substitution_rate_RSVB = 7.7*10^-4
tree_RSVB$edge.length = tree_RSVB$edge.length*genome_length_RSVB
idx = match(unlist(lapply(tree_RSVB$tip.label, function(x)str_split(x, '\\|')[[1]][2])), metadata_RSVB$Accession.ID)
tree_RSVB = drop.tip(tree_RSVB, tree_RSVB$tip.label[which(is.na(idx))])
######################################

#########################################################################################################################################
## Compute ratios by age, for each pathogen
#########################################################################################################################################
traits_SC2 = metadata_SC2$Patient.age[match(unlist(lapply(tree_SC2$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_SC2$Accession.ID)]
age_breaks_SC2 = c(0,20,40,60,80,120)
traits_SC2 = cut(as.numeric(traits_SC2), age_breaks_SC2)
times_SC2 = metadata_SC2$Collection.date[match(unlist(lapply(tree_SC2$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_SC2$Accession.ID)]
demes_SC2 = sort(unique(traits_SC2))
median.and.ci = function(v){ return( c(median(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))}
hists_SC2 = compute_histograms(tree = tree_SC2,
                               traits = traits_SC2,
                               times = times_SC2,
                               demes = demes_SC2, 
                               
                               hist_dist_min_max = c(0,5),
                               hist_break_method = c('Manual'),
                               first_bin_npairs = seq(0,5,1),
                               max_n_bins = 1000,
                               
                               collection_time_diff_max = 0.08, ## 1 month interval (~4 generations)

                               bootstrap = 200)

traits_h1n1 = metadata_h1n1$Host_Age[match(unlist(lapply(tree_h1n1$tip.label, function(x)str_split(x,'\\|')[[1]][1])), metadata_h1n1$Isolate_Id)]
age_breaks_h1n1 = c(0,20,40,60,80,120)
traits_h1n1 = cut(as.numeric(traits_h1n1), age_breaks_h1n1)
times_h1n1 = metadata_h1n1$Collection_Date[match(unlist(lapply(tree_h1n1$tip.label, function(x)str_split(x,'\\|')[[1]][1])), metadata_h1n1$Isolate_Id)]
demes_h1n1 = sort(unique(traits_h1n1))
hists_h1n1 = compute_histograms(tree = tree_h1n1,
                                traits = traits_h1n1,
                                times = times_h1n1,
                                demes = demes_h1n1,
                                
                                hist_dist_min_max = c(0,5),
                                hist_break_method = c('Manual'),
                                first_bin_npairs = seq(0,5,1),
                                
                                collection_time_diff_max = 3*4/365, ## 12 days (~ 4 generations)

                                bootstrap = 200)

traits_h3n2 = metadata_h3n2$Host_Age[match(unlist(lapply(tree_h3n2$tip.label, function(x)str_split(x,'\\|')[[1]][1])), metadata_h3n2$Isolate_Id)]
age_breaks_h3n2 = c(0,20,40,60,80,120)
traits_h3n2 = cut(as.numeric(traits_h3n2), age_breaks_h3n2)
times_h3n2 = metadata_h3n2$Collection_Date[match(unlist(lapply(tree_h3n2$tip.label, function(x)str_split(x,'\\|')[[1]][1])), metadata_h3n2$Isolate_Id)]
demes_h3n2 = sort(unique(traits_h3n2))
hists_h3n2 = compute_histograms(tree = tree_h3n2,
                                traits = traits_h3n2,
                                times = times_h3n2,
                                demes = demes_h3n2,
                                
                                hist_dist_min_max = c(0,5),
                                hist_break_method = c('Manual'),
                                first_bin_npairs = seq(0,5,1),
                                
                                collection_time_diff_max = 3*4/365, ## 12 days (~ 4 generations)

                                bootstrap = 200)


traits_RSVA = metadata_RSVA$Patient.age[match(unlist(lapply(tree_RSVA$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_RSVA$Accession.ID)]
age_breaks_RSVA = c(0,20,60,120)
traits_RSVA = cut(as.numeric(traits_RSVA), age_breaks_RSVA)
times_RSVA = metadata_RSVA$Collection.date[match(unlist(lapply(tree_RSVA$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_RSVA$Accession.ID)]
demes_RSVA = sort(unique(traits_RSVA))
hists_RSVA = compute_histograms(tree = tree_RSVA,
                                traits = traits_RSVA,
                                times = times_RSVA,
                                demes = demes_RSVA,
                                
                                hist_dist_min_max = c(0,5),
                                hist_break_method = c('Manual'),
                                first_bin_npairs = seq(0,5,1),
                                
                                collection_time_diff_max = 0.08, ## 1 month interval (~4 generations)

                                bootstrap = 200)

traits_RSVB = metadata_RSVB$Patient.age[match(unlist(lapply(tree_RSVB$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_RSVB$Accession.ID)]
age_breaks_RSVB = c(0,20,60,120)
traits_RSVB = cut(as.numeric(traits_RSVB), age_breaks_RSVB)
times_RSVB = metadata_RSVB$Collection.date[match(unlist(lapply(tree_RSVB$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_RSVB$Accession.ID)]
demes_RSVB = sort(unique(traits_RSVB))
hists_RSVB = compute_histograms(tree = tree_RSVB,
                                traits = traits_RSVB,
                                times = times_RSVB,
                                demes = demes_RSVB,
                                
                                hist_dist_min_max = c(0,5),
                                hist_break_method = c('Manual'),
                                first_bin_npairs = seq(0,5,1),
                                
                                collection_time_diff_max = 0.08, ## 1 month interval (~4 generations)

                                bootstrap = 200)
#########################################################################################################################################

#########################################################################################################################################
## Fit ratios
#########################################################################################################################################
bootstrap = 200
fit_SC2 = fit_hist_model(
  tree = tree_SC2,
  traits = traits_SC2,
  times = times_SC2,
  demes = demes_SC2,
  model = 'SYM', 
  parallel_chains = 3, chains = 3,
  n_bins_to_fit = rep(5, 10),
  iter_warmup = 2000, iter_sampling = 2000, refresh = 100,
  
  precomputed_hists = hists_SC2, stan_data_override = list(rate_dispersion_prior = 1),
  
  bootstrap = 0
)

fit_h3n2 = fit_hist_model(
  tree = tree_h3n2,
  traits = traits_h3n2,
  times = times_h3n2,
  demes = demes_h3n2,
  model = 'SYM',
  parallel_chains = 3, chains = 3,
  n_bins_to_fit = rep(5, 10),
  iter_warmup = 2000, iter_sampling = 2000,
  
  precomputed_hists = hists_h3n2, stan_data_override = list(rate_dispersion_prior = 1),
  
  bootstrap = 0
)
fit_h1n1 = fit_hist_model(
  tree = tree_h1n1,
  traits = traits_h1n1,
  times = times_h1n1,
  demes = demes_h1n1,
  model = 'SYM',
  parallel_chains = 3, chains = 3,
  n_bins_to_fit = rep(5, 10),
  iter_warmup = 2000, iter_sampling = 2000,
  
  precomputed_hists = hists_h1n1, stan_data_override = list(rate_dispersion_prior = 1),
  
  bootstrap = 0
)
fit_RSVA = fit_hist_model(
  tree = tree_RSVA,
  traits = traits_RSVA,
  times = times_RSVA,
  demes = demes_RSVA,
  model = 'SYM',
  parallel_chains = 3, chains = 3,
  n_bins_to_fit = rep(5, 3),
  iter_warmup = 2000, iter_sampling = 2000,
  
  precomputed_hists = hists_RSVA, stan_data_override = list(rate_dispersion_prior = 1),
  
  bootstrap = 0
)
fit_RSVB = fit_hist_model(
  tree = tree_RSVB,
  traits = traits_RSVB,
  times = times_RSVB,
  demes = demes_RSVB,
  model = 'SYM',
  parallel_chains = 3, chains = 3,
  n_bins_to_fit = rep(5, 3),
  iter_warmup = 2000, iter_sampling = 2000,
  
  precomputed_hists = hists_RSVB, stan_data_override = list(rate_dispersion_prior = 1),
  
  bootstrap = 0
)

saveRDS(fit_SC2, file = '5_age_groups/Data/SC2/fit_SC2_SYM.rds')
saveRDS(fit_h1n1, file = '5_age_groups/Data/H1N1/fit_h1n1_SYM.rds')
saveRDS(fit_h3n2, file = '5_age_groups/Data/H3N2/fit_h3n2_SYM.rds')
saveRDS(fit_RSVA, file = '5_age_groups/Data/RSVA/fit_RSVA_SYM.rds')
saveRDS(fit_RSVB, file = '5_age_groups/Data/RSVB/fit_RSVB_SYM.rds')

saveRDS(hists_SC2, file = '5_age_groups/Data/SC2/hists_SC2.rds')
saveRDS(hists_h1n1, file = '5_age_groups/Data/H1N1/hists_h1n1.rds')
saveRDS(hists_h3n2, file = '5_age_groups/Data/H3N2/hists_h3n2.rds')
saveRDS(hists_RSVA, file = '5_age_groups/Data/RSVA/hists_RSVA.rds')
saveRDS(hists_RSVB, file = '5_age_groups/Data/RSVB/hists_RSVB.rds')
#########################################################################################################################################



