#########################################################################################################################################
## Age group figures, for different pathogens in ARA
#########################################################################################################################################
library(ape)
library(stringr)
library(lubridate)

source('1_code/load_pairate.R', chdir = T)

#########################################################################################################################################
## Load raw data
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

#########################################################################################################################################
## Save Metadata
#########################################################################################################################################
write.csv(metadata_SC2, '5_age_groups/Data/SC2/Table_S2_List_isolates_SC2_age.csv')
write.csv(metadata_h1n1, '5_age_groups/Data/H1N1/Table_S3_List_isolates_h1n1_age.csv')
write.csv(metadata_h3n2, '5_age_groups/Data/H3N2/Table_S4_List_isolates_h3n2_age.csv')
write.csv(metadata_RSVA, '5_age_groups/Data/RSVA/Table_S5_List_isolates_RSVA_age.csv')
write.csv(metadata_RSVB, '5_age_groups/Data/RSVB/Table_S6_List_isolates_RSVB_age.csv')
#########################################################################################################################################

#########################################################################################################################################
## Age breaks
#########################################################################################################################################
traits_SC2 = metadata_SC2$Patient.age[match(unlist(lapply(tree_SC2$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_SC2$Accession.ID)]
age_breaks_SC2 = c(0,20,40,60,80,120) ## Technically 10-20 has 17 sequences
traits_SC2 = cut(as.numeric(traits_SC2), age_breaks_SC2)
times_SC2 = metadata_SC2$Collection.date[match(unlist(lapply(tree_SC2$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_SC2$Accession.ID)]
demes_SC2 = sort(unique(traits_SC2))
median.and.ci = function(v){ return( c(median(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))}

traits_h1n1 = metadata_h1n1$Host_Age[match(unlist(lapply(tree_h1n1$tip.label, function(x)str_split(x,'\\|')[[1]][1])), metadata_h1n1$Isolate_Id)]
age_breaks_h1n1 = c(0,20,40,60,80,120)
traits_h1n1 = cut(as.numeric(traits_h1n1), age_breaks_h1n1)
times_h1n1 = metadata_h1n1$Collection_Date[match(unlist(lapply(tree_h1n1$tip.label, function(x)str_split(x,'\\|')[[1]][1])), metadata_h1n1$Isolate_Id)]
demes_h1n1 = sort(unique(traits_h1n1))

traits_h3n2 = metadata_h3n2$Host_Age[match(unlist(lapply(tree_h3n2$tip.label, function(x)str_split(x,'\\|')[[1]][1])), metadata_h3n2$Isolate_Id)]
age_breaks_h3n2 = c(0,20,40,60,80,120)
traits_h3n2 = cut(as.numeric(traits_h3n2), age_breaks_h3n2)
times_h3n2 = metadata_h3n2$Collection_Date[match(unlist(lapply(tree_h3n2$tip.label, function(x)str_split(x,'\\|')[[1]][1])), metadata_h3n2$Isolate_Id)]
demes_h3n2 = sort(unique(traits_h3n2))

traits_RSVA = metadata_RSVA$Patient.age[match(unlist(lapply(tree_RSVA$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_RSVA$Accession.ID)]
age_breaks_RSVA = c(0,20,60,120)
traits_RSVA = cut(as.numeric(traits_RSVA), age_breaks_RSVA)
times_RSVA = metadata_RSVA$Collection.date[match(unlist(lapply(tree_RSVA$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_RSVA$Accession.ID)]
demes_RSVA = sort(unique(traits_RSVA))

traits_RSVB = metadata_RSVB$Patient.age[match(unlist(lapply(tree_RSVB$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_RSVB$Accession.ID)]
age_breaks_RSVB = c(0,20,60,120)
traits_RSVB = cut(as.numeric(traits_RSVB), age_breaks_RSVB)
times_RSVB = metadata_RSVB$Collection.date[match(unlist(lapply(tree_RSVB$tip.label, function(x)str_split(x,'\\|')[[1]][2])), metadata_RSVB$Accession.ID)]
demes_RSVB = sort(unique(traits_RSVB))
#########################################################################################################################################

#########################################################################################################################################
## Load fit and hist data
#########################################################################################################################################
fit_SC2 = readRDS(file = '5_age_groups/Data/SC2/fit_SC2_SYM.rds')
fit_h1n1 = readRDS(file = '5_age_groups/Data/H1N1/fit_h1n1_SYM.rds')
fit_h3n2 = readRDS(file = '5_age_groups/Data/H3N2/fit_h3n2_SYM.rds')
fit_RSVA = readRDS(file = '5_age_groups/Data/RSVA/fit_RSVA_SYM.rds')
fit_RSVB = readRDS(file = '5_age_groups/Data/RSVB/fit_RSVB_SYM.rds')

hists_SC2 = readRDS(file = '5_age_groups/Data/SC2/hists_SC2.rds')
hists_h1n1 = readRDS(file = '5_age_groups/Data/H1N1/hists_h1n1.rds')
hists_h3n2 = readRDS(file = '5_age_groups/Data/H3N2/hists_h3n2.rds')
hists_RSVA = readRDS(file = '5_age_groups/Data/RSVA/hists_RSVA.rds')
hists_RSVB = readRDS(file = '5_age_groups/Data/RSVB/hists_RSVB.rds')
#########################################################################################################################################

#########################################################################################################################################
## Precompute - estimates in matrices
#########################################################################################################################################
hists = hists_h1n1
rates = fit_h1n1$chains$rates
ts = 365/(substitution_rate_h1n1*genome_length_h1n1)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
estimated_mat = matrix(0, ncol = length(hists$demes_in_tree), nrow = length(hists$demes_in_tree))
estimated_mat[lower.tri(estimated_mat)] = res_ci[1,]/ts
estimated_mat = estimated_mat+t(estimated_mat)
diag(estimated_mat) = NA
colnames(estimated_mat) = hists$demes_in_tree
rownames(estimated_mat) = hists$demes_in_tree
estimated_mat_h1n1 = estimated_mat

hists = hists_h3n2
rates = fit_h3n2$chains$rates
ts = 365/(substitution_rate_h3n2*genome_length_h3n2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
estimated_mat = matrix(0, ncol = length(hists$demes_in_tree), nrow = length(hists$demes_in_tree))
estimated_mat[lower.tri(estimated_mat)] = res_ci[1,]/ts
estimated_mat = estimated_mat+t(estimated_mat)
diag(estimated_mat) = NA
colnames(estimated_mat) = hists$demes_in_tree
rownames(estimated_mat) = hists$demes_in_tree
estimated_mat_h3n2 = estimated_mat

hists = hists_RSVA
rates = fit_RSVA$chains$rates
ts = 365/(substitution_rate_RSVA*genome_length_RSVB)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
estimated_mat = matrix(0, ncol = length(hists$demes_in_tree), nrow = length(hists$demes_in_tree))
estimated_mat[lower.tri(estimated_mat)] = res_ci[1,]/ts
estimated_mat = estimated_mat+t(estimated_mat)
diag(estimated_mat) = NA
colnames(estimated_mat) = hists$demes_in_tree
rownames(estimated_mat) = hists$demes_in_tree
estimated_mat_RSVA = estimated_mat

hists = hists_RSVB
rates = fit_RSVB$chains$rates
ts = 365/(substitution_rate_RSVB*genome_length_RSVB)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
estimated_mat = matrix(0, ncol = length(hists$demes_in_tree), nrow = length(hists$demes_in_tree))
estimated_mat[lower.tri(estimated_mat)] = res_ci[1,]/ts
estimated_mat = estimated_mat+t(estimated_mat)
diag(estimated_mat) = NA
colnames(estimated_mat) = hists$demes_in_tree
rownames(estimated_mat) = hists$demes_in_tree
estimated_mat_RSVB = estimated_mat

hists = hists_SC2
rates = fit_SC2$chains$rates
ts = 365/(substitution_rate_SC2*genome_length_SC2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
estimated_mat = matrix(0, ncol = length(hists$demes_in_tree), nrow = length(hists$demes_in_tree))
diag(estimated_mat) = NA
estimated_mat[lower.tri(estimated_mat)] = res_ci[1,]/ts
estimated_mat = estimated_mat+t(estimated_mat)
colnames(estimated_mat) = hists$demes_in_tree
rownames(estimated_mat) = hists$demes_in_tree
estimated_mat_SC2 = estimated_mat

## Plot
estimated_mats = list(
  H1N1 = estimated_mat_h1n1,
  H3N2 = estimated_mat_h3n2,                                                                                                                                                        
  RSVA = estimated_mat_RSVA,
  RSVB = estimated_mat_RSVB,                                                                                                                                                        
  SC2  = estimated_mat_SC2                                                                                                                                                        
)
#########################################################################################################################################

#########################################################################################################################################
## Precompute - matrices N seqs
#########################################################################################################################################
## SC2
N_seqs_SC2 = matrix(0, ncol = length(hists_SC2$demes_in_tree), nrow = length(hists_SC2$demes_in_tree))
colnames(N_seqs_SC2) = hists_SC2$demes_in_tree
rownames(N_seqs_SC2) = hists_SC2$demes_in_tree
for(i in 1:length(hists_SC2$demes_in_tree)){
  for(j in 1:length(hists_SC2$demes_in_tree)){
    N = length(which(traits_SC2 == hists_SC2$demes_in_tree[i])) + length(which(traits_SC2 == hists_SC2$demes_in_tree[j]))
    N_seqs_SC2[i,j] = N
  }
}

## H1N1
N_seqs_h1n1 = matrix(0, ncol = length(hists_h1n1$demes_in_tree), nrow = length(hists_h1n1$demes_in_tree))
colnames(N_seqs_h1n1) = hists_h1n1$demes_in_tree
rownames(N_seqs_h1n1) = hists_h1n1$demes_in_tree
for(i in 1:length(hists_h1n1$demes_in_tree)){
  for(j in 1:length(hists_h1n1$demes_in_tree)){
    N = length(which(traits_h1n1 == hists_h1n1$demes_in_tree[i])) + length(which(traits_h1n1 == hists_h1n1$demes_in_tree[j]))
    N_seqs_h1n1[i,j] = N
  }
}

## H3N2
N_seqs_h3n2 = matrix(0, ncol = length(hists_h3n2$demes_in_tree), nrow = length(hists_h3n2$demes_in_tree))
colnames(N_seqs_h3n2) = hists_h3n2$demes_in_tree
rownames(N_seqs_h3n2) = hists_h3n2$demes_in_tree
for(i in 1:length(hists_h3n2$demes_in_tree)){
  for(j in 1:length(hists_h3n2$demes_in_tree)){
    N = length(which(traits_h3n2 == hists_h3n2$demes_in_tree[i])) + length(which(traits_h3n2 == hists_h3n2$demes_in_tree[j]))
    N_seqs_h3n2[i,j] = N
  }
}

## RSVA
N_seqs_RSVA = matrix(0, ncol = length(hists_RSVA$demes_in_tree), nrow = length(hists_RSVA$demes_in_tree))
colnames(N_seqs_RSVA) = hists_RSVA$demes_in_tree
rownames(N_seqs_RSVA) = hists_RSVA$demes_in_tree
for(i in 1:length(hists_RSVA$demes_in_tree)){
  for(j in 1:length(hists_RSVA$demes_in_tree)){
    N = length(which(traits_RSVA == hists_RSVA$demes_in_tree[i])) + length(which(traits_RSVA == hists_RSVA$demes_in_tree[j]))
    N_seqs_RSVA[i,j] = N
  }
}

## RSVB
N_seqs_RSVB = matrix(0, ncol = length(hists_RSVB$demes_in_tree), nrow = length(hists_RSVB$demes_in_tree))
colnames(N_seqs_RSVB) = hists_RSVB$demes_in_tree
rownames(N_seqs_RSVB) = hists_RSVB$demes_in_tree
for(i in 1:length(hists_RSVB$demes_in_tree)){
  for(j in 1:length(hists_RSVB$demes_in_tree)){
    N = length(which(traits_RSVB == hists_RSVB$demes_in_tree[i])) + length(which(traits_RSVB == hists_RSVB$demes_in_tree[j]))
    N_seqs_RSVB[i,j] = N
  }
}
#########################################################################################################################################

#########################################################################################################################################
## Precompute - contact matrices
#########################################################################################################################################
library(contactsurveys)
library(socialmixr)
library(data.table)

# Beraud survey
fr_files <- contactsurveys::download_survey('https://doi.org/10.5281/zenodo.3886590')
fr_survey <- socialmixr::load_survey(fr_files)

## Load population data from INSEE
pop_insee = read.csv("5_age_groups/Data/contact_survey/Pop_Insee/DS_RP_TD_POPULATION_AGESEX_PRINC_2023_CSV_FR/DS_RP_TD_POPULATION_AGESEX_PRINC_2023_data.csv", sep = ';', header = T)
pop_insee = pop_insee = pop_insee[which(pop_insee$SEX == '_T' & pop_insee$GEO_OBJECT == 'FRANCE' & pop_insee$GEO == 'FM' & pop_insee$AGE != '_T'),]
pop_insee$AGE = sapply(pop_insee$AGE, function(x)str_split(x, 'Y')[[1]][2])
pop_insee$AGE[which(pop_insee$AGE == '_GE100')] = 100
pop_insee$AGE = as.numeric(pop_insee$AGE)

# Functions to compute the contact matrix
get_cm = function(age_limits){
  cm_fr <- socialmixr::contact_matrix(
    survey      = fr_survey,
    age_limits  = age_limits,
    symmetric   = FALSE,   # set FALSE if you want the raw (non-reciprocal) matrix
    counts      = FALSE   # FALSE = mean contacts per person per day
  )
  ## Load population data from INSEE
  breaks  = c(age_limits, Inf)
  pop_insee$age_group = cut(pop_insee$AGE, 
                          breaks = breaks, right = FALSE, include.lowest = TRUE)
  N = tapply(pop_insee$OBS_VALUE, pop_insee$age_group, sum)
  N = as.numeric(N)
  
  # Population-weighted mean: (N_i * Q_contact[i,j] + N_j * Q_contact[j,i]) / (N_i + N_j)
  N_i = outer(N, N, function(ni, nj) ni)  # matrix of N_i values
  N_j = t(N_i)                             # matrix of N_j values
  
  Q_contact = cm_fr$matrix
  Q_contact = (N_i * Q_contact + N_j * t(Q_contact)) / (N_i + N_j)
  
  diag(Q_contact) = NA
  return(Q_contact)
}
get_cm_socialcov = function(age_limits){
  # SocialCov contact matrix from May 2022
  cm = read.csv("5_age_groups/Data/contact_survey/SocialCov/Mat6_regular.csv", row.names = 1)
  Q_contact = as.matrix(cm)
  
  # Population for each 10-year fine bin (to match SocialCov breaks)
  fine_breaks_SocialCov = c(0, 10, 20, 30, 40, 50, 60, 70, Inf)
  pop_insee$age_group_SocialCov = cut(pop_insee$AGE, breaks = fine_breaks_SocialCov, right = FALSE)
  N_fine_SocialCov = as.numeric(tapply(pop_insee$OBS_VALUE, pop_insee$age_group_SocialCov, sum))  # length 8
  
  # Split proportions for 70+ -> [70,80) and [80+)
  N_70_80 = sum(pop_insee$OBS_VALUE[pop_insee$AGE >= 70 & pop_insee$AGE < 80])
  N_80p   = sum(pop_insee$OBS_VALUE[pop_insee$AGE >= 80])
  N_70p   = N_70_80 + N_80p
  p_70_80 = N_70_80 / N_70p
  p_80p   = N_80p   / N_70p
  
  # Expand SocialCov 8x8 matrix -> 9x9 by splitting last group
  Q_contact_exp = matrix(0, 9, 9)
  Q_contact_exp[1:7, 1:7] = Q_contact[1:7, 1:7]
  # Column split (contacts WITH 70+): proportional to population
  Q_contact_exp[1:7, 8]   = Q_contact[1:7, 8] * p_70_80
  Q_contact_exp[1:7, 9]   = Q_contact[1:7, 8] * p_80p
  # Row split (contacts FROM 70+): same per-person rate assumed
  Q_contact_exp[8,   1:7] = Q_contact[8, 1:7]
  Q_contact_exp[9,   1:7] = Q_contact[8, 1:7]
  # Bottom-right 2x2 block
  Q_contact_exp[8, 8] = Q_contact[8, 8] * p_70_80
  Q_contact_exp[8, 9] = Q_contact[8, 8] * p_80p
  Q_contact_exp[9, 8] = Q_contact[8, 8] * p_70_80
  Q_contact_exp[9, 9] = Q_contact[8, 8] * p_80p
  
  # Population for 9 expanded groups
  N_exp = c(N_fine_SocialCov[1:7], N_70_80, N_80p)
  
  fine_limits  = c(0,10,20,30,40,50,60,70,80)
  coarse_breaks = c(age_limits, Inf)
  coarse_idx = findInterval(fine_limits, coarse_breaks)
  groups = lapply(seq_along(age_limits), function(j) which(coarse_idx == j))
  
  n_coarse = length(groups)
  Q_contact_coarse = matrix(0, n_coarse, n_coarse, dimnames = list(names(groups), names(groups)))
  N_coarse = numeric(n_coarse)
  
  for (I in seq_along(groups)) {
    N_coarse[I] = sum(N_exp[groups[[I]]])
    for (J in seq_along(groups)) {
      fine_i   = groups[[I]]
      fine_j   = groups[[J]]
      col_sums = rowSums(Q_contact_exp[fine_i, fine_j, drop = FALSE])
      Q_contact_coarse[I, J] = sum(N_exp[fine_i] * col_sums) / sum(N_exp[fine_i])
    }
  }
  
  # Population-weighted symmetrization
  N_i   = outer(N_coarse, N_coarse, function(ni, nj) ni)
  N_j   = t(N_i)
  Q_contact_sym = (N_i * Q_contact_coarse + N_j * t(Q_contact_coarse)) / (N_i + N_j)
  
  diag(Q_contact_sym) = NA
  return(Q_contact_sym)
}
get_cm_socialcov_asym = function(age_limits){
  # SocialCov contact matrix from May 2022
  cm = read.csv("5_age_groups/Data/contact_survey/SocialCov/Mat6_regular.csv", row.names = 1)
  Q_contact = as.matrix(cm)
  
  # Population for each 10-year fine bin (to match SocialCov breaks)
  fine_breaks_SocialCov = c(0, 10, 20, 30, 40, 50, 60, 70, Inf)
  pop_insee$age_group_SocialCov = cut(pop_insee$AGE, breaks = fine_breaks_SocialCov, right = FALSE)
  N_fine_SocialCov = as.numeric(tapply(pop_insee$OBS_VALUE, pop_insee$age_group_SocialCov, sum))  # length 8
  
  # Split proportions for 70+ -> [70,80) and [80+)
  N_70_80 = sum(pop_insee$OBS_VALUE[pop_insee$AGE >= 70 & pop_insee$AGE < 80])
  N_80p   = sum(pop_insee$OBS_VALUE[pop_insee$AGE >= 80])
  N_70p   = N_70_80 + N_80p
  p_70_80 = N_70_80 / N_70p
  p_80p   = N_80p   / N_70p
  
  # Expand SocialCov 8x8 matrix -> 9x9 by splitting last group
  Q_contact_exp = matrix(0, 9, 9)
  Q_contact_exp[1:7, 1:7] = Q_contact[1:7, 1:7]
  # Column split (contacts WITH 70+): proportional to population
  Q_contact_exp[1:7, 8]   = Q_contact[1:7, 8] * p_70_80
  Q_contact_exp[1:7, 9]   = Q_contact[1:7, 8] * p_80p
  # Row split (contacts FROM 70+): same per-person rate assumed
  Q_contact_exp[8,   1:7] = Q_contact[8, 1:7]
  Q_contact_exp[9,   1:7] = Q_contact[8, 1:7]
  # Bottom-right 2x2 block
  Q_contact_exp[8, 8] = Q_contact[8, 8] * p_70_80
  Q_contact_exp[8, 9] = Q_contact[8, 8] * p_80p
  Q_contact_exp[9, 8] = Q_contact[8, 8] * p_70_80
  Q_contact_exp[9, 9] = Q_contact[8, 8] * p_80p
  
  # Population for 9 expanded groups
  N_exp = c(N_fine_SocialCov[1:7], N_70_80, N_80p)
  
  # Aggregate 9 groups -> 5 coarse groups
  groups = list("[0,20)"  = 1:2,
                "[20,40)" = 3:4,
                "[40,60)" = 5:6,
                "[60,80)" = 7:8,
                "[80+)"   =   9)
  
  n_coarse = length(groups)
  Q_contact_coarse = matrix(0, n_coarse, n_coarse, dimnames = list(names(groups), names(groups)))
  N_coarse = numeric(n_coarse)
  
  for (I in seq_along(groups)) {
    N_coarse[I] = sum(N_exp[groups[[I]]])
    for (J in seq_along(groups)) {
      fine_i   = groups[[I]]
      fine_j   = groups[[J]]
      col_sums = rowSums(Q_contact_exp[fine_i, fine_j, drop = FALSE])
      Q_contact_coarse[I, J] = sum(N_exp[fine_i] * col_sums) / sum(N_exp[fine_i])
    }
  }
  
  # Population-weighted symmetrization
  # N_i   = outer(N_coarse, N_coarse, function(ni, nj) ni)
  # N_j   = t(N_i)
  # Q_contact_sym = (N_i * Q_contact_coarse + N_j * t(Q_contact_coarse)) / (N_i + N_j)
  
  diag(Q_contact_coarse) = NA
  return(Q_contact_coarse)
}
get_cm_asym = function(age_limits){
  cm_fr <- socialmixr::contact_matrix(
    survey      = fr_survey,
    age_limits  = age_limits,
    symmetric   = F,   # set FALSE if you want the raw (non-reciprocal) matrix
    counts      = FALSE   # FALSE = mean contacts per person per day
  )
  Q_contact = cm_fr$matrix
  diag(Q_contact) = NA
  diag(Q_contact) = -rowSums(Q_contact)
  return(Q_contact)
}

Q_contact_sym_SC2 = get_cm(age_breaks_SC2[-length(age_breaks_SC2)])
Q_contact_sym_SC2_socialcov = get_cm_socialcov(age_breaks_SC2[-length(age_breaks_SC2)])
Q_contact_sym_h3n2 = get_cm(age_breaks_h3n2[-length(age_breaks_h3n2)])
Q_contact_sym_h3n2_socialcov = get_cm_socialcov(age_breaks_h3n2[-length(age_breaks_h3n2)])
Q_contact_sym_h1n1 = get_cm(age_breaks_h1n1[-length(age_breaks_h1n1)])
Q_contact_sym_h1n1_socialcov = get_cm_socialcov(age_breaks_h1n1[-length(age_breaks_h1n1)])
Q_contact_sym_RSVA = get_cm(age_breaks_RSVA[-length(age_breaks_RSVA)])
Q_contact_sym_RSVA_socialcov = get_cm_socialcov(age_breaks_RSVA[-length(age_breaks_RSVA)])
Q_contact_sym_RSVB = get_cm(age_breaks_RSVB[-length(age_breaks_RSVB)])
Q_contact_sym_RSVB_socialcov = get_cm_socialcov(age_breaks_RSVB[-length(age_breaks_RSVB)])

Q_contact_asym_SC2 = get_cm_asym(age_breaks_SC2[-length(age_breaks_SC2)])
Q_contact_asym_SC2_socialcov = get_cm_socialcov_asym(age_breaks_SC2[-length(age_breaks_SC2)])
#########################################################################################################################################


#########################################################################################################################################
## Figure 4: SC2 between age groups
#########################################################################################################################################
# pdf(file = '5_age_groups/Figures/Figure4/Figure4AB_Data_counts.pdf', width = 8, height = 2.5)
par(mfrow = c(1,2), mar = c(2,2,1,2), oma = c(1,1,1,1))
# hist(metadata_SC2$Collection.date[which(!is.na(metadata_SC2$Patient.age))], breaks = seq(2024.45,2025.6,1/52), main = 'SARS-CoV-2', xlab = 'Age', ylab = 'Count', yaxt = 'n', xaxt = 'n', ylim = c(0,200))
hist(metadata_SC2$Collection.date, breaks = seq(2024.45,2025.6,1/52), main = 'SARS-CoV-2', xlab = 'Age', ylab = 'Count', yaxt = 'n', xaxt = 'n', ylim = c(0,200))
axis(1, at = decimal_date(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01'))),
     paste0(Months[month(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01')))], ' ', c(24,24,24,25,25,25,25)))
axis(2, las = 2, at = seq(0,200,50), labels = seq(0,200,50))
hist(metadata_SC2$Patient.age, breaks = seq(0,125,5), main = 'SARS-CoV-2', xlab = 'Age', ylab = 'Count', yaxt = 'n', xlim = c(0,100),
     col = c(rep('grey30', 4), rep('grey', 4), rep('grey30', 4), rep('grey', 4), rep('grey30', 9)))
axis(2, las = 2, at = seq(0,120,40))
# dev.off()


# pdf(file = '5_age_groups/Figures/Figure4/Figure4CD_Ratios_13102026.pdf', width = 8, height = 2.5)
par(mfrow = c(1,3), mar = c(2,2,1,2), oma = c(1,1,1,1))
hists = hists_SC2
res_within = res_between = matrix(0, nrow = dim(hists$pair_hists_boot$within_a)[2], ncol = 200)
for(i in 1:dim(hists$pair_hists_boot$within_a)[1]){
  for(b in 1:200){
    res_within[,b] = res_within[,b] + (hists$pair_hists_boot$within_a[i,,b]/sum(hists$pair_hists_boot$within_a[i,,b])) + (hists$pair_hists_boot$within_b[i,,b]/sum(hists$pair_hists_boot$within_b[i,,b]))
    res_between[,b] = res_between[,b] + (hists$pair_hists_boot$between[i,,b]/sum(hists$pair_hists_boot$between[i,,b]))
  }
}
res_within = apply(res_within, MARGIN = 2, function(x)x/sum(x))
res_between = apply(res_between, MARGIN = 2, function(x)x/sum(x))
tmp = sapply(1:ncol(res_within), function(x)res_within[,x]/res_between[,x])
tmp[which(is.infinite(tmp))] = NA
tmp = apply(tmp, MAR = 1, median_and_ci)
tmp[2,which(tmp[2,] == 0)] = 0.001
tmp[3,which(is.infinite(tmp[3,]))] = 1000
plot((hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2, tmp[1,], log = 'y', pch = 16, bty = 'n', yaxt = 'n',  xaxt = 'n',
     xlim = c(0,5), ylim = c(0.5, 100), ylab = 'Ratio', xlab = 'Genetic distance', main = 'SARS-CoV-2')
axis(1, at = seq(0.5, 4.5, 1), labels = paste0('[', 0:4, '-', 1:5, '['))
axis(2, at = c(0.5, 1, 10, 100), labels = c('<0.5', 1, 10, 100), las = 2, tick = c(0.1, 0.5, 1, 5, 10, 50, 100))
arrows(x0 = (hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2, y0 = tmp[2,],
       x1 = (hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2, y1 = tmp[3,],
       length = 0.05, angle = 90, code = 0, lwd = 1.5)
abline(h = 1, lty = 2)

hists = hists_SC2
plot(NULL, log = 'y', pch = 16, bty = 'n', yaxt = 'n', xaxt = 'n',
     xlim = c(0,5), ylim = c(0.5, 100), ylab = 'Ratio', xlab = 'Genetic distance', main = 'SARS-CoV-2')
axis(1, at = seq(0.5, 4.5, 1), labels = paste0('[', 0:4, '-', 1:5, '['))
axis(2, at = c(0.5, 1, 10, 100), labels = c('<0.5', 1, 10, 100), las = 2, tick = c(0.1, 0.5, 1, 5, 10, 50, 100))
abline(h = 1, lty = 2)
for(p in 1:length(hists$pair_labels)){
  res_within = res_between = matrix(0, nrow = dim(hists$pair_hists_boot$within_a)[2], ncol = 200)
  for(b in 1:200){
    res_within[,b] = res_within[,b] + (hists$pair_hists_boot$within_a[p,,b]/sum(hists$pair_hists_boot$within_a[p,,b])) + (hists$pair_hists_boot$within_b[p,,b]/sum(hists$pair_hists_boot$within_b[p,,b]))
    res_between[,b] = res_between[,b] + (hists$pair_hists_boot$between[p,,b]/sum(hists$pair_hists_boot$between[p,,b]))
  }
  res_within = apply(res_within, MARGIN = 2, function(x)x/sum(x))
  res_between = apply(res_between, MARGIN = 2, function(x)x/sum(x))
  tmp = sapply(1:ncol(res_within), function(x)res_within[,x]/res_between[,x])
  tmp[which(is.infinite(tmp))] = NA
  tmp = apply(tmp, MAR = 1, median_and_ci)
  tmp[2,which(tmp[2,] < 0.5 )] = 0.5
  tmp[3,which(is.infinite(tmp[3,]))] = 1000
  points((hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2 + ((p-length(hists$pair_labels)/2))/(length(hists$pair_labels)/2)*0.25 , 
         tmp[1,], pch = 16, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
  arrows(x0 = (hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y0 = tmp[2,],
         x1 = (hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y1 = tmp[3,],
         length = 0.05, angle = 90, code = 0, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]),
         lwd = 1.5)
}
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', xlim = c(0,1), ylim = c(1,1))
legend('center', col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))), pch = 15, legend = hists$pair_labels, bty = 'n')
# dev.off()

# pdf(file = '5_age_groups/Figures/Figure4/Figure4E_Fits_13102026.pdf', width = 8, height = 2.5)
par(mfrow = c(1,3), mar = c(2,2,1,2), oma = c(1,1,1,1))
hists = hists_SC2
fit = fit_SC2
fit_summary = plot_proportion_fits(fit, plot = F)
p_names = which(lower.tri(Q_contact_sym_SC2)==T, arr.ind = T)
p_names = paste0(p_names[,2], p_names[,1])
pair_colors = c()
for(p in 1:fit$stan_data$N_sets){
  pr = fit$stan_data$set_pair[p]
  p_id = paste0(sort(c(fit$stan_data$pair_a[pr], fit$stan_data$pair_b[pr])), collapse = '')
  p_id = match(p_id, p_names)
  pair_colors[p] = RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p_id]
}
plot(NULL, xlim = c(10^-(2),1), ylim = c(10^-(2),1), bty = 'n', xaxt = 'n',  yaxt = 'n', log = 'xy', main = 'P same, SARS-CoV-2', xlab = 'Observed', ylab = 'Predicted (median)')
axis(1, at = c(10^(seq(-2,0, 1))), labels = c(10^(seq(-2,0, 1))))
axis(2, las = 2, at = c(10^(seq(-2,0, 1))), labels = c(10^(seq(-2,0, 1))))
abline(a=0, b=1, lty= 2)
for(p in 1:fit$stan_data$N_sets){
  bins = which(fit_summary$fitted[p,])
  arrows(x0 = fit_summary$observed_same[p,bins], y0 = fit_summary$predicted_same_lo[p,bins], x1 = fit_summary$observed_same[p,bins], y1 = fit_summary$predicted_same_hi[p,bins], length = 0, col = pair_colors[p])
  points(fit_summary$observed_same[p,bins], fit_summary$predicted_same[p,bins], pch = 16, col = pair_colors[p])
}
plot(NULL, xlim = c(10^-(2),1), ylim = c(10^-(2),1), bty = 'n', xaxt = 'n',  yaxt = 'n', log = 'xy', main = 'P diff, SARS-CoV-2', xlab = 'Observed', ylab = 'Predicted (median)')
axis(1, at = c(10^(seq(-2,0, 1))), labels = c(10^(seq(-2,0, 1))))
axis(2, las = 2, at = c(10^(seq(-2,0, 1))), labels = c(10^(seq(-2,0, 1))))
abline(a=0, b=1, lty= 2)
for(p in 1:fit$stan_data$N_sets){
  bins = which(fit_summary$fitted[p,])
  arrows(x0 = fit_summary$observed_diff[p,bins], y0 = fit_summary$predicted_diff_lo[p,bins], x1 = fit_summary$observed_diff[p,bins], y1 = fit_summary$predicted_diff_hi[p,bins], length = 0, col = pair_colors[p])
  points(fit_summary$observed_diff[p,bins], fit_summary$predicted_diff[p,bins], pch = 16, col = pair_colors[p])
}
plot(NULL, xlim = c(0.5,10), ylim = c(0.5,10), bty = 'n', xaxt = 'n',  yaxt = 'n', log = 'xy', main = 'Ratio, SARS-CoV-2', xlab = 'Observed', ylab = 'Predicted (median)')
axis(1, at = c(0.5, 1, 5, 10), labels = c('<0.5', 1, 5, 10), las = 1, tick = c(0.5, 1, 5, 10))
axis(2, las = 2, at = c(0.5, 1, 5, 10), labels = c('<0.5', 1, 5, 10), las = 1, tick = c(0.5, 1, 5, 10))
abline(a=0, b=1, lty= 2)
for(p in 1:fit$stan_data$N_sets){
  pr = fit$stan_data$set_pair[p]
  w = fit$stan_data$ref_weight[p,]
  bins = which(fit_summary$fitted[p,])
  ref = 0
  for(i in 1:fit$stan_data$N_boot){
    same = (w[1]*fit$stan_data$hist_within_a[pr,bins,i]/fit$stan_data$n_within_a[pr,i] + w[2]*fit$stan_data$hist_within_b[pr,bins,i]/fit$stan_data$n_within_b[pr,i])/(w[1]+w[2])
    ref = ref + ((1-w[3])*same + w[3]*fit$stan_data$hist_diff[pr,bins,i]/fit$stan_data$n_diff[pr,i])/fit$stan_data$N_boot
  }
  pred_same = t(apply(1 - w[3]*fit$chains$cr[,p,bins], MAR = 1, function(x)(x*ref)/sum(x*ref)))
  pred_diff = t(apply(fit$chains$cr[,p,bins], MAR = 1, function(x)(x*ref)/sum(x*ref)))
  ci = apply(pred_same/pred_diff, MARGIN = 2, FUN = median_and_ci)
  data = fit_summary$observed_same[p,bins]/fit_summary$observed_diff[p,bins]
  arrows(x0 = data, y0 = ci[2,], x1 = data, y1 = ci[3,], length = 0, col = pair_colors[p])
  points(data, ci[1,], pch = 16, col = pair_colors[p])
}
# dev.off()

# pdf(file = '5_age_groups/Figures/Figure4/Figure4F_Matrices_13102026.pdf', width = 8, height = 2.5)
par(mfrow = c(1,3), mar = c(2,2,1,2), oma = c(1,1,1,1))
cols = hcl.colors(50, "YlOrRd", rev = TRUE)
zlim = c(0, range(sapply(estimated_mats, range, na.rm = TRUE))[2])
mat = t((Q_contact_sym_SC2))
image(1:nrow(mat), 1:ncol(mat), mat,
      axes = FALSE, xlab = "", ylab = "",
      col = cols,
      main = 'Contact matrix')
axis(1, at = 1:nrow(mat), labels = rownames(mat), las = 2, cex.axis = 0.7)
axis(2, at = 1:ncol(mat), labels = colnames(mat), las = 2, cex.axis = 0.7)

mat = (estimated_mat_SC2)
image(1:nrow(mat), 1:ncol(mat), mat,
      axes = FALSE, xlab = "", ylab = "",
      col = cols, #zlim = zlim,
      main = names(estimated_mat_SC2)[i])
axis(1, at = 1:nrow(mat), labels = rownames(mat), las = 2, cex.axis = 0.7)
axis(2, at = 1:ncol(mat), labels = colnames(mat), las = 2, cex.axis = 0.7)

plot.new()
par(usr = c(0, 1, zlim[1], zlim[2]))
rasterImage(as.raster(matrix(rev(cols), ncol = 1)), 0.2, zlim[1], 0.4, zlim[2])
axis(4, at = pretty(zlim), las = 1, cex.axis = 0.8)
# dev.off()

# pdf(file = '5_age_groups/Figures/Figure4/Figure4G_Rates_13102026.pdf', width = 8, height = 2.5)
par(mfrow = c(1,3), mar = c(2,2,1,2), oma = c(1,1,1,1))
plot.new()
plot.new()

hists = hists_SC2
Q_contact_sym = Q_contact_sym_SC2_socialcov
Q_contact_sym[upper.tri(Q_contact_sym)] = NA
rates = fit_SC2$chains$rates
ts = 365/(substitution_rate_SC2*genome_length_SC2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
idx = which(!is.na(Q_contact_sym), arr.ind = T)
plot(Q_contact_sym[idx], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(0.001), max(1)), xlim = c(0.2, 5), xaxt = 'n',
     bty = 'n', yaxt = 'n', log = 'xy', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'Social contact rates in the population (per day)', ylab = 'Rates from genetic distances (per day)', main = 'SARS-CoV-2')
axis(1, at = c(0.2, 0.5, 1, 2, 5))
axis(2, las = 2, at = c(0.001, 0.01, 0.1, 1))
arrows(x0 = Q_contact_sym[idx], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[idx], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
o = order(Q_contact_sym[lower.tri(Q_contact_sym)])
x = Q_contact_sym[lower.tri(Q_contact_sym)][o]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)[o]
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-2, 2, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the age groups
mat = Q_contact_sym_SC2_socialcov
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value
# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Figure 5: age group spread across different pathogens
#########################################################################################################################################
# pdf(file = '5_age_groups/Figures/Figure5/Figure5_13102026.pdf', width = 8, height = 7)
par(mfcol = c(4,4), mar = c(2,2,2,2), oma = c(1,1,1,1))
hist(metadata_h1n1$Collection_Date[which(!is.na(metadata_h1n1$Host_Age))], breaks = seq(2024.45,2025.5,1/52), 
     main = 'H1N1', xlab = 'Age', ylab = 'Count', yaxt = 'n', xaxt = 'n', ylim = c(0,150))
axis(2, las = 2,
     at = seq(0,150,50), 
     labels = seq(0,150,50))
axis(1, at = decimal_date(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01'))),
     paste0(Months[month(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01')))], ' ', c(24,24,24,25,25,25,25)))
hist(metadata_h3n2$Collection_Date[which(!is.na(metadata_h3n2$Host_Age))], breaks = seq(2024.45,2025.5,1/52), 
     main = 'H3N2', xlab = 'Age', ylab = 'Count', yaxt = 'n', xaxt = 'n', ylim = c(0,150))
axis(2, las = 2, at = seq(0,150,50), labels = seq(0,150,50))
axis(1, at = decimal_date(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01'))),
     paste0(Months[month(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01')))], ' ', c(24,24,24,25,25,25,25)))
hist(metadata_RSVA$Collection.date[which(!is.na(metadata_RSVA$Patient.age))], breaks = seq(2024.45,2025.6,1/52), main = 'RSV A', xlab = 'Age', ylab = 'Count', yaxt = 'n', xaxt = 'n', ylim = c(0,20))
axis(2, las = 2, at = seq(0,20,5), labels = seq(0,20,5))
axis(1, at = decimal_date(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01'))),
     paste0(Months[month(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01')))], ' ', c(24,24,24,25,25,25,25)))
hist(metadata_RSVB$Collection.date[which(!is.na(metadata_RSVB$Patient.age))], breaks = seq(2024.45,2025.6,1/52), main = 'RSV B', xlab = 'Age', ylab = 'Count', yaxt = 'n', xaxt = 'n', ylim = c(0,50))
axis(2, las = 2, at = seq(0,40,10), labels = seq(0,40,10))
axis(1, at = decimal_date(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01'))),
     paste0(Months[month(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01')))], ' ', c(24,24,24,25,25,25,25)))

## Histogram by age
hist(metadata_h1n1$Host_Age, breaks = seq(0,125,5), main = 'H1N1', xlab = 'Age', ylab = 'Count', yaxt = 'n', xlim = c(0,100),
     col = c(rep('grey30', 4), rep('grey', 4), rep('grey30', 4), rep('grey', 4), rep('grey30', 9)))
axis(2, las = 2)
hist(metadata_h3n2$Host_Age, breaks = seq(0,125,5), main = 'H3N2', xlab = 'Age', ylab = 'Count', yaxt = 'n', xlim = c(0,100),
     col = c(rep('grey30', 4), rep('grey', 4), rep('grey30', 4), rep('grey', 4), rep('grey30', 9)))
axis(2, las = 2)
hist(metadata_RSVA$Patient.age, breaks = seq(0,125,5), main = 'RSV A', xlab = 'Age', ylab = 'Count', yaxt = 'n', xlim = c(0,100),
     col = c(rep('grey30', 4), rep('grey', 9), rep('grey30', 9)))
axis(2, las = 2)
hist(metadata_RSVB$Patient.age, breaks = seq(0,125,5), main = 'RSV B', xlab = 'Age', ylab = 'Count', yaxt = 'n', xlim = c(0,100),
     col = c(rep('grey30', 4), rep('grey', 9), rep('grey30', 9)))
axis(2, las = 2)


hists = hists_h1n1
plot(NULL, log = 'y', pch = 16, bty = 'n', yaxt = 'n', xaxt = 'n',
     xlim = c(0,5), ylim = c(0.5, 5), ylab = 'Ratio', xlab = 'Genetic distance', main = 'H1N1')
axis(1, at = seq(0.5, 4.5, 1), labels = paste0('[', 0:4, '-', 1:5, '['))
axis(2, at = c(0.5, 1, 2, 5), labels = c(0.5, 1, 2, 5), las = 2, tick = c(0.1, 0.5, 1, 5, 10, 50, 100))
abline(h = 1, lty = 2)
for(p in 1:length(hists$pair_labels)){
  tmp = lapply(1:hists$bootstrap, function(x)(((hists$pair_hists_boot$within_a[p,1:5,x]+1)/sum(hists$pair_hists_boot$within_a[p,1:5,x]+1) + (hists$pair_hists_boot$within_b[p,1:5,x]+1)/sum(hists$pair_hists_boot$within_b[p,1:5,x]+1))/2)/((hists$pair_hists_boot$between[p,1:5,x]+1)/sum((hists$pair_hists_boot$between[p,1:5,x]+1))))  
  tmp = do.call(rbind, tmp)
  tmp = apply(tmp, MAR = 2, median.and.ci)
  tmp[2,which(tmp[2,] == 0)] = 0.001
  tmp[3,which(is.infinite(tmp[3,]))] = 1000
  points((hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-length(hists$pair_labels)/2))/(length(hists$pair_labels)/2)*0.25 , tmp[1,], pch = 16, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
  arrows(x0 = (hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y0 = tmp[2,],
         x1 = (hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y1 = tmp[3,],
         length = 0.05, angle = 90, code = 0, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
}

hists = hists_h3n2
plot(NULL, log = 'y', pch = 16, bty = 'n', yaxt = 'n', xaxt = 'n',
     xlim = c(0,5), ylim = c(0.5, 5), ylab = 'Ratio', xlab = 'Genetic distance', main = 'H3N2')
axis(1, at = seq(0.5, 4.5, 1), labels = paste0('[', 0:4, '-', 1:5, '['))
axis(2, at = c(0.5, 1, 2, 5), labels = c(0.5, 1, 2, 5), las = 2, tick = c(0.1, 0.5, 1, 5, 10, 50, 100))
abline(h = 1, lty = 2)
for(p in 1:length(hists$pair_labels)){
  tmp = lapply(1:hists$bootstrap, function(x)(((hists$pair_hists_boot$within_a[p,1:5,x]+1)/sum(hists$pair_hists_boot$within_a[p,1:5,x]+1) + (hists$pair_hists_boot$within_b[p,1:5,x]+1)/sum(hists$pair_hists_boot$within_b[p,1:5,x]+1))/2)/((hists$pair_hists_boot$between[p,1:5,x]+1)/sum((hists$pair_hists_boot$between[p,1:5,x]+1))))  
  tmp = do.call(rbind, tmp)
  tmp = apply(tmp, MAR = 2, median.and.ci)
  tmp[2,which(tmp[2,] == 0)] = 0.001
  tmp[3,which(is.infinite(tmp[3,]))] = 1000
  points((hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-length(hists$pair_labels)/2))/(length(hists$pair_labels)/2)*0.25 , tmp[1,], pch = 16, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
  arrows(x0 = (hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y0 = tmp[2,],
         x1 = (hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y1 = tmp[3,],
         length = 0.05, angle = 90, code = 0, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
}
# legend('topright', col = 1:p, pch = 15, legend = hists$pair_labels, bty = 'n')

hists = hists_RSVA
plot(NULL, log = 'y', pch = 16, bty = 'n', yaxt = 'n', xaxt = 'n',
     xlim = c(0,5), ylim = c(0.25, 50), ylab = 'Ratio', xlab = 'Genetic distance', main = 'RSV A')
axis(1, at = seq(0.5, 4.5, 1), labels = paste0('[', 0:4, '-', 1:5, '['))
axis(2, at = c(0.25, 1, 10, 50), labels = c('<0.25', 1, 10, 50), las = 2, tick = c(0.25, 0.5, 1, 5, 10, 50, 100))
abline(h = 1, lty = 2)
for(p in 1:length(hists$pair_labels)){
  tmp = lapply(1:hists$bootstrap, function(x)(((hists$pair_hists_boot$within_a[p,1:5,x]+1)/sum(hists$pair_hists_boot$within_a[p,1:5,x]+1) + (hists$pair_hists_boot$within_b[p,1:5,x]+1)/sum(hists$pair_hists_boot$within_b[p,1:5,x]+1))/2)/((hists$pair_hists_boot$between[p,1:5,x]+1)/sum((hists$pair_hists_boot$between[p,1:5,x]+1))))  
  tmp = do.call(rbind, tmp)
  tmp = apply(tmp, MAR = 2, median.and.ci)
  tmp[2,which(tmp[2,] == 0)] = 0.001
  tmp[3,which(is.infinite(tmp[3,]))] = 1000
  points((hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-length(hists$pair_labels)/2))/(length(hists$pair_labels)/2)*0.25 , tmp[1,], pch = 16, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
  arrows(x0 = (hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y0 = tmp[2,],
         x1 = (hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y1 = tmp[3,],
         length = 0.05, angle = 90, code = 0, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
}
# legend('topright', col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))), pch = 15, legend = hists$pair_labels, bty = 'n')

hists = hists_RSVB
plot(NULL, log = 'y', pch = 16, bty = 'n', yaxt = 'n', xaxt = 'n',
     xlim = c(0,5), ylim = c(0.25, 50), ylab = 'Ratio', xlab = 'Genetic distance', main = 'RSV B')
axis(1, at = seq(0.5, 4.5, 1), labels = paste0('[', 0:4, '-', 1:5, '['))
axis(2, at = c(0.25, 1, 10, 50), labels = c('<0.25', 1, 10, 50), las = 2, tick = c(0.25, 0.5, 1, 5, 10, 50, 100))
abline(h = 1, lty = 2)
for(p in 1:length(hists$pair_labels)){
  tmp = lapply(1:hists$bootstrap, function(x)(((hists$pair_hists_boot$within_a[p,1:5,x]+1)/sum(hists$pair_hists_boot$within_a[p,1:5,x]+1) + (hists$pair_hists_boot$within_b[p,1:5,x]+1)/sum(hists$pair_hists_boot$within_b[p,1:5,x]+1))/2)/((hists$pair_hists_boot$between[p,1:5,x]+1)/sum((hists$pair_hists_boot$between[p,1:5,x]+1))))  
  tmp = do.call(rbind, tmp)
  tmp = apply(tmp, MAR = 2, median.and.ci)
  tmp[2,which(tmp[2,] == 0)] = 0.001
  tmp[3,which(is.infinite(tmp[3,]))] = 1000
  points((hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-length(hists$pair_labels)/2))/(length(hists$pair_labels)/2)*0.25 , tmp[1,], pch = 16, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
  arrows(x0 = (hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y0 = tmp[2,],
         x1 = (hists$common_breaks[1:5]+hists$common_breaks[2:6])/2 + ((p-(length(hists$pair_labels)/2))/(length(hists$pair_labels)/2))*0.25, y1 = tmp[3,],
         length = 0.05, angle = 90, code = 0, col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))[p]))
}
# legend('topright', col = adjustcolor(RColorBrewer::brewer.pal(name = 'Paired', n = length(hists$pair_labels))), pch = 15, legend = hists$pair_labels, bty = 'n')


hists = hists_h1n1
Q_contact_sym = Q_contact_sym_h1n1_socialcov
rates = fit_h1n1$chains$rates
ts = 365/(substitution_rate_h1n1*genome_length_h1n1)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(Q_contact_sym[lower.tri(Q_contact_sym)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(1e-3), max(1)), xlim = c(0.2, 5), xaxt = 'n',
     bty = 'n', yaxt = 'n', log = 'xy', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'Social contact rates in the population (per day)', ylab = 'Rates from genetic distances (per day)', main = 'H1N1')
axis(1, at = c(0.2, 0.5, 1, 2, 5))
axis(2, las = 2)
arrows(x0 = Q_contact_sym[lower.tri(Q_contact_sym)], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[lower.tri(Q_contact_sym)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = Q_contact_sym[lower.tri(Q_contact_sym)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-2, 2, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the age groups
mat = Q_contact_sym_h1n1_socialcov
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value
slopes_H1N1 = slopes
pvalue_H1N1 = length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)


hists = hists_h3n2
Q_contact_sym = Q_contact_sym_h3n2_socialcov
rates = fit_h3n2$chains$rates
ts = 365/(substitution_rate_h3n2*genome_length_h3n2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(Q_contact_sym[lower.tri(Q_contact_sym)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(0.001), max(1)), xlim = c(0.2, 5), xaxt = 'n',
     bty = 'n', yaxt = 'n', log = 'xy', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'Social contact rates in the population (per day)', ylab = 'Rates from genetic distances (per day)', main = 'H3N2')
axis(1, at = c(0.2, 0.5, 1, 2, 5))
axis(2, las = 2)
arrows(x0 = Q_contact_sym[lower.tri(Q_contact_sym)], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[lower.tri(Q_contact_sym)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = Q_contact_sym[lower.tri(Q_contact_sym)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-2, 2, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the age groups
mat = Q_contact_sym_h3n2_socialcov
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value
slopes_H3N2 = slopes
pvalue_H3N2 = length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)

hists = hists_RSVA
Q_contact_sym = Q_contact_sym_RSVA_socialcov
rates = fit_RSVA$chains$rates
ts = 365/(substitution_rate_RSVA*genome_length_RSVB)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(Q_contact_sym[lower.tri(Q_contact_sym)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(1e-3), max(1)), xlim = c(1.5, 5), xaxt = 'n',
     bty = 'n', yaxt = 'n', log = 'xy', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'Social contact rates in the population (per day)', ylab = 'Rates from genetic distances (per day)', main = 'RSV A')
axis(1, at = c(1.5, 3, 5))
axis(2, las = 2)
arrows(x0 = Q_contact_sym[lower.tri(Q_contact_sym)], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[lower.tri(Q_contact_sym)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = Q_contact_sym[lower.tri(Q_contact_sym)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 2, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the age groups
mat = Q_contact_sym_RSVA_socialcov
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value
slopes_RSVA = slopes
pvalue_RSVA = length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)


hists = hists_RSVB
Q_contact_sym = Q_contact_sym_RSVB_socialcov
rates = fit_RSVB$chains$rates
ts = 365/(substitution_rate_RSVB*genome_length_RSVB)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(Q_contact_sym[lower.tri(Q_contact_sym)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(1e-3), max(1)), xlim = c(1.5, 5), xaxt = 'n',
     bty = 'n', yaxt = 'n', log = 'xy', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'Social contact rates in the population (per day)', ylab = 'Rates from genetic distances (per day)', main = 'RSV B')
axis(1, at = c(1.5, 3, 5))
axis(2, las = 2)
arrows(x0 = Q_contact_sym[lower.tri(Q_contact_sym)], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[lower.tri(Q_contact_sym)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))

## Compute model fit
slopes = NULL
intercept = NULL
x = Q_contact_sym[lower.tri(Q_contact_sym)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 2, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the age groups
mat = Q_contact_sym_RSVB_socialcov
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value
slopes_RSVB = slopes
pvalue_RSVB = length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)
# dev.off()
#########################################################################################################################################


#########################################################################################################################################
## Slopes regressions
#########################################################################################################################################
median_and_ci(slopes_H1N1)
median_and_ci(slopes_H3N2)
median_and_ci(slopes_RSVA)
median_and_ci(slopes_RSVB)
pvalue_H1N1
pvalue_H3N2
pvalue_RSVA
pvalue_RSVB
#########################################################################################################################################


#########################################################################################################################################
## Figure SXX: Model fit
#########################################################################################################################################
## Observed vs predicted P same, P diff and ratio (distributions normalised over the fitted bins; median and 95% CrI) for each pair of age groups and fitted bin
plot_fits_virus = function(fit, virus) {
  fit_summary = plot_proportion_fits(fit, plot = F)
  p_names = which(lower.tri(diag(fit$stan_data$K))==T, arr.ind = T)
  p_names = paste0(p_names[,2], p_names[,1])
  pair_colors = c()
  for(p in 1:fit$stan_data$N_sets){
    pr = fit$stan_data$set_pair[p]
    p_id = paste0(sort(c(fit$stan_data$pair_a[pr], fit$stan_data$pair_b[pr])), collapse = '')
    p_id = match(p_id, p_names)
    pair_colors[p] = RColorBrewer::brewer.pal(name = 'Paired', n = fit$stan_data$N_sets)[p_id]
  }

  plot(NULL, xlim = c(10^-(2),1), ylim = c(10^-(2),1), bty = 'n', xaxt = 'n',  yaxt = 'n', log = 'xy', main = paste0('P same, ', virus), xlab = 'Observed', ylab = 'Predicted (median)')
  axis(1, at = 10^(seq(-2,0, 1)), labels = 10^(seq(-2,0, 1)))
  axis(2, las = 2, at = 10^(seq(-2,0, 1)), labels = 10^(seq(-2,0, 1)))
  abline(a=0, b=1, lty= 2)
  for(p in 1:fit$stan_data$N_sets){
    bins = which(fit_summary$fitted[p,])
    arrows(x0 = fit_summary$observed_same[p,bins], y0 = fit_summary$predicted_same_lo[p,bins], x1 = fit_summary$observed_same[p,bins], y1 = fit_summary$predicted_same_hi[p,bins], length = 0, col = pair_colors[p])
    points(fit_summary$observed_same[p,bins], fit_summary$predicted_same[p,bins], pch = 16, col = pair_colors[p])
  }

  plot(NULL, xlim = c(10^-(2),1), ylim = c(10^-(2),1), bty = 'n', xaxt = 'n',  yaxt = 'n', log = 'xy', main = paste0('P diff, ', virus), xlab = 'Observed', ylab = 'Predicted (median)')
  axis(1, at = 10^(seq(-2,0, 1)), labels = 10^(seq(-2,0, 1)))
  axis(2, las = 2, at = 10^(seq(-2,0, 1)), labels = 10^(seq(-2,0, 1)))
  abline(a=0, b=1, lty= 2)
  for(p in 1:fit$stan_data$N_sets){
    bins = which(fit_summary$fitted[p,])
    arrows(x0 = fit_summary$observed_diff[p,bins], y0 = fit_summary$predicted_diff_lo[p,bins], x1 = fit_summary$observed_diff[p,bins], y1 = fit_summary$predicted_diff_hi[p,bins], length = 0, col = pair_colors[p])
    points(fit_summary$observed_diff[p,bins], fit_summary$predicted_diff[p,bins], pch = 16, col = pair_colors[p])
  }

  plot(NULL, xlim = c(0.1,10), ylim = c(0.1,10), bty = 'n', xaxt = 'n',  yaxt = 'n', log = 'xy', main = paste0('Ratio, ', virus), xlab = 'Observed', ylab = 'Predicted (median)')
  axis(1, at = c(0.1, 1, 10), labels = c('<0.1', 1, 10), las = 1, tick = c(0.1, 0.5, 1, 5, 10, 50, 100))
  axis(2, at = c(0.1, 1, 10), labels = c('<0.1', 1, 10), las = 2, tick = c(0.1, 0.5, 1, 5, 10, 50, 100))
  abline(a=0, b=1, lty= 2)
  for(p in 1:fit$stan_data$N_sets){
    pr = fit$stan_data$set_pair[p]
    w = fit$stan_data$ref_weight[p,]
    bins = which(fit_summary$fitted[p,])
    ref = 0
    for(i in 1:fit$stan_data$N_boot){
      same = (w[1]*fit$stan_data$hist_within_a[pr,bins,i]/fit$stan_data$n_within_a[pr,i] + w[2]*fit$stan_data$hist_within_b[pr,bins,i]/fit$stan_data$n_within_b[pr,i])/(w[1]+w[2])
      ref = ref + ((1-w[3])*same + w[3]*fit$stan_data$hist_diff[pr,bins,i]/fit$stan_data$n_diff[pr,i])/fit$stan_data$N_boot
    }
    pred_same = t(apply(1 - w[3]*fit$chains$cr[,p,bins], MAR = 1, function(x)(x*ref)/sum(x*ref)))
    pred_diff = t(apply(fit$chains$cr[,p,bins], MAR = 1, function(x)(x*ref)/sum(x*ref)))
    ci = apply(pred_same/pred_diff, MARGIN = 2, FUN = median_and_ci)
    data = fit_summary$observed_same[p,bins]/fit_summary$observed_diff[p,bins]
    arrows(x0 = data, y0 = ci[2,], x1 = data, y1 = ci[3,], length = 0, col = pair_colors[p])
    points(data, ci[1,], pch = 16, col = pair_colors[p])
  }
}

# pdf(file = '5_age_groups/Figures/Figure_SXX_Age_group_fits_proportions_ratio_13102026.pdf', width = 5.7, height = 9)
par(mfrow = c(5,3), mar = c(2.5,2.5,2.5,2.5), oma = c(1,1,1,1))
plot_fits_virus(fit_h1n1, 'H1N1')
plot_fits_virus(fit_h3n2, 'H3N2')
plot_fits_virus(fit_RSVA, 'RSVA')
plot_fits_virus(fit_RSVB, 'RSVB')
plot_fits_virus(fit_SC2, 'SARS-CoV-2')
# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: Rates estimated vs N seqs
#########################################################################################################################################
# pdf(file = '5_age_groups/Figures/Figure_SXX_Estimated_rates_vs_Nseqs_13102026.pdf', 
    # width = 5, height = 6)
par(mfcol = c(3,2), mar = c(3,3,3,3), oma = c(1,1,1,1))

hists = hists_h1n1
N_seq = N_seqs_h1n1
rates = fit_h1n1$chains$rates
ts = 365/(substitution_rate_h1n1*genome_length_h1n1)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(N_seq[lower.tri(N_seq)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(1e-3), max(1)), 
     bty = 'n', yaxt = 'n', log = 'y', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'N seq', ylab = 'Rates from genetic distances (per day)', main = 'H1N1')
axis(2, las = 2)
arrows(x0 = N_seq[lower.tri(N_seq)], y0 = res_ci[2,]/ts, 
       x1 = N_seq[lower.tri(N_seq)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = N_seq[lower.tri(N_seq)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 7, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)


hists = hists_h3n2
N_seq = N_seqs_h3n2
rates = fit_h3n2$chains$rates
ts = 365/(substitution_rate_h3n2*genome_length_h3n2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(N_seq[lower.tri(N_seq)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(1e-3), max(1)), 
     bty = 'n', yaxt = 'n', log = 'y', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'N seq', ylab = 'Rates from genetic distances (per day)', main = 'H3N2')
axis(2, las = 2)
arrows(x0 = N_seq[lower.tri(N_seq)], y0 = res_ci[2,]/ts, 
       x1 = N_seq[lower.tri(N_seq)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = N_seq[lower.tri(N_seq)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 7, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)


hists = hists_SC2
N_seq = N_seqs_SC2
rates = fit_SC2$chains$rates
ts = 365/(substitution_rate_SC2*genome_length_SC2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(N_seq[lower.tri(N_seq)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(1e-3), max(1)), 
     bty = 'n', yaxt = 'n', log = 'y', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'N seq', ylab = 'Rates from genetic distances (per day)', main = 'SARS-CoV-2')
axis(2, las = 2)
arrows(x0 = N_seq[lower.tri(N_seq)], y0 = res_ci[2,]/ts, 
       x1 = N_seq[lower.tri(N_seq)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = N_seq[lower.tri(N_seq)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 7, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)


hists = hists_RSVA
N_seq = N_seqs_RSVA
rates = fit_RSVA$chains$rates
ts = 365/(substitution_rate_RSVA*genome_length_RSVA)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(N_seq[lower.tri(N_seq)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(1e-3), max(1)), 
     bty = 'n', yaxt = 'n', log = 'y', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'N seq', ylab = 'Rates from genetic distances (per day)', main = 'RSV A')
axis(2, las = 2)
arrows(x0 = N_seq[lower.tri(N_seq)], y0 = res_ci[2,]/ts, 
       x1 = N_seq[lower.tri(N_seq)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = N_seq[lower.tri(N_seq)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 7, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

hists = hists_RSVB
N_seq = N_seqs_RSVB
rates = fit_RSVB$chains$rates
ts = 365/(substitution_rate_RSVB*genome_length_RSVB)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(N_seq[lower.tri(N_seq)], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(1e-3), max(1)), 
     bty = 'n', yaxt = 'n', log = 'y', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'N seq', ylab = 'Rates from genetic distances (per day)', main = 'RSV B')
axis(2, las = 2)
arrows(x0 = N_seq[lower.tri(N_seq)], y0 = res_ci[2,]/ts, 
       x1 = N_seq[lower.tri(N_seq)], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = N_seq[lower.tri(N_seq)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 7, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: SC2 rates estimated vs Beraud
#########################################################################################################################################
# pdf(file = '5_age_groups/Figures/Figure_SXX_Estimated_rates_vs_Beraud_13102026.pdf', 
    # width = 6, height = 3)
par(mfcol = c(1,2), mar = c(2,2,2,2), oma = c(0,0,1,0))

hists = hists_SC2
Q_contact_sym = Q_contact_sym_SC2
Q_contact_sym[upper.tri(Q_contact_sym)] = NA
rates = fit_SC2$chains$rates
ts = 365/(substitution_rate_SC2*genome_length_SC2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
idx = which(!is.na(Q_contact_sym), arr.ind = T)
plot(Q_contact_sym[idx], res_ci[1,]/ts, pch = 16, 
     ylim = c(min(0.001), max(1)), main = 'Beraud', xlim = c(0.2, 8), xaxt = 'n',
     bty = 'n', yaxt = 'n', log = 'xy', col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)),
     xlab = 'Social contact rates in the population (per day)', ylab = 'Rates from genetic distances (per day)')
axis(1, at = c(0.2, 0.5, 1,3, 8))
axis(2, las = 2, at = c(0.001, 0.01, 0.1, 1))
arrows(x0 = Q_contact_sym[idx], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[idx], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = Q_contact_sym[lower.tri(Q_contact_sym)]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-2, 3, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the age groups
mat = Q_contact_sym_SC2
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value

# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: SC2 rates estimated vs asymetric
#########################################################################################################################################
# pdf(file = '5_age_groups/Figures/Figure_SXX_Estimated_rates_vs_Asym_13102026.pdf', 
    # width = 6, height = 3)
par(mfcol = c(1,2), mar = c(2,2,2,2), oma = c(0,0,1,0))

hists = hists_SC2
Q_contact_sym = t(Q_contact_asym_SC2)
rates = fit_SC2$chains$rates
ts = 365/(substitution_rate_SC2*genome_length_SC2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(NULL,  ylim = c(min(0.001), max(1)), xlim = c(min(Q_contact_sym, na.rm = T), max(Q_contact_sym, na.rm = T)), main = 'Asym Beraud',  bty = 'n', yaxt = 'n', log = 'xy', 
     xlab = 'Social contact rates in the population (per day)', ylab = 'Rates from genetic distances (per day)')
axis(2, las = 2, at = c(0.001, 0.01, 0.1, 1))
idx = which(lower.tri(Q_contact_sym), arr.ind = T)
points(Q_contact_sym[idx], res_ci[1,]/ts, pch = 16, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
arrows(x0 = Q_contact_sym[idx], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[idx], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
Q_contact_sym = t(Q_contact_sym)
points(Q_contact_sym[idx], res_ci[1,]/ts, pch = 16, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
arrows(x0 = Q_contact_sym[idx], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[idx], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))

## Compute model fit
slopes = NULL
intercept = NULL
x = c(t(Q_contact_asym_SC2)[which(lower.tri(Q_contact_asym_SC2), arr.ind = T)], (Q_contact_asym_SC2)[which(lower.tri(Q_contact_asym_SC2), arr.ind = T)])
for(ii in 1:nrow(rates)){
  y =  rep(rates[ii,]/ts, 2)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-3, 3, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the age groups
mat = t(Q_contact_asym_SC2)
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(cbind(rates[draws,], rates[draws,])/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  x = c(t(mat_tmp)[which(lower.tri(mat_tmp), arr.ind = T)], (mat_tmp)[which(lower.tri(mat_tmp), arr.ind = T)])
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value



hists = hists_SC2
Q_contact_sym = t(Q_contact_asym_SC2_socialcov)
rates = fit_SC2$chains$rates
ts = 365/(substitution_rate_SC2*genome_length_SC2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
plot(NULL,  ylim = c(min(0.001), max(1)), xlim = c(min(Q_contact_sym, na.rm = T), max(Q_contact_sym, na.rm = T)), main = 'Asym SocialCov',  bty = 'n', yaxt = 'n', log = 'xy', 
     xlab = 'Social contact rates in the population (per day)', ylab = 'Rates from genetic distances (per day)')
axis(2, las = 2, at = c(0.001, 0.01, 0.1, 1))
idx = which(lower.tri(Q_contact_sym), arr.ind = T)
points(Q_contact_sym[idx], res_ci[1,]/ts, pch = 16, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
arrows(x0 = Q_contact_sym[idx], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[idx], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
Q_contact_sym = t(Q_contact_sym)
points(Q_contact_sym[idx], res_ci[1,]/ts, pch = 16, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
arrows(x0 = Q_contact_sym[idx], y0 = res_ci[2,]/ts, 
       x1 = Q_contact_sym[idx], y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = RColorBrewer::brewer.pal(name = 'Paired', n = ncol(res_ci)))
## Compute model fit
slopes = NULL
intercept = NULL
x = c(t(Q_contact_asym_SC2_socialcov)[which(lower.tri(Q_contact_asym_SC2_socialcov), arr.ind = T)], (Q_contact_asym_SC2_socialcov)[which(lower.tri(Q_contact_asym_SC2_socialcov), arr.ind = T)])
for(ii in 1:nrow(rates)){
  y =  rep(rates[ii,]/ts, 2)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-3, 3, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)
## Compute p-value: median slope over draws, for every relabelling of the age groups
mat = t(Q_contact_asym_SC2_socialcov)
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(cbind(rates[draws,], rates[draws,])/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  x = c(t(mat_tmp)[which(lower.tri(mat_tmp), arr.ind = T)], (mat_tmp)[which(lower.tri(mat_tmp), arr.ind = T)])
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value

# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Trees colored by age group
#########################################################################################################################################
# pdf(file = '5_age_groups/Figures/FigureSXX_Tree_influenza_age.pdf', width = 19/2.56, height = 25/2.56)
par(mfrow = c(2,1), mar = c(0,0,2,0))
metadata_h1n1_tmp = metadata_h1n1[which(!is.na(metadata_h1n1$Host_Age)),]
tree_h1n1 = ladderize(tree_h1n1)
age_group_label = cut(as.numeric(metadata_h1n1_tmp$Host_Age), age_breaks_h1n1)
levels(age_group_label) = RColorBrewer::brewer.pal(name = 'Paired', n = 6)
age_group_label = as.character(age_group_label)
plot.phylo(tree_h1n1, show.tip.label = F)
tiplabels(pch = 16, cex = 1, col = adjustcolor(age_group_label, alpha.f = 1))
add.scale.bar()
legend('topright', pch = 16, legend  = age_breaks_h1n1, col = RColorBrewer::brewer.pal(name = 'Paired', n = 6))

metadata_h3n2_tmp = metadata_h3n2[which(!is.na(metadata_h3n2$Host_Age)),]
tree_h3n2 = ladderize(tree_h3n2)
age_group_label = cut(as.numeric(metadata_h3n2_tmp$Host_Age), age_breaks_h3n2)
levels(age_group_label) = RColorBrewer::brewer.pal(name = 'Paired', n = 6)
age_group_label = as.character(age_group_label)
plot.phylo(tree_h3n2, show.tip.label = F)
tiplabels(pch = 16, cex = 1, col = adjustcolor(age_group_label, alpha.f = 1))
add.scale.bar()
# dev.off()

# pdf(file = '5_age_groups/Figures/FigureSXX_Tree_RSV_age.pdf', width = 19/2.56, height = 25/2.56)
par(mfrow = c(2,1), mar = c(0,0,2,0))
metadata_RSVA_tmp = metadata_RSVA[which(!is.na(metadata_RSVA$Patient.age)),]
tree_RSVA = ladderize(tree_RSVA)
age_group_label = cut(as.numeric(metadata_RSVA_tmp$Patient.age), age_breaks_RSVA)
levels(age_group_label) = RColorBrewer::brewer.pal(name = 'Paired', n = 3)
age_group_label = as.character(age_group_label)
plot.phylo(tree_RSVA, show.tip.label = F)
tiplabels(pch = 16, cex = 1, col = adjustcolor(age_group_label, alpha.f = 1))
add.scale.bar()

metadata_RSVB_tmp = metadata_RSVB[which(!is.na(metadata_RSVB$Patient.age)),]
tree_RSVB = ladderize(tree_RSVB)
age_group_label = cut(as.numeric(metadata_RSVB_tmp$Patient.age), age_breaks_RSVB)
levels(age_group_label) = RColorBrewer::brewer.pal(name = 'Paired', n = 3)
age_group_label = as.character(age_group_label)
plot.phylo(tree_RSVB, show.tip.label = F)
tiplabels(pch = 16, cex = 1, col = adjustcolor(age_group_label, alpha.f = 1))
add.scale.bar()
# dev.off()

# pdf(file = '5_age_groups/Figures/FigureSXX_Tree_SC2_age.pdf', width = 19/2.56, height = 25/2.56)
par(mfrow = c(2,1), mar = c(0,0,2,0))
metadata_SC2_tmp = metadata_SC2[which(!is.na(metadata_SC2$Patient.age)),]
idx = match(unlist(lapply(tree_SC2$tip.label, function(x)str_split(x, '\\|')[[1]][2])), metadata_SC2_tmp$Accession.ID)
idx = which(is.na(idx))
if(length(idx) > 1){tree_SC2 = drop.tip(tree_SC2, tree_SC2$tip.label[idx])}
tree_SC2 = ladderize(tree_SC2)
age_group_label = cut(as.numeric(metadata_SC2_tmp$Patient.age), age_breaks_SC2)
levels(age_group_label) = RColorBrewer::brewer.pal(name = 'Paired', n = 6)
age_group_label = as.character(age_group_label)
plot.phylo(tree_SC2, show.tip.label = F)
tiplabels(pch = 16, cex = 1, col = adjustcolor(age_group_label, alpha.f = 1))
add.scale.bar()
# dev.off()
#########################################################################################################################################







