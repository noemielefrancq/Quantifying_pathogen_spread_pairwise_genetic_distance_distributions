###########################################################################################
## Get daily passengers between pairs of countries (NCBI analysis, worldwide countries)
###########################################################################################
data_GTMD2 = read.csv('4_SC2_geography/Data_countries/Mobility_between_countries/GTMD2_trips.csv')
data_GTMD2 = data_GTMD2[which(data_GTMD2$year == 2022),]

## Population sizes
pops = read.csv('4_SC2_geography/Data_countries/Mobility_between_countries/API_SP.POP.TOTL_DS2_en_csv_v2_451027.csv', header = F)  ## https://data.worldbank.org/indicator/SP.POP.TOTL ## https://data.worldbank.org/indicator/SP.POP.TOTL
pops = pops[-c(1,2),]
colnames(pops) = pops[1,]
pops = pops[-c(1),]
m = match(c('DNK', 'JPN', 'SVK', 'CHE', 'THA', 'GBR', 'USA'), pops$`Country Code`)
pop_sizes = pops$`2022`[m]
names(pop_sizes) = c('DK', 'JP', 'SK', 'CH', 'TH', 'GB', 'US')

## Fill matrix
mobility_flow_mat = matrix(0, nrow = 7, ncol = 7)
colnames(mobility_flow_mat) = rownames(mobility_flow_mat) = names(pop_sizes) 
for(i in 1:length(names(pop_sizes))){
  for(j in 1:length(names(pop_sizes))){
    if(i != j){
      idx = which(data_GTMD2$iso2code_i == names(pop_sizes)[i] & data_GTMD2$iso2code_j == names(pop_sizes)[j])
      if(length(idx)>0){
        mobility_flow_mat[i,j] = data_GTMD2$gtmd2_trips_s2[idx]/pop_sizes[names(pop_sizes)[i]]
      }
    }
  }
}
diag(mobility_flow_mat) = NA

## Build symmetric matrix
mobility_flow_mat_sym = mobility_flow_mat
mobility_flow_mat_sym = (mobility_flow_mat + t(mobility_flow_mat))/2
heatmap(log(mobility_flow_mat_sym), scale = 'none', Rowv = NA, Colv = NA)

## Save
saveRDS(file = '4_SC2_geography/Data_countries/Mobility_between_countries/Mobility_data_between_countries_ncbi.rds', object = mobility_flow_mat)
saveRDS(file = '4_SC2_geography/Data_countries/Mobility_between_countries/Mobility_data_between_countries_symmetric_ncbi.rds', object = mobility_flow_mat_sym)
############################################################################



