## Load pairate functions
path_pairate = paste0(getwd(), '/pairate/')
for(f in list.files(paste0(path_pairate, 'R/'), full.names = T)) source(f)
stan_file_path = function(filename) return(paste0(path_pairate, 'inst/stan_models/', filename))
