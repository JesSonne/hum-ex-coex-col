# hum-ex-coex-col
 Supplementary data for the paper: "Extinction, coextinction, and colonization affecting plant-hummingbird networks under climate change"

Sp_clim.RData: a named list of 356 data.frames (i.e. one for each extant hummingbird species). Each data.frame contains two principal component axes used to characterize a species climate volume. The principal componant analysis includes the variation in mean annual temperature, mean annual precipitation, annual temperature seasonality (standard deviation*100), and anual precipitation seasonality (coefficient of variation). All climate data were downloaded from the CHELSA climate database (https://chelsa-climate.org)

pca.RData: prcomp data object containing the principal components analysis perfomred on the contemporary climate data

source_pool_list_10km.RData: a named list of 84 character strings (i.e. one for each hummingbird-plant community). Each character string lists the hummingbird species occuring within a 10 km radius surrounding a focal network.

source_pool_list_100km.RData: a named list of 84 character strings (i.e. one for each hummingbird-plant community). Each character string lists the hummingbird species occuring within a 100 km radius surrounding a focal network.
