## Description of data files

Sp_clim.RData: a named list of 356 data.frames (i.e. one for each extant hummingbird species). Each data.frame contains two principal component axes used to characterize a species climate volume. The principal component analysis includes the variation in mean annual temperature, mean annual precipitation, annual temperature seasonality (standard deviation * 100), and annual precipitation seasonality (coefficient of variation). All climate data were downloaded from the CHELSA climate database (https://chelsa-climate.org; Karger et al., 2017).


**pca.RData**: prcomp data object containing the principal components analysis performed on the contemporary climate data

**source_pool_list_10km.RData**: a named list of 84 character strings (i.e. one for each hummingbird-plant community). Each character string lists the hummingbird species occurring within a 10 km radius surrounding a focal network.

**source_pool_list_100km.RData**: a named list of 84 character strings (i.e. one for each hummingbird-plant community). Each character string lists the hummingbird species occurring within a 100 km radius surrounding a focal network.

**Species-level centrality measures:** data frame listing names, geographical coordinates, and biogeographical region for each hummingbird-plant community. Cross-referenced with the network dataset published in Dalsgaard et al. (2021)

### References
**B. Dalsgaard et al.,** The influence of biogeographical and evolutionary histories on morphological trait-matching and resource specialization in mutualistic hummingbird-plant networks. Functional Ecology (2021). https://doi.org/10.1111/1365-2435.13784.

**D. N. Karger et al.,** Climatologies at high resolution for the earth’s land surface areas. Scientific Data 4, 170122 (2017). http://dx.doi.org/doi:10.5061/dryad.kd1d4


