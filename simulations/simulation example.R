require(terra)

#download the network data from Dalsgaard et el (2021)
#B. Dalsgaard et al., The influence of biogeographical and evolutionary histories on morphological trait-matching and resource specialization in mutualistic hummingbird-plant networks. Functional Ecology (2021). https://doi.org/10.1111/1365-2435.13784.
path="...Downloads/doi_10.5061_dryad.rr4xgxd7n__v3-2 2"

# set working directory to the github folder
setwd("...Downloads/hum-ex-coex-col-main 3")

#ids of central American networks in the database
C_Am_nets=c(4 , 5 , 6,  7 , 8,  9, 10, 11, 12, 21, 22, 23, 24, 25, 79, 80, 87, 88, 89, 90)
#list of temperate migratory hummingbirds occasionally observed in Central American networks
migra=read.csv("data/migratory hummingbirds.csv",sep=";",h=T)

#loading network metadata
dat=read.csv("data/Community data.csv",h=T,sep=";")

#loading climate PCA data
load("data/pca.RData")

#loading data on the hummingbird species' contemporary climate volume climate  
load("data/sp_clim.RData")

#loading contemporary climate variables
clim_PCA1_reg=rast("data/clim_PCA1_reg.tif")
clim_PCA2_reg=rast("data/clim_PCA2_reg.tif")

# downloading future climate scenario (example using ACCESS1 model) from the CHELSA database (https://chelsa-climate.org)
fc1=rast("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/cmip5/2061-2080/bio/CHELSA_bio_mon_ACCESS1-0_rcp45_r1i1p1_g025.nc_1_2061-2080_V1.2.tif")
fc4=rast("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/cmip5/2061-2080/bio/CHELSA_bio_mon_ACCESS1-0_rcp45_r1i1p1_g025.nc_1_2061-2080_V1.2.tif")
fc12=rast("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/cmip5/2061-2080/bio/CHELSA_bio_mon_ACCESS1-0_rcp45_r1i1p1_g025.nc_12_2061-2080_V1.2.tif")
fc15=rast("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/cmip5/2061-2080/bio/CHELSA_bio_mon_ACCESS1-0_rcp45_r1i1p1_g025.nc_15_2061-2080_V1.2.tif")
fut_ras=rast(list(fc1,fc4,fc12,fc15))
fut_ras=crop(fut_ras,clim_PCA1_reg)

#loading simulation algorithms
source("simulations/Rewiring function.R") # function for rewireing
source("simulations/Simulating extinctions.R") # function for climate-driven extinction and coextinctions
source("simulations/Simulating extinctions with rewiring.R") # function for climate-driven extinction and coextinctions with rewireing
source("simulations/Simulating colonization.R") # function for species colonization

#running extinction-coextinction simulations for a brazilian humingbird plant network network id91
res_without_rewireing=coex_model("91_Brazil.csv",network_folder = path,n_reps = 500)
ex=res_without_rewireing[[1]] #proportion of climate-driven extinctions
ex_and_coex=res_without_rewireing[[2]] # proportion of climate-driven extinctions + coextinctions


#running extinction-coextinction with 0.5 rewireing of interactions
res_with_rewireing=coex_model_rewireing("91_Brazil.csv",network_folder = path,n_reps = 500,r_flex = 0.5)
ex_rw=res_with_rewireing[[1]] #proportion of climate-driven extinctions
ex_and_coex_rw=res_with_rewireing[[2]]

#running colonization simulations 
col_res=col_model("91_Brazil.csv",network_folder = path,n_reps = 500)
col_res$col_net # number of colonizing species relative to the number of species in the network
col_res$col_net # number of colonizing species relative to the number of species in the network's source pool
