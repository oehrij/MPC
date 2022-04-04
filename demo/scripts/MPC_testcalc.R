##########################################################
#####
##### project:       multispecies connectivity modelling 
##### title:         MPC_testcalc
##### description:   testing metapopulation capacity function calculation with shapefile, raster and patch area vector/distance matrix input
##### author:        Jacqueline Oehri (JO)
##### date:          30.03.2022
##### last modified: 31.03.2022
##### dependencies:  MPC_genexdata.R output
#########################################################

##########################################################
### clean space
rm(list=ls(all=TRUE))

########################################################
## libraries
library(sf)               # new spatial vector processing package
library(raster)           # old (but still useful) spatial raster processing
library(stars)            # spatial processing (conversion between sf and raster possible, spatiotemporal datacube analyses)
library(metacapa)         # for metapopulation capacity according to Strimas-Mackey & Brodie 2018 

##########################################################
## set directories 
## load MPC function
# dir = getwd()
# source(paste(dir,"R/MPC_functions.R",sep="/"))
library(devtools);
load_all("."); # Working directory should be in the package MPC package

##########################################################
## Start
##########################################################

########################################################
###  Test with raster and shapefile (generated in MPC_genexdata_JO.R)
shpf = st_read(paste(dir,"demo/data/MPC_patchex.shp",sep="/"))
rasf = raster(paste(dir,"demo/data/MPC_patchex.tif",sep="/"))
MPC_fun(x=rasf)
MPC_fun(x=shpf) 
# $mpc
# [1] 369482018
# $pimport
# [1] 0.09928649 0.27038876 0.08307504 0.18079092 0.36645879
# $pid
# [1] "1" "4" "5" "2" "3"
### works!

##########################################################
## Test with pa and mdist data only
##
## Patch Areas - parameter "pa" (unit here=m2)
dfex   = data.frame(patchID=c(1,2,3),aream2=c(500,600,700))
pa     = dfex$aream2
## Pairwise patch distances - parameter "mdist"
mdist  = matrix(ncol=3,nrow=3,byrow=TRUE,c(0.00,150.3465,300.2232,
                                           150.3465,0.00,200.05,
                                           300.2232,200.05,0.00))

## Calculate MPC 
print(MPC_fun(pa=pa,mdist=mdist,dispfop="negex",alpha=317)$mpc)   # JO version
#[1] 30255.36

## Comparison with Strimas-Mackey: define the dispersal function and alpha 
alpha   = 317
dispfun = dispfun= function(d) { exp(-(1/alpha) * abs(d)) }
print(metacapa::meta_capacity(x=mdist,a=pa,f=dispfun,ex=0.5,self=TRUE,patch_mc=TRUE)$capacity) # Strimas-Mackey version
#[1] 30255.36

####
## More and larger areas and distances are now in square meters! (to make sure patch areas and distances are correct,they were originally derived from 9e+06, see below!!)
dfex   = data.frame(patchID=c(1:5),aream2=c(306900, 321300, 175500, 191700, 335700))
pa     = dfex$aream2
mdist = matrix(ncol=5,nrow=5,byrow=TRUE,
               c(0.0000, 1201.4991, 1991.1052, 742.7651, 1209.3387,
                 1201.4991,0.0000, 150.0000, 660.6815, 402.4922,
                 1991.1052, 150.0000, 0.0000, 1290.3488, 792.0227,
                 742.7651, 660.6815, 1290.3488, 0.0000, 30.0000,
                 1209.3387, 402.4922, 792.0227, 30.0000, 0.0000))

colnames(mdist) = c(1:5)
rownames(mdist) = c(1:5)

## Calculate MPC 
print(MPC_fun(pa=pa,mdist=mdist,dispfop="negex",alpha=317)$mpc)   # JO version                                           # JO version
## Comparison with Strimas-Mackey: define the dispersal function and alpha 
alpha   = 317
dispfun = dispfun= function(d) { exp(-(1/alpha) * abs(d)) }
print(metacapa::meta_capacity(x=mdist,a=pa,f=dispfun,ex=0.5,self=TRUE,patch_mc=TRUE)$capacity) # Strimas-Mackey version

####
## One maximum area
dfamax   = data.frame(patchID=1,aream2=9e+06)
pamax    = dfamax$aream2
mdistmax = matrix(ncol=1,nrow=1,byrow=TRUE,c(0.0000))
## Calculate MPC 
print(MPC_fun(pa=pamax,mdist=mdistmax,dispfop="negex",alpha=317)$mpc)   # JO version                                           # JO version
## Comparison with Strimas-Mackey: define the dispersal function and alpha 
alpha   = 317
dispfun = dispfun= function(d) { exp(-(1/alpha) * abs(d)) }
print(metacapa::meta_capacity(x=mdistmax,a=pamax,f=dispfun,ex=0.5,self=TRUE,patch_mc=TRUE)$capacity) # Strimas-Mackey version
#[1] 2.7e+10


########################################################
### End
########################################################





