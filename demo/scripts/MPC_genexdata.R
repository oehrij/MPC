##########################################################
#####
##### project:       multispecies connectivity modelling
##### title:         MPC_genexdata
##### description:   generate example data (a raster and a shapefile with n patches) for MPC_testcalc_JO.R script
##### author:        Jacqueline Oehri (JO)
##### date:          30.03.2022
##### last modified: 31.03.2022
#########################################################

##########################################################
### clean space
rm(list=ls(all=TRUE))

##########################################################
## libraries
library(pascal)           # library(devtools); install_github("pascal-niklaus/pascal/pascal")
library(sf)               # new spatial vector processing package
library(raster)           # old (but still useful) spatial raster processing
library(stars)            # spatial processing (conversion between sf and raster possible, spatiotemporal datacube analyses)
library(dismo)            # species distribution modeling, create circles, circular polygons..
library(metacapa)         # for metapopulation capacity according to Strimas-Mackey & Brodie 2018

##########################################################
## set directory
dir = getwd()

##########################################################
## Start
##########################################################

##########################################################
## extent around Montréal example
MTLex  = extent(1674150,1677150,855840,858840)  # larger extent could also be available: #MTLex = extent(1674150,1704150,855840,885840)
ri     = raster(MTLex, res = 30)
## set EPSG Albers Equal Area projection: this is the old proj.4 string in the future it should be replaced..e.g. with similar notions such as # "EPSG:9822" # st_crs(9822)$wkt
crs(ri) = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=44.75 +lat_2=55.75 +x_0=0 +y_0=0 +a=6378137 +rf=298.257223999999 +units=m +no_defs"

## create some patches of habitat in the landscape
set.seed(16)
## make actual coordinates of circle centers
npatch    = 5    # select number of patches
radius    = 300  # select mean radius of patches
sdradius  = 50   # select sd of radius of patches

## sample locations for circle centers in the virtual landscape
rsampl = sample(c(1:ncell(ri)),npatch,replace=FALSE)
## set all values in raster NA except circle centers
ri[] =NA; ri[rsampl] =1
## plot the circle centers
par(mfrow=c(1,1),oma=c(1,1,1,1),mar=c(4,4,3,3),xpd=NA); plot(ri,col=c(NA,"black")); title(main=sprintf("npatch%d",npatch),line=1)
## get the coordinates of the circle centers
rsamplxy <- xyFromCell(ri,cell=rsampl)
## extract Ns and Es of circle centers
Ns <- rsamplxy[,2]
Es <- rsamplxy[,1]
## make a spatial points dataframe
spg1 <- data.frame(N=Ns,E=Es)
coordinates(spg1) <- ~ E + N
crs(spg1) <- crs(ri)
## plot for better overview
plot(spg1,add=TRUE,col="red",lwd=2)

## actually make circles with different radii
plist    = list()
radiuss1 = c()
index    = 1
for(ii in 1:length(spg1)) {
  print(paste("circle --- ",ii,sep=""))
  radiuss = rnorm(1,mean=radius,sd=sdradius)   # add the SD among radii
  circ <- circles(spg1[ii], d =radiuss,n=3600) # if I increase the n (default 360, the circles are more precise later?) ## add some sd here!!
  c    <- circ@polygons@polygons
  c1   <- c[[1]]@Polygons[[1]]
  coords <- slot(c1,"coords")
  polyA <- list(x=coords[,1],y=coords[,2])
  p  = Polygon(as.matrix(cbind(polyA$x,polyA$y)))
  ps = Polygons(list(p),radius+index)
  #slot(ps,"ID") <- ii
  plist[[index]] <- ps
  radiuss1 =c(radiuss1,radiuss)
  index = index + 1
} # spg1
## make a spatial polygon
sps = SpatialPolygons(plist)
proj4string(sps) <- crs(ri)
plot(sps,add=TRUE,border="red",lwd=2)
## just in case: make an sf object out of it
sps2         <- st_as_sf(sps)
### make sure crs is conserved
st_crs(sps2) <- st_crs(sps)

## make a dataframe and add it to  simple feature (sf) (for sure can be done more beautifully..)
## actually this area is much more exact than the st_area that is calculated afterwards...
IDs1       <- sapply(sps@polygons,function(x){slot(x,"ID")})
areas1     <- sapply(sps@polygons,function(x){slot(x,"area")})
centroids1 <- sapply(sps@polygons,function(x){slot(x,"labpt")})
dfr  <- data.frame(ID = IDs1,area=areas1,centrN=centroids1[2,],centrE=centroids1[1,],radius=radiuss1) # note how we use the same IDs from above!
sps3 <- st_sf(dfr , geometry=sps2$geometry)
## make sure crs is preserved
st_crs(sps3) <- st_crs(sps)

####
## write as shapefile
if(file.exists(paste(dir,"demo/data/MPC_circlex.shp",sep="/"))) { file.remove(paste(dir,"demo/data/MPC_circlex.shp",sep="/"))}
st_write(sps3, paste(dir,"demo/data/MPC_circlex.shp",sep="/"), driver="ESRI Shapefile")  # create to a shapefile

### make a raster clump where circles are!
ri4 = ri
ri4[] =NA
## extract raster cell identities that are covered by circles
vals = unique(unlist(raster::extract(ri4, as_Spatial(sps3$geometry),method="simple",cellnumbers=TRUE)))
vals = vals[is.finite(vals)]
ri4[vals] = 1

## make sure that shapefile and raster are EXACTLY the same area (i.e. along squares and not idealized circles like above...)
lclump   <- landscapemetrics::get_patches(ri4, directions=4, class=1)[[1]][[1]]
x <- st_as_stars(lclump) %>% st_as_sf(merge = TRUE)
x <- st_make_valid(x,reason=TRUE)
## make sure crs is preserved
# compareCRS(x,sps)
# [1] TRUE
# compareCRS(lclump,sps)
# [1] TRUE

plot(lclump,add=TRUE) # works! ri4 is now the new habitat raster
plot(x$geometry,border="magenta",add=TRUE,lwd=2)

####
## write as shapefile
if(file.exists(paste(dir,"demo/data/MPC_patchex.shp",sep="/"))) { file.remove(paste(dir,"demo/data/MPC_patchex.shp",sep="/"))}
st_write(x, paste(dir,"demo/data/MPC_patchex.shp",sep="/"), driver="ESRI Shapefile")  # create to a shapefile
####
## write as raster file
writeRaster(lclump, paste(dir,"demo/data/MPC_patchex.tif",sep="/"), format="GTiff", overwrite=TRUE) # Create a geoTiff file
####

########################################################
### End
########################################################





