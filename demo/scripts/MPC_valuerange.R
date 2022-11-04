##########################################################
#####
##### project:       multispecies connectivity modelling
##### title:         MPC_valuerange
##### description:   test metapopulation capacity (MPC) values in dependence of area and distance
##### author:        Jacqueline Oehri (JO)
##### date:          29.03.2022
##### last modified: 07.04.2022 # changed library loading
##### comments:      -question: if patch areas are in km2 - can distances still be counted in metres? # actually - they should also be noted in km??
#####                -this script can be written much more efficiently..
#########################################################

##########################################################
### clean space
rm(list=ls(all=TRUE))

##########################################################
## libraries
library(sf)               # new spatial vector processing package
library(raster)           # old (but still useful) spatial raster processing
library(terra)            # new spatial raster processing
library(stars)            # spatial processing (conversion between sf and raster possible, spatiotemporal datacube analyses)
library(dismo)            # species distribution modeling, create circles, circular polygons..
library(metacapa)         # for metapopulation capacity according to Strimas-Mackey & Brodie 2018
library(lconnect)         # for landscape connectivity metrics

##########################################################
## set directories
## load MPC function
## 1) simplest way
dir = getwd()
# source(paste(dir,"R/MPC_functions.R",sep="/"))
## 2) local way
# library(devtools)
# load_all(".")      # Working directory should be in the package MPC package
## 3) more elegant local way (install it directly)
#devtools::install()
library(MPC)         # Working directory should be in the package MPC package
## 4) from github (soon available)
# library(devtools)
# install_github("oehrij/MPC")


##########################################################
## Start
##########################################################

##########################################################
## Value range analysis: change in metapopulation capacity across a range of habitat areas, when whole landscape extent is habitat
## test ranges in km^2 (40000 km^ ~ size of Switzerland)
avec = c(0,seq(1,40001,by=10))
## make empty dataframe for adding calculated MPC values
mpdf = data.frame(pa=numeric(0),mpc1=numeric(0),mpc2=numeric(0),mpc3=numeric(0))
## loop through avec for MPC calculation
for(pa in avec) {
  print(pa)
  mdist = matrix(ncol=1,nrow=1,byrow=TRUE,0)
  # parameters are loaded above
  mpc1=MPC_fun(pa=pa,mdist=mdist,alpha=317,dispfop="negex")$mpc
  mpc2=MPC_fun(pa=pa,mdist=mdist,alpha=317,dispfop="linear")$mpc
  mpc3=MPC_fun(pa=pa,mdist=mdist,alpha=317,dispfop="log-sech")$mpc
  mpdf0=data.frame(pa=pa,mpc1=mpc1,mpc2=mpc2,mpc3=mpc3)
  mpdf=rbind(mpdf,mpdf0)
}

## write to file
write.csv(mpdf,paste(dir,"demo/data/MPC_value_range.csv",sep="/"),row.names = FALSE)

### make a plot
pdf(paste(dir,"demo/data/MPC_value_range.pdf",sep="/"),width=6,height=4)
sink(paste(dir,"demo/data/MPC_value_range.txt",sep="/"))
## preparations
par(mfrow=c(1,3))
plotcols = c("black","olivedrab3","coral")
print("---summary of MPC in one patch---")
print(summary(mpdf))
print("---remove the case where pa and mpc's = 0 for model fit---")
mpdf=mpdf[-1,]
# fit a linear model to estimate coefficients! since all mpcs (from different dispun's) are the same because only 1 patch here, do it only with 1
print("---relationship almost linear in log-log space---")
lm2  = lm(log(mpdf$mpc1) ~ log(mpdf$pa))
print(summary(lm2))
print(coef(lm2))

plot(mpdf$mpc1~mpdf$pa,col=plotcols[1],ylab="MPC",xlab="patch area (km^2)",main='untransformed')
#lines(mpdf$pa,mpdf$pa)                  # one to one line is really far away..
points(mpdf$mpc2~mpdf$pa,col=plotcols[2])
points(mpdf$mpc3~mpdf$pa,col=plotcols[3])
# add model predictions
lines(mpdf$pa, exp(predict(lm2, newdata=list(x=log(mpdf$pa)))),col="dodgerblue",lwd=2)
# dispfuns not relevant in case of 1 patch...
legend("topleft",legend=c("negex","linear","log-sech"),
       title="dispfun",fill=plotcols,bty='n')

###
## show the change in log space..
plot(log(mpdf$mpc1)~mpdf$pa,col=plotcols[1],ylab="log(MPC)",xlab="patch area (km^2)",main='log-level')
#lines(mpdf$pa,mpdf$pa)                  # one to one line is really far away..
points(log(mpdf$mpc2)~mpdf$pa,col=plotcols[2])
points(log(mpdf$mpc3)~mpdf$pa,col=plotcols[3])
# add model predictions
lines(mpdf$pa, predict(lm2, newdata=list(x=log(mpdf$pa))),col="dodgerblue",lwd=2)

###
## show the change in log-log space..
plot(log(mpdf$mpc1)~log(mpdf$pa),col=plotcols[1],ylab="log(MPC)",xlab="log[patch area (km^2)]",main='log-log')
points(log(mpdf$mpc2)~log(mpdf$pa),col=plotcols[2])
points(log(mpdf$mpc3)~log(mpdf$pa),col=plotcols[3])
# add model predictions
lines(log(mpdf$pa), predict(lm2, newdata=list(x=log(mpdf$pa))),col="dodgerblue",lwd=2)
lines(log(mpdf$pa),log(mpdf$pa)) # one to one line..
legend("topleft",legend=c(paste("int:",round(coef(lm2),3)[1],"coef:",round(coef(lm2),3)[2],sep=" "),"y=x"),
       title="log-log model",fill=c("dodgerblue","black"),bty='n')

sink()
dev.off()


##########################################################
## Dependence on distance analysis: change in metapopulation capacity across a range of distances - with 2 equal habitat patch areas
#
alpha = 317       # Assuming constant alpha = 317m
pa    = c(50.5,50.5)  # Assuming a constant patch area of two patches with each 50 km2 (i.e. together = 101km2 area)
dvec  = c(seq(1,100*alpha,by=10)) # dist = 0 makes no sense, is covered in mdf with 1 patch of 100km2
compareval1 = mpdf[mpdf$pa==101,"mpc1"] # same area as the two patches
# compareval2 = mpdf[mpdf$pa==51,]  # ~ area of 1 patch: actually make it more exact:
####
## One maximum area
dfamax   = data.frame(patchID=1,aream2=50.5)
pamax    = dfamax$aream2
mdistmax = matrix(ncol=1,nrow=1,byrow=TRUE,c(0.0000))
## Calculate MPC
compareval2 = MPC_fun(pa=pamax,mdist=mdistmax,dispfop="negex",alpha=317)$mpc


#####
## make empty dataframe for adding calculated MPC values
mddf = data.frame(dist=numeric(0),mpc1=numeric(0),mpc2=numeric(0),mpc3=numeric(0))

## loop through avec for MPC calculation
for(dist in dvec) {
  print(dist)
  ## Pairwise patch distances - parameter "mdist"
  mdist  = matrix(ncol=2,nrow=2,byrow=TRUE,c(0.00,dist,
                                             dist,0.00))

  # parameters are loaded above
  mpc1=MPC_fun(pa=pa,mdist=mdist,alpha=alpha,dispfop="negex")$mpc
  mpc2=MPC_fun(pa=pa,mdist=mdist,alpha=alpha,dispfop="linear")$mpc
  mpc3=MPC_fun(pa=pa,mdist=mdist,alpha=alpha,dispfop="log-sech")$mpc
  mddf0=data.frame(dist=dist,mpc1=mpc1,mpc2=mpc2,mpc3=mpc3)
  mddf=rbind(mddf,mddf0)
}

## write to file
write.csv(mddf,paste(dir,"demo/data/MPC_value_range_dist.csv",sep="/"),row.names = FALSE)

### make a plot
var  = 'dist'    # for easier plotting below..
data = mddf

pdf(paste(dir,"demo/data/MPC_value_range_dist.pdf",sep="/"),width=6,height=4)
sink(paste(dir,"demo/data/MPC_value_range_dist.txt",sep="/"))
## preparations
par(mfrow=c(1,3))
plotcols = c("black","olivedrab3","coral")
print("---summary of MPC in one patch---")
print(summary(data))
# fit a linear model to estimate coefficients! since all mpcs (from different dispun's) are the same because only 1 patch here, do it only with 1
print("---relationship not linear in log-log space---")
# lm2  = lm(log(data$mpc1) ~ log(data[[var]]))
# print(summary(lm2))
# print(coef(lm2))

plot(data$mpc1~data[[var]],col=plotcols[1],ylab="MPC",xlab="distance (m)",main='untransformed')
#lines(mpdf$pa,mpdf$pa)                  # one to one line is really far away..
points(data$mpc2~data[[var]],col=plotcols[2])
points(data$mpc3~data[[var]],col=plotcols[3])
# # add model predictions
# lines(data[[var]], exp(predict(lm2, newdata=list(x=log(data[[var]])))),col="dodgerblue",lwd=2)
# dispfuns not relevant in case of 1 patch...
##abline with comparison value:
##abline(h=compareval$mpc1,col="red")
legend("topleft",legend=c("negex","linear","log-sech"),
       title="dispfun",fill=plotcols,bty='n')

###
## show the change in log space..
plot(log(data$mpc1)~data[[var]],col=plotcols[1],ylab="log(MPC)",xlab="distance (m)",main='log-level')
points(log(data$mpc2)~data[[var]],col=plotcols[2])
points(log(data$mpc3)~data[[var]],col=plotcols[3])
# # add model predictions
# lines(data[[var]], predict(lm2, newdata=list(x=log(data[[var]]))),col="dodgerblue",lwd=2)
legend("topleft",legend=c(paste("101=", round(compareval1,1),sep=" "),
                          paste("101 log=", round(log(compareval1),1),sep=" "),
                          paste("50.5=", round(compareval2,1),sep=" "),
                          paste("50.5 log=", round(log(compareval2),1),sep=" ")
                          ),
       title="MPC for d=0:",fill=c("red"),bty='n',cex=0.8)

###
## show the change in log-log space..
plot(log(data$mpc1)~log(data[[var]]),col=plotcols[1],ylab="log(MPC)",xlab="log[distance (m)]",main='log-log')
points(log(data$mpc2)~log(data[[var]]),col=plotcols[2])
points(log(data$mpc3)~log(data[[var]]),col=plotcols[3])
# # add model predictions
# lines(log(data[[var]]), predict(lm2, newdata=list(x=log(data[[var]]))),col="dodgerblue",lwd=2)
# lines(log(data[[var]]),log(data[[var]])) # one to one line..
# legend("topleft",legend=c(paste("int:",round(coef(lm2),3)[1],"coef:",round(coef(lm2),3)[2],sep=" "),"y=x"),
#        title="log-log model",fill=c("dodgerblue","black"),bty='n')

sink()
dev.off()


########################################################
### End
########################################################





