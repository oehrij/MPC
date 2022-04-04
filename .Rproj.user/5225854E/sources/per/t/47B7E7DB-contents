#' @name MPC
#' @title Metapopulation capacity (MPC) core function
#' @description Calculates MPC for a given set of habitat patches. The function relies on the 'sf' and 'stars'packages
#' @param pa A numeric vector of areas of habitat patches 1-n 
#' @param mdist A square matrix of pairwise distances between habitat patches 1-n 
#' @param alpha A species-specific dispersal distance (e.g. average gap-crossing distance), default=317m
#' @param dispfun A species-specific dispersal survival function, currently 3 options can be chosen by parameter dispfop: "negex", "linear" and "log-sech"
#' @param ex A stochasticity parameter: if x >1, it becomes very unlikely for populations to go extinct after a critical patch size has been reached. If x <1, there exists no critical patch size and populations can even go extinct if they are large (Hanski 1994)
#' @param self Logical: should self-colonization of patches be modeled (default, TRUE) or not (FALSE)
#' @param evec Logical: should the dominant eigenvalue associated eigenvector be returned (default, TRUE) or not (FALSE)
#' @return An 'MPC' object: a list containing a numeric MPC value (mpc), the MPC-associated squared dominant eigenvector indicating patch importance (pimport) and the patch IDs in case they are indicated (pid)
#' @export
MPC         = function(pa,mdist,alpha,dispfun,ex,self,evec) {
  
  #matrix of dispersal survival probabilities
  dispr       = as.matrix(dispfun(mdist))                                                # Haung et al. 2020: Convert distances measurements to probability of dispersal
  
  #add self-colonization option (1=on, 0=off;Strimas-Mackey & Brodie 2018)
  diag(dispr) = ifelse(self, 1, 0)                                                       # Strimas-Mackey & Brodie 2018: the self is only if you explicitly want to exclude the self; because all current distance functions yield 1 for dist =0...
  
  #create 2 matrices for patch areas (Huang et al. 2020)         
  fragArea = matrix(pa, ncol=1)                                                           # Huang et al. 2020: actually take areas in km^2
  j.area   = matrix(fragArea,nrow=length(fragArea),ncol = length(fragArea),byrow = TRUE)
  i.area   = matrix(fragArea,nrow=length(fragArea),ncol = length(fragArea),byrow = FALSE)
  
  #create the colonization matrix (Huang et al. 2020)
  col.matrix = dispr*j.area 
  
  #create the metapopulation matrix (Huang et al. 2020)
  M = col.matrix*i.area^ex                                                               # Cf. Eqn 1 in Huang et al. 2020
  
  ## calculate metapopulation capacity (Strimas-Mackey & Brodie 2018)
  MPC0   = eigen(M, only.values = !evec)                                                 # Strimas-Mackey & Brodie 2018: create associated eigenvector if evec=TRUE
  if (evec) {                                                                            # Strimas-Mackey & Brodie 2018: note: Since the 'eigen()' function sorts values in descending order, the dominant eigenvalue appears first (check also the Perron–Frobenius Theorem: all real positive square matrices will have a positive dominant eigenvalue)
    MPC <- list(mpc = Re(MPC0$values[1]), pimport = Re(MPC0$vectors[, 1]) ^ 2, 
                pid = colnames(mdist))                                                   # Strimas-Mackey & Brodie 2018: Re(x) only return the real part of the number
    } else {                                                                             # Strimas-Mackey & Brodie 2018: exponent 2 to highlight the absolute weight of the eigenvector
    MPC <- list(mpc = Re(MPC0$values[1]), pimport = NA, pid=NA)                            
    }
  return(MPC)
}

#' @name MPC_fun
#' @title Metapopulation capacity (MPC) wrapper function
#' @description Calculates MPC for a given set of habitat patches. The function relies on the 'sf' and 'stars'packages
#' @param x A raster object where habitat patches are set to 1 and matrix area is set to NA OR a shapefile object containing habitat patches
#' @param pa A numeric vector of areas of habitat patches 1-n 
#' @param mdist A square matrix of pairwise distances between habitat patches 1-n 
#' @param alpha A species-specific dispersal distance (e.g. average gap-crossing distance), default=317m
#' @param dispfop An option to choose the dispersal survival function (dispfun parameter), currently 3 options implemented: "negex", "linear" and "log-sech"
#' @param ex A stochasticity parameter: if x >1, it becomes very unlikely for populations to go extinct after a critical patch size has been reached. If x <1, there exists no critical patch size and populations can even go extinct if they are large (Hanski 1994)
#' @param self Logical: should self-colonization of patches be modeled (default, TRUE) or not (FALSE)
#' @param evec Logical: should the dominant eigenvalue associated eigenvector be returned (default, TRUE) or not (FALSE)
#' @param distfun A function defining how interpatch distances are calculated in case the input is a raster or shapefile (default: euclidean, nearest, edge-to-edge: sf::st_distance)
#' @param areafun A function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @return the 'MPC' object created by the function 'MPC' 
#' @export
MPC_fun     = function(x=NULL,pa=NULL,mdist=NULL,alpha=317,dispfop="log-sech",ex=0.5,self=TRUE,evec=TRUE,
                       distfun=sf::st_distance,areafun=sf::st_area) {
  ## prepare data: make objects pa (numeric vector of patch areas) and mdist (square matrix of interpatch distances)
  if(class(x)[1]=="RasterLayer") {
    x <- stars::st_as_stars(x) %>% sf::st_as_sf(merge = TRUE) # if it is a raster, make a shapefile out of it
    x <- sf::st_make_valid(x,reason=TRUE)                        # make sure polygons are correctly extracted
  } # if raster
  if(class(x)[1]=="sf") {                                        # if it is a shapefile, extract the patch areas (numeric vector) & distance matrix
    ID    = x[[1]]
    pa    = as.numeric(areafun(x))
    mdist = distfun(x)
    mdist = matrix(mdist, dim(mdist)[1], dim(mdist)[2])
    rownames(mdist) <- ID
    colnames(mdist) <- ID
  } # if sf
  
  ## if x is not specified but mdist and pa are, we can directly go ahead:
  if(class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]) {
    
    ## define the dispersal function
    if(dispfop=="negex")   {dispfun= function(d) { exp(-(1/alpha) * abs(d)) }}                                # alpha=usually average dispersal distance
    else if(dispfop=="linear") {dispfun= function(d) { d <- abs(d); ifelse(d < alpha, 1 - d/alpha, 0) }}      # alpha=usually maximum dispersal distance
    else if(dispfop=="log-sech") {dispfun= function(d) { beta=1.77; (2*atan((alpha/d)^(1/(1/(beta-1)))))/pi}} # beta =thickness of tail in fat-tailed distribution # alpha=usually average dispersal distance
    
    ## actually calculate the metapopulation capacity
    MPC_obj = MPC(pa=pa,mdist=mdist,alpha=alpha,dispfun=dispfun,ex=ex,self=self,evec=evec)
    
  } else {
    print("--- error --- pa or mdist do not conform to one of the following ---")
    print('class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]')
    MPC_obj = "error"
  } # else
  
  return(MPC_obj)
}


##########################################################
##### project:       multispecies connectivity modelling 
##### author:        Jacqueline Oehri (JO)
##### date:          29.03.2022
##### last modified: 29.03.2022
##### comments:      -references: 
#####                 1) Hanski, Ilkka, and Otso Ovaskainen. 2000. “The Metapopulation Capacity of a Fragmented Landscape.” Nature 404 (6779): 755–58. doi:10.1038/35008063.
#####                 2) Schnell, Jessica K., Grant M. Harris, Stuart L. Pimm, and Gareth J. Russell. 2013. “Estimating Extinction Risk with Metapopulation Models of Large-Scale Fragmentation.” Conservation Biology 27 (3): 520–30. doi:10.1111/cobi.12047.
#####                 3) Strimas‐Mackey, M., & Brodie, J. F. (2018). Reserve design to optimize the long‐term persistence of multiple species. Ecological Applications, 28(5), 1354-1361. https://doi.org/10.1002/eap.1739; https://github.com/mstrimas/metacapa;  
#####                 4) Huang, R., Pimm, S. L., & Giri, C. (2020). Using metapopulation theory for practical conservation of mangrove endemic birds. Conservation Biology, 34(1), 266-275. https://doi.org/10.1111/cobi.13364
#####                 5) Hanski, I. (1994). A practical model of metapopulation dynamics. Journal of animal ecology, 151-162.
#########################################################


