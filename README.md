## Content
The MPC R-package essentially contains a wrapper function for the derivation of metapopulation capacity (MPC) from binary habitat distribution maps (shapefile, raster file or patch-area vector/patch-distance matrix). Additionally it enables the individual setting of species-specific parameters, such as dispersal capacity and distance-decay function. 
The implemented MPC formula is based on Hanski and Ovaskainen 2000, Hanski 1994, as well as the modifications suggested by Schnell et al. 2013.
The R-code is largely based on the R-codes of Huang et al. 2020 and Strimas-Mackey & Brodie 2018.
The MPC R-package was developed for the project "Linking multispecies connectivity modelling and ecosystem services in the context of landscape urbanization", a postdoc project collaboration of McGill University and the environmental research firm Habitat (https://www.habitat-nature.com/) that is co-supervised by Dr. Andrew Gonzalez and Dr. Brian Leung.

## Full MPC report
A detailed report on the metapopulation capacity indicator and related indices can be accessed here: https://oehrij.shinyapps.io/MPC_report/

## References
1) Hanski, Ilkka, and Otso Ovaskainen. 2000. The Metapopulation Capacity of a Fragmented Landscape. Nature 404 (6779): 755-58.  https://doi.org/10.1038/35008063
2) Schnell, Jessica K., Grant M. Harris, Stuart L. Pimm, and Gareth J. Russell. 2013. Estimating Extinction Risk with Metapopulation Models of Large-Scale Fragmentation. Conservation Biology 27(3): 520-30  https://doi.org/10.1111/cobi.12047
3) Strimas-Mackey, M., & Brodie, J. F. (2018). Reserve design to optimize the long-term persistence of multiple species. Ecological Applications, 28(5), 1354-1361. https://doi.org/10.1002/eap.1739; https://github.com/mstrimas/metacapa  
4) Huang, R., Pimm, S. L., & Giri, C. (2020). Using metapopulation theory for practical conservation of mangrove endemic birds. Conservation Biology, 34(1), 266-275. https://doi.org/10.1111/cobi.13364
5) Hanski, I. (1994). A practical model of metapopulation dynamics. Journal of animal ecology, 151-162.https://doi.org/10.2307/5591
