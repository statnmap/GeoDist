# GeoDist

**Work in progress**

## Constrained distance calculation and associated geotools

This R-package allows the calculation of distances that are contrained by frontiers, islands, mountains,... Distances are calculated from a raster grid of the area.  
_From A to B_
       
       A     XXX
       |   XXXXXX
       | XXXXXXXX
       |XXXXXXXX
       \ XXXXXX
        \_XX_ _ B

These distances can then be implemented in classical geotools like kriging with a modified version of `geoR` functions that accept custom calculated distances (`variog.dist`, `likfit.dist`, `krige.conv.dist`).

Currently distances are calculated using library `igraph`, but library `Rvcg` would be much more efficient (I keep that in mind for one day...). 

