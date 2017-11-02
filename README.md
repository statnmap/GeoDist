# GeoDist: Constrained distance calculation and associated geotools

**Work in progress**

This R-package allows the calculation of distances that are contrained by frontiers, islands, mountains,... Distances are calculated from a raster grid of the area.  
_From A to B_
       
       A     XXX
       |   XXXXXX
       | XXXXXXXX
       |XXXXXXXX
       \ XXXXXX
        \_XX_ _ B

These distances can then be implemented in classical geotools like inverse distance interpolation (`idw.dist`) or kriging with a modified version of `geoR` functions that accept custom calculated distances (`variog.dist`, `likfit.dist`, `krige.conv.dist`).

Currently distances are calculated using library `igraph`, but library `Rvcg` would be much more efficient (I keep that in mind for one day...). 

# Download and Install

To download the development version of the `GeoDist` package, type the following at the R command line:

```r
install.packages("devtools")  
devtools::install_github("statnmap/GeoDist")
```

Note that spatial libraries like `rgdal` and `sp` may require additional softwares to be installed on your computer if you work with Mac or Linux. Look how to install `proj4`, `geos` and `gdal` on your system.
