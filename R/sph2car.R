#' Computes cartesian coordinates from long,lat geographical coordinates
#' Use of geoid with radius and flattening to calculate correct x,y,z coordinates
#'
#' @param long longitude values, can also contain a matrix of (long, lat), or
#' (long, lat and and radius) in that order.
#' @param lat	latitude values.
#' @param radius major (equatorial) radius of the ellipsoid. The default value is for WGS84
#' @param f numeric. Ellipsoid flattening. The default value is for WGS84
#' @param deg Specifies if input is in degrees (default) or radians.
#' @export
#'
#' @details
#' Use geosphere::refEllipsoids() for other radius and f possibilites.
#'
sph2car <- function(long, lat, radius = 6378137, f = 1/298.257223563, deg = TRUE)
{
  if (is.matrix(long) || is.data.frame(long)) {
    if (ncol(long) == 1) {
      long = long[, 1]
    }
    else if (ncol(long) == 2) {
      lat = long[, 2]
      long = long[, 1]
    }
    else if (ncol(long) == 3) {
      radius = long[, 3]
      lat = long[, 2]
      long = long[, 1]
    }
  }
  if (missing(long) | missing(lat)) {
    stop("Missing full spherical 3D input data.")
  }
  if (deg) {
    long = long * pi/180
    lat = lat * pi/180
  }

  #https://de.mathworks.com/help/aeroblks/llatoflatearth.html?requestedDomain=www.mathworks.com
  #R.n Prime vertical
  Rn <- radius / sqrt(1 - (2*f - f*f) * (sin(lat)^2))
  # Rm Radius meridian
  Rm <- Rn * (1 - (2*f - f*f) ) / (1 - (2*f - f*f) * (sin(lat)^2))

  cbind(x = Rm * cos(long) * cos(lat), y = Rm *
          sin(long) * cos(lat), z = Rm * sin(lat))
}
