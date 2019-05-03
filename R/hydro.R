#' Extraterrestrial solar radiation, per Shuttleworth 1993
#'
#' @param psi latitude, in degrees
#' @param J Julian day
#' @return solar radiation in mm/day
ETSR <- function(psi, J){

      # Shuttleworth's test:
      # the following should output c(15.0, 15.1, 11.2)
      # round(ETSR(psi=c(30,0,-30), J=105), 1)

      # convert degrees latitude to radians
      psi <- psi * pi / 180

      # solar declination
      delta <- 0.4093 * sin((J*2*pi/365) - 1.405)

      # sunset hour angle
      omega <- acos(-tan(psi) * tan(delta))

      # relative earth-to-sun distance
      dr <- 1 + 0.033 * cos(J*2*pi/365)

      # extraterrestrial solar radiation (mm / day)
      S0 <- 15.392 * dr * (omega * sin(psi) * sin(delta) + cos(psi) * cos(delta) * sin(omega))

      return(S0)
}



#' Mean S0 for a month of the year
#'
#' @param month integer from 1 to 12
#' @param latitude decimal degrees
#' @return mean S0 for a month of the year
monthly_S0 <- function(month, latitude){

      # julian days for target month
      j <- rep(1:12, times=c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
      jdays <- which(j==month)

      # average ETSR across all days of the month
      mean(ETSR(latitude, jdays))
}


#' Hargreaves equation for evapotranspiration
#'
#' Per Hargreaves & Samani 1985. (Shuttleworth 1993 seems to incorrectly omit the exponent)
#'
#' @param S0 numeric
#' @param tmean degrees C
#' @param tmin degrees C
#' @param tmax degrees C
#' @return evapotranspiration (ETP) in mm/day
hargreaves <- function(S0, tmean, tmin, tmax){
      0.0023 * S0 * (tmax - tmin)^.5 * (tmean + 17.8)
}



#' Annual water balance variables for one site
#'
#' @param latitude degrees
#' @param ppt vector of length 12, mm
#' @param tmean vector of length 12, degrees C
#' @param tmax vector of length 12, degrees C
#' @param tmin vector of length 12, degrees C
#' @param annual logical, if FALSE will return monthly values
#' @return A named vector of hydrologic variables: PPT, PET, AET, CWD, RAR
water_balance <- function(latitude, ppt, tmean, tmax, tmin, annual=TRUE){

      if(is.na(tmean[1])) return(rep(NA, 5))

      # [note: Wang et al 2012 page 21 use a latitude correction;
      # this note is a placeholder for that calculation if we decide to use it]

      # get S0 values for all 12 months
      S0 <- sapply(1:12, monthly_S0, latitude=latitude)

      # monthly water balance variables
      pet <- hargreaves(S0, tmean, tmin, tmax) *
            c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      pet[tmean<0] <- 0 # per Wang et al 2012
      aet <- pmin(pet, ppt)
      cmd <- pmax(0, pet - ppt)
      rr <- pmax(0, ppt - pet)

      m <- cbind(PPT=ppt, PET=pet, AET=aet, CWD=cmd, RAR=rr)
      if(annual) return(colSums(m)) else(return(m))
      #return(c(PPT=sum(ppt), PET=sum(pet), AET=sum(aet), CWD=sum(cmd), RAR=sum(rr)))
}




#' Create a latitude raster
#'
#' @param template a raster layer
#' @param already_latlong leave as F unless rasters are already in lat-long projection
#' @return A raster layer with values representing degrees latitude
latitude <- function(template, already_latlong=FALSE){

      requireNamespace("raster")
      requireNamespace("sp")

      lat <- template * 1 # force into RAM
      coords <- coordinates(lat)

      if(!already_latlong){
            coords <- as.data.frame(coords)
            coordinates(coords) <- c("x", "y")
            crs(coords) <- crs(lat)
            coords <- spTransform(coords, crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
            coords <- coordinates(coords)
      }

      lat[] <- coords[,2]
      return(lat)
}


#' Derive annual hydrologic variables
#'
#' This is the main external function in the package, which uses the Hargreaves
#' equation to generate a set of 5 annual water balance variables from a 48
#' monthly temperature and precipitation variables. The function calculates
#' cumulative annual totals for precipitation (PPT), potential
#' evapotranspiration (PET), actual evapotranspiration (AET), climatic water
#' deficit (CWD), and reacharge and runoff (RAR), all in mm. Note that it does
#' not work inside the arctic and antarctic circles.
#'
#' @param rasters stack of 48 rasters in this order: ppt1-12, tmean1-12,
#'   tmax1-12, tmin1-12
#' @param temp_scalar multiplier to convert input temperatures to deg C
#' @param ppt_scalar multiplier to convert input precipitation to mm
#' @param ncores number of computing cores to use for parallel processing
#' @param already_latlong leave as F unless rasters are already in lat-long
#'   projection
#' @param filename output filename
#' @param m see clusterR
#' @param ... additional arguments to clusterR or writeRaster
#' @return A raster stack with 5 layers: PPT, PET, AET, CWD, RAR
#' @references Wang et al 2012 -- ClimateWNA -- Journal of Applied Meteorology
#'   and Climatology. Shuttleworth 1993 -- Evaporation (chapter 4 in Handbook of
#'   Hydrology). Hargreaves and Samani 1985 -- Reference Crop Evaporation from
#'   Ambient Air Temperature
hydro <- function(rasters, temp_scalar=1, ppt_scalar=1,
                  ncores=1, already_latlong=F, filename, m){

      requireNamespace("raster")

      # latitude raster
      lat <- latitude(rasters[[1]], already_latlong=already_latlong)
      tempfile <- paste0(dirname(filename), "/temporary.tif")
      lat <- writeRaster(lat, tempfile, overwrite=T)
      rasters <- stack(rasters, lat)

      # compute annual water balance variables
      w <- function(x, ...) water_balance(latitude=x[49],
                                          ppt=x[1:12],
                                          tmean=x[13:24],
                                          tmax=x[25:36],
                                          tmin=x[37:48])
      if(temp_scalar != 1 | ppt_scalar != 1) w <- function(x, ...) water_balance(latitude=x[49],
                                                                                 ppt=x[1:12] * ppt_scalar,
                                                                                 tmean=x[13:24] * temp_scalar,
                                                                                 tmax=x[25:36] * temp_scalar,
                                                                                 tmin=x[37:48] * temp_scalar)

      beginCluster(ncores, type="SOCK")
      wb <- clusterR(rasters, calc, args=list(fun=w),
                     export=c("ETSR", "hargreaves", "water_balance", "monthly_S0"),
                     filename=filename, m=m)
      endCluster()

      names(wb) <- c("PPT", "PET", "AET", "CWD", "RAR")
      file.remove(tempfile)
      return(wb)
}


