#' Data from Homicidios (murder and murder rate)
#'
#' This data contains the total number of homicides and the intentional homicide
#' rate in Colombia for the years 1960 to 2018.
#' The data were drawn from the National Planning Department (DNP) and the
#' criminal statistics of the National Police.
#'
#' @format A data frame with 50 observations and 4 variables:
#' \describe{
#'   \item{Year}{Year in which the homicide was reported.}
#'   \item{Murder}{Total number of homicides registered for a given year.}
#'   \item{Population}{Total population in Colombia for a given year.}
#'   \item{Rate}{Homicide rate per year per 100.000 inhabitants.}
#'   ...
#' }
"homicides"

#' Data from Incendios Forestales (wildfires)
#'
#' These data describe the number of large forest fires that occurred in Colombia,
#' from January 2002 to December 2016. The data was extracted from the Institute
#' of Hydrology, Meteorology and Environmental Studies (IDEAM), and the periodicity
#' of the data is monthly.
#'
#' @format A data frame with 180 rows and 2 variables:
#' \describe{
#'   \item{Date}{Year and month in which the wildfire was reported.}
#'   \item{GIF}{Grandes Incendios Forestales. Large wildfires are defined as those fires that exceed 500 acres affected forest.}
#'   ...
#' }
"wildfires"
