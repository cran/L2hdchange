#' U.S. COVID-19 Data
#'
#' Daily number of COVID-19 cases for 58 areas in the United States (including 50 states,
#' Washington D.C., 5 territories and 2 cruise ships)
#' for 812 days from 22 Jan 2020 to 12 April 2022.
#'
#' @format ## `covid_data`
#' A data matrix with p = 58 rows and n = 812 columns.
#'
#' @source U.S. CDC \url{https://covid.cdc.gov/covid-data-tracker/#maps_new-admissions-rate-county}
"covid_data"

#' U.S. COVID-19 Data Neighbourhood Information
#'
#' @format ## `covid_nbd_info`
#' A list containing five arrays indicating the constituents of five U.S. regions:
#'
#' \strong{Northeast}: Connecticut, Maine, Massachusetts, New Hampshire, Rhode Island, Vermont,
#' New Jersey, New York, and Pennsylvania.
#'
#' \strong{Midwest}: Illinois, Indiana, Michigan, Ohio, Wisconsin, Iowa, Kansas, Minnesota, Missouri,
#' Nebraska, North Dakota, and South Dakota.
#'
#' \strong{South}: Delaware, Florida, Georgia, Maryland, North Carolina, South Carolina, Virginia,
#' District of Columbia, West Virginia, Alabama, Kentucky, Mississippi, Tennessee, Arkansas, Louisiana,
#' Oklahoma, and Texas.
#'
#' \strong{West}: Arizona, Colorado, Idaho, Montana, Nevada, New Mexico, Utah, Wyoming, Alaska,
#' California, Hawaii, Oregon, and Washington.
#'
#' \strong{Others}: American Samoa, Diamond Princess, Grand Princess, Guam, Northern Mariana
#' Islands, Puerto Rico, and Virgin Islands.
#'
#' @source U.S. Census Bureau, W. (2000). \emph{List of regions of the United States}.
"covid_nbd_info"
