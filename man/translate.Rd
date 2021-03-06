% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{translate}
\alias{translate}
\title{Translate common climate data abbreviations.}
\usage{
translate(key, to = "words")
}
\arguments{
\item{key}{A character vector.}

\item{to}{A character vector, must be either "words", "units", "delta",
"abbv", or "letter". These are not all available for every type of key.}
}
\value{
A character vector of translated results.
}
\description{
This is a simple translation dictionary that can be used to look up things
like names and units, for climate variables, statistics, and months.
}
