% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_prep_data.R
\name{filter_outliers}
\alias{filter_outliers}
\title{Filter outliers}
\usage{
filter_outliers(
  tissue = NULL,
  sex = NULL,
  outliers = MotrpacRatTraining6moData::OUTLIERS$viallabel
)
}
\arguments{
\item{tissue}{optional character, tissue abbreviation, one of \link[MotrpacRatTraining6moData:TISSUE_ABBREV]{MotrpacRatTraining6moData::TISSUE_ABBREV}}

\item{sex}{optional character, one of 'male' or 'female'}

\item{outliers}{vector of vial labels to consider as outliers.
Defaults to vial labels in \link[MotrpacRatTraining6moData:OUTLIERS]{MotrpacRatTraining6moData::OUTLIERS}.}
}
\value{
character vector, subset of \code{outliers} that correspond to the
specified tissue and sex
}
\description{
Filter a list of outliers to those belonging to the specified dataset.
Used to specify sex-specific outliers within differential analysis functions.
}
\examples{
curr_outliers = filter_outliers(tissue="HIPPOC")
curr_outliers = filter_outliers(tissue="HIPPOC", sex="male")
}
