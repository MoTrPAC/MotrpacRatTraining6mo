% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/immunoassay_differential_analysis.R
\name{immuno_training_da}
\alias{immuno_training_da}
\title{Immunoassay training differential analysis}
\usage{
immuno_training_da()
}
\value{
a data frame with one row per feature:
\describe{
\item{\code{feature_ID}}{character, MoTrPAC feature identifier}
\item{\code{assay}}{character, assay abbreviation, one of \link[MotrpacRatTraining6moData:ASSAY_ABBREV]{MotrpacRatTraining6moData::ASSAY_ABBREV}}
\item{\code{assay_code}}{character, assay code used in data release. See \link[MotrpacBicQC:assay_codes]{MotrpacBicQC::assay_codes}.}
\item{\code{dataset}}{character, LUMINEX panel}
\item{\code{tissue}}{character, tissue abbreviation, one of \link[MotrpacRatTraining6moData:TISSUE_ABBREV]{MotrpacRatTraining6moData::TISSUE_ABBREV}}
\item{\code{tissue_code}}{character, tissue code used in data release. See \link[MotrpacBicQC:bic_animal_tissue_code]{MotrpacBicQC::bic_animal_tissue_code}.}
\item{\code{removed_samples_male}}{character, comma-separated list of male samples (vial labels) removed from differential analysis}
\item{\code{removed_samples_female}}{character, comma-separated list of female samples (vial labels) removed from differential analysis}
\item{\code{lrt_male}}{double, likelihood ratio test statistic for males}
\item{\code{lrt_female}}{double, likelihood ratio test statistic for females}
\item{\code{p_value_male}}{double, nominal LRT p-value for males}
\item{\code{p_value_female}}{double, nominal LRT p-value for females}
\item{\code{full_model}}{character, full model used in LRT for males and females}
\item{\code{reduced_model}}{character, reduced model used in LRT for males and females}
\item{\code{p_value}}{double, combined male and female nominal p-value using the sum of logs}
}
}
\description{
For each panel and tissue, perform a likelihood ratio test to test the effect of training
across time points. Analysis is performed separately for males and females.
}
\examples{
res = immuno_training_da()
}
