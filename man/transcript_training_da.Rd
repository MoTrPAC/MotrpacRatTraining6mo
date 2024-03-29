% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transcript_differential_analysis.R
\name{transcript_training_da}
\alias{transcript_training_da}
\title{RNA-seq training differential analysis}
\usage{
transcript_training_da(
  tissue,
  covariates = c("pct_globin", "RIN", "pct_umi_dup", "median_5_3_bias"),
  outliers = na.omit(MotrpacRatTraining6moData::OUTLIERS$viallabel),
  rdata_outfile = NULL,
  overwrite = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{tissue}{character, tissue abbreviation, one of \link[MotrpacRatTraining6moData:TISSUE_ABBREV]{MotrpacRatTraining6moData::TISSUE_ABBREV}}

\item{covariates}{character vector of covariates that correspond to column names of \link[MotrpacRatTraining6moData:TRNSCRPT_META]{MotrpacRatTraining6moData::TRNSCRPT_META}.
Defaults to covariates that were used for the manuscript.}

\item{outliers}{vector of viallabels to exclude during differential analysis. Defaults
to \code{[MotrpacRatTraining6moData::OUTLIERS]$viallabel}}

\item{rdata_outfile}{NULL or path in which to save DESeq2 objects in an RData file}

\item{overwrite}{boolean, whether to overwrite the file if \code{rdata_outfile} exists}

\item{verbose}{boolean, whether to print messages}
}
\value{
a data frame with one row per gene:
\describe{
\item{\code{feature_ID}}{character, MoTrPAC feature identifier}
\item{\code{assay}}{character, assay abbreviation, one of \link[MotrpacRatTraining6moData:ASSAY_ABBREV]{MotrpacRatTraining6moData::ASSAY_ABBREV}}
\item{\code{assay_code}}{character, assay code used in data release. See \link[MotrpacBicQC:assay_codes]{MotrpacBicQC::assay_codes}.}
\item{\code{tissue}}{character, tissue abbreviation, one of \link[MotrpacRatTraining6moData:TISSUE_ABBREV]{MotrpacRatTraining6moData::TISSUE_ABBREV}}
\item{\code{tissue_code}}{character, tissue code used in data release. See \link[MotrpacBicQC:bic_animal_tissue_code]{MotrpacBicQC::bic_animal_tissue_code}.}
\item{\code{removed_samples_male}}{character, comma-separated list of male outliers (vial labels) removed from differential analysis}
\item{\code{removed_samples_female}}{character, comma-separated list of female outliers (vial labels) removed from differential analysis}
\item{\code{lrt_male}}{double, likelihood ratio test statistic for males}
\item{\code{lrt_female}}{double, likelihood ratio test statistic for females}
\item{\code{p_value_male}}{double, nominal LRT p-value for males}
\item{\code{p_value_female}}{double, nominal LRT p-value for females}
\item{\code{full_model_male}}{character, full model used in LRT for males}
\item{\code{full_model_female}}{character, full model used in LRT for females}
\item{\code{reduced_model_male}}{character, reduced model used in LRT for males}
\item{\code{reduced_model_female}}{character, reduced model used in LRT for females}
\item{\code{p_value}}{double, combined male and female nominal p-value using the sum of logs}
}
}
\description{
Use DESeq2 to perform a likelihood ratio test to test the effect of training
across time points. Analysis is performed separately for males and females.
}
\examples{
\dontrun{
# Perform differential analysis for expressed genes in brown adipose tissue 
# with default parameters, i.e., outliers and covariates used for the manuscript
da = transcript_training_da("BAT")

# Same as above but save the [DESeq2::DESeq2()] DESeqResults objects in an RData file 
da = transcript_training_da("BAT", 
                            rdata_outfile = "~/test/BAT_RNA_training-da.RData", 
                            overwrite = TRUE)
}
}
