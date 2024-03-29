% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_graphical_clustering.R
\name{extract_main_clusters}
\alias{extract_main_clusters}
\title{Extract main graphical clusters}
\usage{
extract_main_clusters()
}
\value{
a data frame with 5 columns and one row per combination of
feature ID and cluster:
\describe{
\item{\code{feature}}{character, unique feature identifier in the format '\link[MotrpacRatTraining6moData:ASSAY_ABBREV]{MotrpacRatTraining6moData::ASSAY_ABBREV};\link[MotrpacRatTraining6moData:TISSUE_ABBREV]{MotrpacRatTraining6moData::TISSUE_ABBREV};feature_ID' only for training-regulated features at 5\% IHW FDR. For redundant differential features, 'feature_ID' is prepended with the specific platform to make unique identifiers. See \link[MotrpacRatTraining6moData:REPEATED_FEATURES]{MotrpacRatTraining6moData::REPEATED_FEATURES} for details.}
\item{\code{cluster}}{character, cluster label}
\item{\code{ome}}{character, assay abbreviation, one of \link[MotrpacRatTraining6moData:ASSAY_ABBREV]{MotrpacRatTraining6moData::ASSAY_ABBREV}}
\item{\code{tissue}}{character, tissue abbreviation, one of \link[MotrpacRatTraining6moData:TISSUE_ABBREV]{MotrpacRatTraining6moData::TISSUE_ABBREV}}
\item{\code{feature_ID}}{character, MoTrPAC feature identifier}
}
}
\description{
Return a data frame with features from the 2 largest nodes, 2 largest edges, 10 largest non-null paths,
and all 8-week nodes from the graphical representation of training-regulated
features in each tissue. This code replicates the graphical clusters for which
pathway enrichment was performed for the landscape manuscript.
}
\details{
Notes about cluster labels:
\itemize{
\item All clusters are prefixed with the tissue abbreviation and a colon, e.g. "SKM-GN:"
\item Nodes are defined by the time point and state in each sex, where state is 1 for up,
0 for null, and -1 for down. For example, "1w_F-1_M-1" is a node that characterizes molecules
at the "1w" time point that are down-regulated in females ("F-1") and down-regulated in males ("M-1").
These three pieces of information (time point, female state, male state) are separated by underscores ("_")
\item Edges contain "---" and connect a pair of nodes
\item Paths contain "->" and connect four nodes
}
}
\examples{
cluster_df = extract_main_clusters()

}
