% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathway_enrichment_networks.R
\name{cleanup}
\alias{cleanup}
\title{Clean up pathway names}
\usage{
cleanup(x)
}
\arguments{
\item{x}{character, pathway names}
}
\value{
vector of candidate words to determine overlap
}
\description{
Remove punctuation, numbers, and common words to determine qualitative
overlap between pathways enriched by metabolites versus other pathways.
Function used internally in \code{\link[=enrichment_network_vis]{enrichment_network_vis()}}.
}
\seealso{
\code{\link[=enrichment_network_vis]{enrichment_network_vis()}}
}
\keyword{internal}
