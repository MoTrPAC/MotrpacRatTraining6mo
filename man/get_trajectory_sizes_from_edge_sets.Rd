% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_graphical_clustering.R
\name{get_trajectory_sizes_from_edge_sets}
\alias{get_trajectory_sizes_from_edge_sets}
\title{Get trajectory sizes}
\usage{
get_trajectory_sizes_from_edge_sets(
  edge_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$edge_sets,
  min_size = 10
)
}
\arguments{
\item{edge_sets}{A named list of string vectors. The name of an edge is \code{[node_id]---[node_id]}
edges with no analytes have a NULL set (a set of size zero, but are still represented),
node ids are \code{[timepoint]_F[x]_M[y]} where \code{x} and \code{y} represent the up/down state in each sex.
\code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$edge_sets} by default.}

\item{min_size}{An integer specifying the minimal path size to be considered.}
}
\value{
NULL if no paths of size of at least \code{min_size} were found, otherwise
return a data frame that represents all paths of size min_size or greater,
ranked from the largest path to the smallest one.
}
\description{
Auxiliary function for getting the largest paths in a graphical solution.
}
\details{
This is implemented using a dynamic programming approach where
we iteratively add the data of the next edge.

When analyzing an edge (x,y) with a set of analytes s, we extend all
trajectories that end with x using the new edge, but also all the
trajectories that start with y.

At the end, because we examine all edges and all extensions we are
guaranteed to have covered all full paths.

The min_size parameter is important, since we are interested in paths
and edges with at least this number of analytes, then we know that if a
current trajectory does not have at least min_size analytes, then since the
set of any extension can only be the same or smaller then we can ignore
such paths moving forward.

However, this parameter has to be considered with care, as specifying a number
that is too high will result in no paths in the output.

See bayesian_graphical_clustering for more details about the graphical analysis.
}
\examples{
\dontrun{
### Example 1: Simulate data with a single cluster
zcolnames = c(
  paste("female",c("1w","2w","4w","8w"),sep="_"),
  paste("male",c("1w","2w","4w","8w"),sep="_")
)
zscores = matrix(rnorm(80000),ncol=8,dimnames = list(1:10000,zcolnames))
# now add a cluster with a strong signal and rerun
zscores[1:500,1:4] = zscores[1:500,1:4] + 5

# run the clustering solution wrapper
clustering_sol = bayesian_graphical_clustering(zscores)

# extract the top full trajectories in the data
# these should be the clusters with at least 10 features
min_cluster_size=10
get_trajectory_sizes_from_edge_sets(clustering_sol$edge_sets, 
                                    min_size = min_cluster_size)

# extract the edges of the top two full trajectories
# this step "cleans" the edge sets by removing edges of trajectories with very few features
top2traj_edge_sets = filter_edge_sets_by_trajectories(
  clustering_sol$edge_sets, 
  topk = 2, 
  min_path_size = 10
)
# examine the new edge set sizes, excluded edges should have zero size
sapply(top2traj_edge_sets,length)
}

### Example 2: Use published data 
get_trajectory_sizes_from_edge_sets()
}
