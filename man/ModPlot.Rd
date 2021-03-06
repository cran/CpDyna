% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModPlot.R
\name{ModPlot}
\alias{ModPlot}
\title{Function to plot the results.}
\usage{
ModPlot(Gnu, res, type = "graphs")
}
\arguments{
\item{Gnu}{A (3 x Nb.events) matrix. The first row contains the ordered interaction times, the second row the source nodes, the third row the destination nodes. Currently, only undirected interactions are taken into account. Thus, for each interaction time, the source node is stricly smaller than the destination node. No self loops allowed.}

\item{res}{ModFit function's output.}

\item{type}{String that indicates the type of graph to be plotted. It must be either "graphs" or "adjacency".}
}
\description{
Function to plot the results.
}
\examples{

    data("Gnu")
    Th <- 10  # final time T
    step <- 0.1
    myptn <- seq(from = step, to = Th, by = step)
    res <- ModFit(Gnu = Gnu, kmin = 3, kmax = 3, MinimalPartition = FALSE, custom_partition = myptn)
  	ModPlot(Gnu, res, type = "adjacency")
  
}
