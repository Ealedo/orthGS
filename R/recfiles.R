rm(list = ls())
library(igraph)

## --- Selected
x <- orthG()
As_ <- x[[1]]
gs_ <- x[[2]]

x <- mapTrees(paste(system.file("extdata", package = "orthGS"), "/selected", sep  = ""))
As <- x[[3]]
gs <- graph_from_adjacency_matrix(As, "undirected")

identical_graphs(gs, gs_) # Same graph

## --- Selected

c <- paste(system.file("extdata", package = "orthGS"), "/conifers", sep  = "")
mc <- mapTrees(path2rec = c)
mcc <- mc[[3]]
gc <- igraph::graph_from_adjacency_matrix(mcc, "undirected")
plot(gc)
degree(gc)

gy <- paste(system.file("extdata", package = "orthGS"), "/gymnosperms", sep  = "")
mgy<- mapTrees(path2rec = gy)
mgyy <- mgy[[3]]
ggy <- igraph::graph_from_adjacency_matrix(mgyy, "undirected")
plot(ggy)
degree(ggy)
