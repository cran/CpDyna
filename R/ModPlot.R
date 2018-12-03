# function to plot the results
#' Function to plot the results.
#' 
#' @param Gnu A (3 x Nb.events) matrix. The first row contains the ordered interaction times, the second row the source nodes, the third row the destination nodes. Currently, only undirected interactions are taken into account. Thus, for each interaction time, the source node is stricly smaller than the destination node. No self loops allowed.
#' @param res ModFit function's output.
#' @param type String that indicates the type of graph to be plotted. It must be either "graphs" or "adjacency".
ModPlot <- function(Gnu, res, type = "graphs"){
  N <- length(res$est_z)
  if (type == "graphs"){
    pal <- brewer.pal(12, 'Paired')
    graphs <- GnuToX(Gnu, res$est_cp, N)
    D <- res$est_D
    for (d in 1:D){
      ig <- igraph::graph.adjacency(graphs[,,d], mode = "undirected", weighted = TRUE)
      
      # nodes
      igraph::V(ig)$color = res$est_z
      igraph::V(ig)$label = " "
      igraph::V(ig)$label.color = rgb(0,0,.2,1)
      igraph::V(ig)$label.cex = .9
      #igraph::V(ig)$size = degree(ig)/2.5 # ou V(g)$size = degree(g) etc, idem pour les edges avec edge betweenness, closeness etc.
      igraph::V(ig)$size = 2.5
      igraph::V(ig)$frame.color <- igraph::V(ig)$frame.color <-rgb(160/255,160/255,160/255,1)
      
      
      #edges
      igraph::E(ig)$color <- rgb(204/255,204/255,204/255,.5)
      igraph::E(ig)$width <- .8
      igraph::E(ig)$arrow.size <- .025
      
      # layout
      if (d == 1) dl <- layout.fruchterman.reingold(ig)
      
      # plotting
      plot(ig, layout=dl, main=paste(" Snapshot ", d))
    }
  }
  if (type == "adjacency"){
    op <- par(bty = "n")
    graphs <- GnuToX(Gnu, res$est_cp, N)
    N <- dim(graphs)[1]
    K <- res$est_K
    for (d in 1:res$est_D){
      I <- graphs[,,d]
      I[I != 0] <- 1
      sZ <- seq(1:N)
      Zclus.copy <- res$est_z
      names(Zclus.copy) <- sZ
      Zclus.copy <- sort(Zclus.copy)
      posZ <- as.numeric(names(Zclus.copy))
      J <- I[posZ, posZ]
      image(t(J)[,N:1], col=gray(255:0/255), useRaster = TRUE, 
            xaxt = 'n', xlim = c(0,1), ylim = c(0,1), yaxt = 'n', main = paste("Snapshot ", d))  
      # hlines
      lev_hl <- cumsum(rev(table(Zclus.copy)))/N
      lev_hl <- lev_hl[-K]
      clip(0,1,0,1)
      abline(h = lev_hl, col = "red", lwd = 2)
      # vlines
      lev_vl <- cumsum(table(Zclus.copy))/N
      lev_vl <- lev_vl[-K]
      clip(0,1,0,1)
      abline(v = lev_vl, col = "red", lwd = 2)
    }
  }
}
