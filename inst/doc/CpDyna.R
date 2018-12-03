## ------------------------------------------------------------------------
set.seed(1)

## ---- echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE-----------------
library(CpDyna)

## ------------------------------------------------------------------------
data("Gnu")

## ------------------------------------------------------------------------
tail(Gnu[1,])

## ------------------------------------------------------------------------
custom_ptn <- c(1:100)

## ---- echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE-----------------
res <- ModFit(Gnu, 
              4, 
              4, 
              MinimalPartition = FALSE, 
              N_initializations = 20, 
              custom_partition = custom_ptn)

## ---- echo=TRUE----------------------------------------------------------
res$est_z
res$est_cp

## ---- eval=FALSE, message=FALSE, warning=FALSE---------------------------
#  par(mfrow = c(1,4))
#  ModPlot(Gnu, res, type = "adjacency")

## ----echo=FALSE, out.width = "90%", out.extra='style="display: block;margin: auto;"', fig.cap=""----
library(knitr)    # For knitting document and include_graphics function

include_graphics("figures/AdjSnaps.png")

## ------------------------------------------------------------------------
data(wtab)
data(stations_df)

## ---- echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE----------------
#  
#  par(mfrow=c(1,1), bty="n", cex.axis=1.5, cex.lab=1.5)
#  hist(wtab$Start.Date, breaks = 60, xaxt="n", yaxt="n", main = "", xlab = "time (hours)")
#  last <- 86340      # roughly 24h
#  time.seq<-seq(from=1, to=last, by=3600)
#  axis(1, at=time.seq, lab=c(0:23), lwd=.2)
#  axis(2, lwd=.2)
#  

## ----echo=FALSE, out.width = "75%", out.extra='style="display: block;margin: auto;"', fig.cap=""----
library(knitr)    # For knitting document and include_graphics function

include_graphics("figures/historiginal.png")

## ---- echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE----------------
#  Gnu <- as.matrix(wtab)
#  Gnu <- t(Gnu)

## ---- echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE----------------
#  N <- max(Gnu[2,], Gnu[3,])   # number of nodes/stations
#  
#  step <- 900      # in seconds, corresponds to 15 minutes
#  custom_ptn <- seq(from = step, to = max(Gnu[1,]), by = step)  # user defind partition
#  
#  res <- ModFit(Gnu, 4, 4, eps=10^(-1),MinimalPartition = FALSE, custom_partition = custom_ptn, N_initializations = 1,  verbose = TRUE)
#  

## ---- echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE----------------
#  cp <- res$est_cp
#  abline(v = cp, col="red", lwd=1.5)

## ----echo=FALSE, out.width = "75%", out.extra='style="display: block;margin: auto;"', fig.cap=""----
library(knitr)    # For knitting document and include_graphics function

include_graphics("figures/histCpDyna.png")

## ---- echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE----------------
#  
#  library(mapview)
#  library(sp)
#  
#  
#  new_df <- stations_df[, c(1,2,4,5)]
#  
#  #converting columns from factor to proper formats
#  sq<-c(1,3,4)
#  for(i in 1:length(sq)) new_df[,sq[i]]<-as.numeric(levels(new_df[,sq[i]]))[new_df[,sq[i]]]
#  new_df[,2]<-as.character(levels(new_df[,2]))[new_df[,2]]
#  
#  WhoIsWho <- seq(1:N)
#  
#  match_pos <- res$est_z[match(new_df$id, WhoIsWho)]  # for each station (id) in new_df, i look for its position in WhoIsWho and take outz at this position
#  new_df<-cbind(new_df, match_pos)
#  pos_na <- which(is.na(new_df$match_pos))
#  new_df <- new_df[-pos_na, ]
#  
#  
#  
#  tav <-c ("RoyalBlue3","red", "green3","gold")
#  sub_df <- new_df[,c(3,4)]
#  coordinates(sub_df) <- ~long+lat
#  proj4string(sub_df) <- CRS("+init=epsg:4326")
#  mapview(sub_df, color=tav[new_df$match_pos], cex= 0.5, alpha = 0.8, lwd=6)
#  

## ----echo=FALSE, out.width = "75%", out.extra='style="display: block;margin: auto;"', fig.cap=""----
library(knitr)    # For knitting document and include_graphics function

include_graphics("figures/map.png")

