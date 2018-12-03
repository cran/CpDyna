EqualTo<-function(v1, v2) return((adjustedRandIndex(v1,v2)==1)) 
CleanLines<-function(M){
  i1<-1
  while (i1<nrow(M)){
    i2<-i1+1
    ToClean<-c()
    while(i2<=nrow(M)){
      if (EqualTo(M[i1,], M[i2,])) {
        ToClean<-c(ToClean, i2)
        i2<-i2+1
      }
      else i2<-i2+1
    }
    if (!is.null(ToClean)) M<-M[-ToClean,]
    i1<-i1+1
    if (!is.matrix(M)) break
  }
  return(M)
}