
#' Node clustering and multiple change points detection on a dynamic graph.
#' @param Gnu A (3 x Nb.events) matrix. The first row contains the ordered interaction times, the second row the source nodes, the third row the destination nodes. Currently, only undirected interactions are taken into account. Thus, for each interaction time, the source node is stricly smaller than the destination node. No self loops allowed.
#' @param kmin A positive integer. The minimum number of clusters the BIC criterion will look for.
#' @param kmax A positive integer. The maximum number of clusters the BIC criterion will look for.
#' @param eps A small double. When the difference between two consecutive values of the variational lower bound is smaller than eps, the VEM algorithm stops.
#' @param iter.vem An integer. The maximum number of iterations the VEM algorithm will perform.
#' @param N_initializations An Integer. The number of different k-means initialisations to adopt.
#' @param MinimalPartition A Boolean. If TRUE (the default), the algorithm will look for change points among the time points in the first row of Gnu.
#' @param custom_partition A regular user defined partition. It is ignored if MinimalPartition is TRUE. 
#' @param verbose A Boolean. 
#' @return \item{est_K}{An integer representing the detected number of clusters.}
#' @return \item{est_z}{An integer vector of length N (the number of nodes). The i-th entry of est_z is the estimated cluster of node i.}
#' @return \item{est_D}{An integer representing the detected number of segments.}
#' @return \item{est_cp}{A vector of doubles containing the D estimated change points (the last entry is equal to the last point of the partition).}
#' @return \item{est_Lbd}{ A (K x K x D) tensor. Its entries are the estimated non-homogeneous Poisson intensities.}
#' @return \item{est_pi}{A vector of length K containing the estimated mixing proportions.}

ModFit <- function(Gnu,
                   kmin,
                   kmax,
                   eps = 10^(-9),
                   iter.vem = 50,
                   N_initializations = 1,
                   MinimalPartition = TRUE,
                   custom_partition = NULL, 
                   verbose = FALSE
                   ){
    InitN <- N <- max(max(Gnu[2,]), max(Gnu[3,]))    # for the moment InitN and N are the same thing.. 
    M <- ncol(Gnu)
    if (MinimalPartition==TRUE){
      X <- array(0, dim = c(InitN,InitN,M))
      for (i in 1:M) X[as.integer(Gnu[2,i]), as.integer(Gnu[3,i]), i]<-X[as.integer(Gnu[3,i]), as.integer(Gnu[2,i]), i]<-1
    }
    else {
      ptn <- custom_partition
      X <- GnuToX(EGnu=Gnu, Eptn = ptn, N = InitN)
    }
    tmp_D<-c()
    tmp_LB<-c()
    
    # Init (Fabrice k-means)
    # Building an initial partition (ipart) possibly having nothing in common with custom partition. This partition
    # will be used for initialization purposes only
    istep <- round((Gnu[1,M] - Gnu[1,1])/100)
    ipart <- seq(from = Gnu[1,1], to = Gnu[1,M], by = istep)
    X_ <- GnuToX(EGnu = Gnu, Eptn=ipart, N=InitN)
    X2<-X_
    dim(X2)<-c(InitN, InitN*length(ipart))
    lock_tau <- vector("list", length = kmax)
    lock_ptn <- vector("list", length = kmax)
    lock_Lbd <- vector("list", length = kmax)
    lock_pi <- vector("list", length = kmax)
    ## loop in K
    for (K in kmin:kmax){
    if (verbose == TRUE) cat(" testing for number of clusters K equal to: ", K, "\n")
    if (K==1){ 
      inz<-rep(1, times=N)
      tau<-matrix(0, N, K)  
      for (k in 1:K){
        store<-which(inz==k)
        tau[store,k]<-1
      }
      tau[tau<.Machine$double.xmin] <- .Machine$double.xmin
      
      # firt Maximization step
      if (MinimalPartition==TRUE){
        out<-VM_MP(Etau =tau, ptn=Gnu[1,], EGnu = Gnu)
        out$Lbd[out$Lbd<.Machine$double.xmin]<- .Machine$double.xmin
        LB<-LowerBound_MP(tau, Gnu, out$ptn, out$Pi, out$Lbd)
        if (verbose==TRUE) print(LB)
      }
      if (MinimalPartition==FALSE){
        out<-VM_(Etau = tau, ptn=ptn, EX=X)
        out$Lbd[out$Lbd<.Machine$double.xmin]<-.Machine$double.xmin
        LL<-log(out$Lbd)
        LB<-LowerBound_(tau, X, out$ptn[-1], out$Pi, out$Lbd,LL,
                        ValuesAtPtn = ptn[out$ptn])
        if (verbose==TRUE) cat(LB)
      }
      
      counter<-0
      old_LB=LB-1
      while (abs(LB-old_LB)>eps && counter<iter.vem){
        if (LB - old_LB < 0.001) print("R::warning: decreasing lower bound!!!")
        if (LB == .Machine$double.xmin) print("R::Warning: NaN LB")
        counter <- counter+1
        old_LB <- LB
        store_pi <- out$Pi
        
        # variational Expectation
        Iptn<-out$ptn[-1]  # initial changepoints
        ID<-length(Iptn)   # initial D
        IL<-out$Lbd    
        ipi<-out$Pi
        
        if (MinimalPartition==FALSE){
          tau<-VE_(tau, Iptn,ValuesAtPtn = Iptn, ipi, IL, X)
          tau[tau<.Machine$double.xmin] <-.Machine$double.xmin
          LB<-LowerBound_(tau, X, out$ptn[-1], out$Pi, out$Lbd, LL,ValuesAtPtn = ptn[out$ptn])
          if (verbose == TRUE) print(LB)
          
          # Variational Maximization
          out<-VM_(Etau = tau, ptn=ptn, EX = X)
          out$Lbd[out$Lbd<.Machine$double.xmin]<-.Machine$double.xmin
          LL<-log(out$Lbd)
          LB<-LowerBound_(tau, X, out$ptn[-1], out$Pi, out$Lbd,LL,ValuesAtPtn = ptn[out$ptn])
          if (verbose == TRUE) print(LB)
        }
        if (is.nan(LB)) LB<-.Machine$double.xmin
        if (MinimalPartition==TRUE){
          # No Variational Expetctation: there is a single cluster
          tau[tau<.Machine$double.xmin] <- .Machine$double.xmin
          # Variational Maximization
          out<-VM_MP(Etau = tau, ptn=Gnu[1,], EGnu = Gnu)
          out$Lbd[out$Lbd<.Machine$double.xmin]<- .Machine$double.xmin
          LB<-LowerBound_MP(tau, Gnu, out$ptn, out$Pi, out$Lbd)
          if (verbose == TRUE) print(LB)
        }
      }
      if (MinimalPartition==TRUE) tmpo<-LowerBound_MP(tau, Gnu, out$ptn, out$Pi, out$Lbd)
      else{
        LL<-log(out$Lbd)
        tmpo<-LowerBound_(tau,X, out$ptn[-1], out$Pi, out$Lbd,LL, ValuesAtPtn = ptn[out$ptn])
      }
      
      # stocking the estimated values of D and the lower bound for K=1
      tmp_D<-c(tmp_D, length(out$ptn[-1]))
      tmp_LB<-c(tmp_LB, tmpo)
      
      # stocking tau, the estimated changepoints and the model parameters
      lock_tau[[K-kmin+1]] <- tau
      lock_ptn[[K-kmin+1]] <- out$ptn
      lock_Lbd[[K-kmin+1]] <- out$Lbd
      lock_pi[[K-kmin+1]] <- out$Pi
    }
    else {
      # we produce 10 different initializations with the spectral clustering
      INz<-matrix(nrow=N_initializations, ncol=N)
      for (ctr in 1:N_initializations) INz[ctr, ]<-kmeans(X2, K)$cluster
      # I am selecting the equivalent (due to label switching) initial taus
      INz<-CleanLines(INz) 
      if (!is.matrix(INz)) INz<-t(as.matrix(INz))
      # I am running the algorithm for each initialization
      Intau<-vector("list", length = nrow(INz))
      Inout<-vector("list", length = nrow(INz))
      InLB<-vector(length=nrow(INz))
      # Running VEM for each initialization
      for (RowCounter in 1:nrow(INz)){
        inz<-INz[RowCounter,]
        ## init tau
        tau<-matrix(0, N, K)  
        for (k in 1:K){
          store<-which(inz==k)
          tau[store,k]<-1
        }
        tau[tau<.Machine$double.xmin] <- .Machine$double.xmin
        
        ## first M step
        if (MinimalPartition==FALSE) out<-VM_(Etau = tau, ptn=ptn, EX=X)
        else out<-VM_MP(Etau = tau, ptn=Gnu[1,], EGnu = Gnu)
        out$Lbd[out$Lbd<.Machine$double.xmin]<-.Machine$double.xmin
        LL<-log(out$Lbd)
        if (MinimalPartition==FALSE) LB<-LowerBound_(tau, X, out$ptn[-1], out$Pi, out$Lbd,LL,
                        ValuesAtPtn = ptn[out$ptn])
        else LB <- LowerBound_MP(tau, Gnu, out$ptn, out$Pi, out$Lbd)
        if (verbose == TRUE) cat("M step: ", LB, "\n")
        
        ## VEM loop
        counter<-0
        old_LB=LB-1
        while (abs(LB-old_LB)>eps && counter< iter.vem){
          if (LB<old_LB && LB!=.Machine$double.xmin) { 
            print("R::warning: decreasing lower bound !!!")
          }
          if (LB==.Machine$double.xmin) print("R::Warning: NaN LB")
          counter<-counter+1
          old_LB<-LB
          store_pi<-out$Pi
          
          # Variational E step
          Iptn<-out$ptn[-1]
          ID<-length(Iptn)
          IL<-out$Lbd
          ipi<-out$Pi
          if (MinimalPartition==FALSE) tau <- VE_(tau, Iptn, ValuesAtPtn = ptn[Iptn], ipi, IL, X)
          if (MinimalPartition==TRUE) tau <- VE_MP(tau, Iptn, ipi, IL, X, Gnu)
          tau[tau<.Machine$double.xmin] <-.Machine$double.xmin;
          if (MinimalPartition==FALSE) LB<-LowerBound_(tau, X, out$ptn[-1], out$Pi, out$Lbd, LL,ValuesAtPtn = ptn[out$ptn])
          else LowerBound_MP(tau, Gnu, out$ptn, out$Pi, out$Lbd)
          if (verbose==TRUE) cat("E step: ",LB, "\n")
          
          # Variational M step
          if (MinimalPartition==FALSE) out<-VM_(Etau = tau, ptn=ptn, EX = X)
          else out <- VM_MP(tau, ptn=Gnu[1,], EGnu = Gnu)
          out$Lbd[out$Lbd<.Machine$double.xmin]<-.Machine$double.xmin
          LL<-log(out$Lbd)
          if (MinimalPartition == FALSE) LB<-LowerBound_(tau, X, out$ptn[-1], out$Pi, out$Lbd,LL,ValuesAtPtn =ptn[out$ptn])
          else LB <- LowerBound_MP(tau, Gnu, out$ptn, out$Pi, out$Lbd)
          if (is.nan(LB)){
            print(c("Warning: NaN lower bound "))
            LB<-.Machine$double.xmin
          }
          if (verbose == TRUE) cat("M step: ",LB, "\n")
        }
        
        ## Last lower bound
        LL<-log(out$Lbd)
        if (MinimalPartition==FALSE) tmpo<-LowerBound_(tau,X, out$ptn[-1], out$Pi, out$Lbd,LL, ValuesAtPtn = ptn[out$ptn])
        else tmpo <- LowerBound_MP(tau, Gnu, out$ptn, out$Pi, out$Lbd)
        ## I'm saving what needed
        Intau[[RowCounter]]<-tau
        Inout[[RowCounter]]<-out
        InLB[RowCounter]<-tmpo
      }
      BestInit<-which.max(InLB) 
      tmp_D <- c(tmp_D, length(Inout[[BestInit]]$ptn[-1]))
      tmp_LB <- c(tmp_LB, InLB[BestInit])
      #
      lock_tau[[K - kmin + 1]] <- Intau[[BestInit]]
      lock_ptn[[K - kmin + 1]] <- Inout[[BestInit]]$ptn
      lock_Lbd[[K - kmin + 1]] <- Inout[[BestInit]]$Lbd
      lock_pi[[K - kmin +1]] <- Inout[[BestInit]]$Pi
    }
    }
    
    Kcount <- which.max(tmp_LB)              # selecting the best run
    estK <- Kcount + kmin - 1                # corresponding number of clusters
    estD <- tmp_D[Kcount]                    # corresponding segmentation
    estLbd <- lock_Lbd[[Kcount]]             # correspondint Lambda
    estPi <- lock_pi[[Kcount]]               # corresponding Pi
    N <- nrow(lock_tau[[Kcount]])
    
    # Building outz (the estimated z) 
    newMM<-MM <- matrix(0, nrow=N, ncol=estK)
    max_tau<-apply(lock_tau[[Kcount]], 1, max)
    MM[which(lock_tau[[Kcount]]==max_tau, arr.ind=1)] <- 1
    for (k in 1:estK) newMM[,k]<-k*MM[,k]
    outz<-apply(newMM, 1, sum)
        ## managing the inactive nodes if any
        #     if (length(inactive.nodes!=0)){
        #       tmp1<-seq(1:InitN);
        #       tmp2<-match(tmp1, WhoIsWho);
        #       tmp2[inactive.nodes]<-1;
        #       tmp2[WhoIsWho]<-outz
        #       outz<-tmp2
        #     }
    
    # Estimated change points
    tmp <- lock_ptn[[Kcount]]
    if (MinimalPartition == TRUE) est_eta <- Gnu[1, tmp]
    else est_eta <- ptn[tmp]
    return(list(est_K = estK, est_z = outz, est_D = estD, est_cp = est_eta, est_Lbd = estLbd, est_pi = estPi))
}