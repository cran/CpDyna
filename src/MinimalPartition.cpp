# include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace std;

///  SAMPLING

// to create an event like (\nu_m, i_m, j_m)
class Event{
public:  
  // instances
  double value;
  int from;
  int to;
public:
  // mehtod
  Event();
  Event(double time, int departure, int arrival): value(time), from(departure), to(arrival){}
  bool operator>= (const Event& E) const {return this->value >= E.value;}
  bool operator<(const Event& E) const { return !(*this>=E);}
};

// stepwise function
double lbd(const arma::vec& breaks, const arma::vec lambdas, double t){
  int N(lambdas.size());
  int ctr(0);
  while (t>breaks.at(ctr) && t<breaks.at(N-1)) ctr++;
  return lambdas(ctr);
}


// rejection method to simulate a NHPP trajectory
void One_Traj_MP(vector<Event>& nu,
                           const arma::vec& breaks, 
                           const arma::vec& lambdas, 
                           const int i, 
                           const int j){
 double U(breaks.at(breaks.size()-1));
 double lbar(lambdas.max());
 Event T(0,i,j);
 while (T.value<U){
   double tmp(R::rexp(1/lbar));
   T.value+=tmp;
   if (R::runif(0,1) < lbd(breaks,lambdas, T.value)/lbar) nu.push_back(T);
 }
}


//[[Rcpp::export]]
arma::mat GetGnu(const arma::vec z,
                 const arma::vec breaks,
                 const arma::vec lIN,
                 const arma::vec lOut){
  int N(z.size());
  vector<Event> Gnu;
  for (int i(1); i<=N; i++){
    for (int j(i+1); j<=N; j++){
      vector<Event> nu;
      if (z(i-1)==z(j-1)) One_Traj_MP(nu, breaks, lIN,i,j);
      else One_Traj_MP(nu,breaks, lOut,i,j);
      Gnu.insert(Gnu.end(), nu.begin(), nu.end());
    }
  }
  std::sort(Gnu.begin(), Gnu.end());
  arma::mat M(arma::zeros(3,Gnu.size()));
  for (auto it(Gnu.begin()); it!=Gnu.end(); ++it){
    arma::colvec tmp(3);
    tmp(0)=it->value;
    tmp(1)=it->from;
    tmp(2)=it->to;
    M.col(it-Gnu.begin())=tmp;
  }
return M;
} 


//NHPP increments
Rcpp::NumericVector& BuildCounts_C(Rcpp::NumericVector& out, // it will crash if out is shorter than ptn
                                   const arma::rowvec& nu, 
                                   const arma::vec& ptn){
  const int U(ptn.size());
  if (nu.empty()) out.fill(0);
  else{
    auto it=nu.begin();
    auto old_it=it;
    for (int i(0); i<U;i++){
      if ((*it)>ptn(i)) out[i]=0;
      else{
        while((*it)<ptn(i) && it<nu.end()) it++;
        out[i]=it-old_it;  
      }
      old_it=it;
    }
  }
  return out;
}


// script to get a tensor from Gnu, when switching to the minimal partition structure to the wider partition
//[[Rcpp::export]]
arma::cube GnuToX(Rcpp::NumericVector EGnu, 
                  Rcpp::NumericVector Eptn,
                  int N
){
 // Gnu
 Rcpp::IntegerVector dimEGnu(EGnu.attr("dim"));
 arma::mat Gnu(EGnu.begin(), dimEGnu[0], dimEGnu[1], false);
 // ptn
 arma::vec ptn(Eptn.begin(), Eptn.size(), false);
 
 int Uptn(ptn.size());
 Rcpp::NumericVector cts(Uptn);
 BuildCounts_C(cts,Gnu.row(0),ptn);
 arma::cube out(arma::zeros<arma::cube>(N,N,Uptn));
 unsigned int counter(0);
 for (int i(0); i<Uptn; i++){
   if (cts[i]!=0){
     for (int j(1); j<=cts[i]; j++){
       out(Gnu(1, counter+j-1)-1, Gnu(2, counter+j-1)-1, i)++;
       out(Gnu(2, counter+j-1)-1, Gnu(1, counter+j-1)-1, i)++;
     }
     counter+=cts[i];
   }
 }
 return out;
}


//------ MAXIMIZATION -------

// function to move (for each (k,g)) from the X tensor (N*N*U) to a vector Y_{kg}^u=\sum_{j>i}^N (\tau^{i,k}\tau^{j,g} + \tau^{i,g}\tau^{j,k})X[i,j,u]
// same function ToYkg for the "Minimal partition case
void ToYkg(Rcpp::NumericVector& out,
           const arma::mat& Events, 
           const arma::mat& tau, 
           const int k, 
           const int g){
 int U(Events.n_cols);
 if (g>k){
  for (int i(0); i<U; i++) out(i)=tau(Events.at(1,i)-1,k)*tau(Events.at(2,i)-1,g) + tau(Events.at(1,i)-1,g)*tau(Events.at(2,i)-1,k);
 }
 else{
  for (int i(0); i<U; i++) out(i)=tau(Events.at(1,i)-1,k)*tau(Events.at(2,i)-1,k);
 }
}



// function to create, for each pair (k,g) a matrix (à la LEBARBIER) containing at position (i,j) Y_{k,g}^{]i,j]} for j \geq than i
void ToMkg(arma::mat& out,
           Rcpp::NumericVector& Ykg){
  int U(Ykg.size());
  int i(0);
  out(0,0)=Ykg(0);
  for (int j(1); j<U; j++) out(i,j)=out(i,j-1)+Ykg(j);
  for (int i(1); i<U; i++){
    for (int j(i); j<U; j++) out(i,j)=out(i-1,j)-out(i-1, i-1);//   sum(Ykg(arma::span(i,j)));
  }
}

// function that return for a fixed pair (k,g) the "variational size" \sum_i \sum_{j>i} [\tau^{i,k}\tau^{j,g} + \tau^{i,g}\tau^{j,k}

double VS(const arma::mat& tau, const int k, const int g){
  double out(0);
  if (g>k){
    for (int i(0); i<tau.n_rows; i++){
      for (int j(i+1); j<tau.n_rows;j++) out+=tau(i,k)*tau(j,g)+ tau(i,g)*tau(j,k);
    }
  }
  else if (g==k){
    for (int i(0); i<tau.n_rows; i++){
      for (int j(i+1); j<tau.n_rows;j++) out+=tau(i,k)*tau(j,k);
    }
  }
  else Rcpp::Rcout<<"Wrong value of g, greater than k in function ToYkg"<<endl; 
  if (out==0) Rcpp::Rcout<<" warning! Null variational size "<<k<<" "<<g<<endl;
  return out;
}




//[[Rcpp::export]]
Rcpp::List VM_MP(Rcpp::NumericVector Etau,  // responsibilities output from VE
                 arma::vec ptn,             //a priori (1:U)
                 Rcpp::NumericVector EGnu,    // the Events matrix
                 bool PeltOrNot=true
){
  // tau
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  // Gnu
  Rcpp::IntegerVector dimEGnu(EGnu.attr("dim"));
  arma::mat Gnu(EGnu.begin(), dimEGnu[0], dimEGnu[1], true);
  //
  int K(tau.n_cols), N(tau.n_rows), U(Gnu.n_cols);

  ///---- step 1: estimation of break points
  
  // some useful data structures: for each pair (k,g) a matrix (à la LEBARBIER) containing at position (i,j) Y_{k,g}^{]i,j]} for j \geq than i
  // T: a map linking the pair (k,g) to its variational size
  std::map<pair<int, int>, shared_ptr<Rcpp::NumericVector> > Str;
  std::map<pair<int, int>, double> T;
   //cout<<"start building"<<endl;
  for (int k(0); k<K; k++){
    for (int g(k); g<K; g++){
      Rcpp::NumericVector Ykg(U);
      ToYkg(Ykg,Gnu,tau,k,g);
      pair<int, int> tmp(k,g);
      //arma::mat tmp_M(arma::zeros(U,U));
      shared_ptr<Rcpp::NumericVector> tmp_M(new Rcpp::NumericVector(Rcpp::cumsum(Ykg)));
      // old code
      //shared_ptr<arma::mat> tmp_M(new arma::mat(arma::zeros(U,U)));
      //ToMkg(*(tmp_M), Ykg);
      //ToMkg(tmp_M, Ykg);
      //shared_ptr<arma::sp_mat> p_tmp_M(new arma::sp_mat(tmp_M));
      //Str.insert(pair<pair<int,int>, shared_ptr<arma::mat> >(tmp, tmp_M));
      Str.insert(pair<pair<int,int>, shared_ptr<Rcpp::NumericVector> >(tmp, tmp_M));
      T.insert(pair< pair<int, int>, double >(tmp, VS(tau,k,g)));
    }
  }
  //cout<<"end building"<<endl;
  // opt --
  double beta(.25*K*(K+1)*log((N*(N-1)/2)*U));
  arma::vec F(U+1);
  // change points
  map<int, Rcpp::NumericVector> cp;
  cp[0].push_back(0);
  list<pair<int, double> > ptr;
  ptr.push_back(pair<int, double>(0,0.01));
  int eta_star;
  for (eta_star=1; eta_star<=U; eta_star++){
    //cout<<" "<<eta_star;
    if (ptr.empty()){ 
      Rcpp::Rcout<< " Warning!!!!!!!!!!!!!!!!!!!!!" <<endl;
      break;
    }
    double eta(0);
    double delta(ptn(eta_star-1)-eta);
    F(eta_star)=0;
    for (int k(0); k<K; k++){
      for (int g(k); g<K; g++){
        pair<int, int> tmp(k,g);  // to have access to the map
        if ((Str[tmp]->at(eta_star-1))!=0.0 && T[tmp]!=0){ // normally this should NOT arrive...
          double lbd_kgd(Str[tmp]->at(eta_star-1)/(delta*T[tmp]));
          F(eta_star)+=(-lbd_kgd*delta*T[tmp] + Str[tmp]->at(eta_star-1)*log(lbd_kgd));
        }
      }
    }
    double eta_prime(eta);  
    cp[eta_star]=cp[eta_prime];
    cp[eta_star].push_back(eta_prime);
    auto it1(ptr.begin());
    it1->second=F(eta_star);
    if ((eta_star-1)!=0) ptr.push_back(pair<int, double>((eta_star-1), 0.01));
    it1++;
    while (it1!=ptr.end()){
      delta=ptn(eta_star-1)-ptn((it1->first) - 1);
      double ll(0);
      //arma::mat StoreToPar(arma::zeros(K,K));
# pragma omp parallel for reduction (+:ll) num_threads(1)
      for (int k=0; k<K; k++){
        for (int g=k; g<K; g++){
          pair<int, int> tmp(k,g);  // to have access to the map
          if (it1->first==0){ // to avoid an at(-1)
            if (Str[tmp]->at(eta_star-1)==0 || T[tmp]==0) {ll+=0;}
            else{
              double lbd_kgd(Str[tmp]->at(eta_star-1)/(delta*T[tmp]));
              ll+=- lbd_kgd*delta*T[tmp] + (Str[tmp]->at(eta_star-1)*log(lbd_kgd));
              //StoreToPar(k,g)= -lbd_kgd*delta*T[tmp] + (Str[tmp]->at(eta_star-1)*log(lbd_kgd));
            }
          }
          else{
           if (Str[tmp]->at(eta_star-1)-Str[tmp]->at(it1->first -1)==0 || T[tmp]==0) {ll+=0;}
           else{
             double lbd_kgd((Str[tmp]->at(eta_star-1)- Str[tmp]->at(it1->first-1))/(delta*T[tmp]));
             ll+=- lbd_kgd*delta*T[tmp] + (Str[tmp]->at(eta_star-1)- Str[tmp]->at(it1->first-1))*log(lbd_kgd);
             //StoreToPar(k,g)=- lbd_kgd*delta*T[tmp] + (Str[tmp]->at(eta_star-1)- Str[tmp]->at(it1->first-1))*log(lbd_kgd);
             
           }
          }
        }
      }
      //ll+=arma::accu(StoreToPar);
      ll+=F(it1->first)-beta;
      it1->second=ll;  
      if (ll>F(eta_star)){
        F(eta_star)=ll;
        eta_prime=it1->first;                                
        cp[eta_star]=cp[eta_prime];
        cp[eta_star].push_back(eta_prime);        
      }
     it1++;  
     }
    if (ptr.size()>1 && PeltOrNot==true){
     bool idx(false);
     auto it2(ptr.begin());
     ++it2;  
     while (it2!=ptr.end()){
       if ((it2->second + beta) <= F(eta_star)){
         idx=true;
         it2=ptr.erase(it2); // minore o uguale perghé io sto selezinando i punti da levare!
       }
       else ++it2;
     }
    }
  }
  Rcpp::NumericVector bp(cp[eta_star-1]);
  bp.erase(0,1);
  bp.push_back(U);
  
  //--Lambdas
  //cout<<" Lambdas"<<endl;
  int D(bp.length()-1);
  arma::cube Lbd(K,K,D);
  for (int d(0); d<D; d++){
    double delta(ptn(bp[d+1]-1)-ptn(bp[d]));
    for (int k(0); k<K; k++){
      for (int g(k); g<K; g++){ 
        pair<int,int> tmp(k,g);
        if (k!=g){
          if (T[tmp]==0) Lbd(g,k,d)=Lbd(k,g,d)=0;
          else{
            if (bp[d]!=0) Lbd(g,k,d)=Lbd(k,g,d)=(Str[tmp]->at(bp[d+1]-1)- Str[tmp]->at(bp[d]-1))/(delta*T[tmp]);
            else Lbd(g,k,d)=Lbd(k,g,d)=(Str[tmp]->at(bp[d+1]-1))/(delta*T[tmp]);
          }
        }
        else{
          if (T[tmp]==0) Lbd(k,g,d)=0;
          else{
            if (bp[d]!=0) Lbd(k,g,d)=(Str[tmp]->at(bp[d+1]-1)-Str[tmp]->at(bp[d]-1))/(delta*T[tmp]);
            else Lbd(k,g,d)=(Str[tmp]->at(bp[d+1]-1))/(delta*T[tmp]);
          }
        }
      }
    }
  }
  //--- step 2: estimation of pi's
  
  Rcpp::NumericVector Pi;
  for (int k(0); k<K; k++) Pi.push_back(mean(tau.col(k)));  
  
  return Rcpp::List::create(Rcpp::Named("ptn")=bp,
                            Rcpp::Named("Pi")=Pi,
                            Rcpp::Named("Lbd")=Lbd
  );
  
  
}


// Nota: the version above (VM_MP_) is equivalent to VM_MP execpt for the use of a differnt
// data structure (tmp_M is a matrix) that is too heavy (pesante). Hence the above version
// works with simulated small data frame but induce an overflow if used with real data

//[[Rcpp::export]]
Rcpp::List VM_MP_(Rcpp::NumericVector Etau,  // responsibilities output from VE
                  arma::vec ptn,             //a priori (1:U)
                  Rcpp::NumericVector EGnu,    // the Events matrix
                  bool PeltOrNot=true
){
  // tau
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  // Gnu
  Rcpp::IntegerVector dimEGnu(EGnu.attr("dim"));
  arma::mat Gnu(EGnu.begin(), dimEGnu[0], dimEGnu[1], true);
  //
  int K(tau.n_cols), N(tau.n_rows), U(Gnu.n_cols);
  
  ///---- step 1: estimation of break points
  
  // some useful data structures: for each pair (k,g) a matrix (à la LEBARBIER) containing at position (i,j) Y_{k,g}^{]i,j]} for j \geq than i
  // T: a map linking the pair (k,g) to its variational size
  std::map<pair<int, int>, shared_ptr<arma::mat> > Str;
  std::map<pair<int, int>, double> T;
  Rcpp::Rcout<<"start building"<<endl;
  for (int k(0); k<K; k++){
    for (int g(k); g<K; g++){
      Rcpp::NumericVector Ykg(U);
      ToYkg(Ykg,Gnu,tau,k,g);
      pair<int, int> tmp(k,g);
      //arma::mat tmp_M(arma::zeros(U,U));
      //shared_ptr<Rcpp::NumericVector> tmp_M(new Rcpp::NumericVector(Rcpp::cumsum(Ykg)));
      // old code
      shared_ptr<arma::mat> tmp_M(new arma::mat(arma::zeros(U,U)));
      ToMkg(*(tmp_M), Ykg);
      //ToMkg(tmp_M, Ykg);
      //shared_ptr<arma::sp_mat> p_tmp_M(new arma::sp_mat(tmp_M));
      Str.insert(pair<pair<int,int>, shared_ptr<arma::mat> >(tmp, tmp_M));
      T.insert(pair< pair<int, int>, double >(tmp, VS(tau,k,g)));
    }
  }
  Rcpp::Rcout<<"end building"<<endl;
  // opt --
  double beta(.25*K*(K+1)*log((N*(N-1)/2)*U));
  arma::vec F(U+1);
  // change points
  map<int, Rcpp::NumericVector> cp;
  cp[0].push_back(0);
  list<pair<int, double> > ptr;
  ptr.push_back(pair<int, double>(0,0.01));
  int eta_star;
  for (eta_star=1; eta_star<=U; eta_star++){
    //cout<<" "<<eta_star;
    if (ptr.empty()){ 
      Rcpp::Rcout<< " Warning!!!!!!!!!!!!!!!!!!!!!" <<endl;
      break;
    }
    double eta(0);
    double delta(ptn(eta_star-1)-eta);
    F(eta_star)=0;
    for (int k(0); k<K; k++){
      for (int g(k); g<K; g++){
        pair<int, int> tmp(k,g);  // to have access to the map
      if ((Str[tmp]->at(0,eta_star-1))!=0.0 && T[tmp]!=0){ // normally this should NOT arrive...
        double lbd_kgd(Str[tmp]->at(0,eta_star-1)/(delta*T[tmp]));
        F(eta_star)+=(-lbd_kgd*delta*T[tmp] + Str[tmp]->at(0,eta_star-1)*log(lbd_kgd));
      }
    }
  }
    double eta_prime(eta);  
    cp[eta_star]=cp[eta_prime];
    cp[eta_star].push_back(eta_prime);
    auto it1(ptr.begin());
    it1->second=F(eta_star);
    if ((eta_star-1)!=0) ptr.push_back(pair<int, double>((eta_star-1), 0.01));
    it1++;
    while (it1!=ptr.end()){
      delta=ptn(eta_star-1)-ptn((it1->first) - 1);
      double ll(0);
      for (int k(0); k<K; k++){
        for (int g(k); g<K; g++){
          pair<int, int> tmp(k,g);  // to have access to the map
          if (Str[tmp]->at(it1->first,eta_star-1)==0 || T[tmp]==0) {ll+=0;}
          else{
            double lbd_kgd(Str[tmp]->at(it1->first,eta_star-1)/(delta*T[tmp]));
            ll+=- lbd_kgd*delta*T[tmp] + Str[tmp]->at(it1->first,eta_star-1)*log(lbd_kgd);
          }
        }
      }
      ll+=F(it1->first)-beta;
      it1->second=ll;  
      if (ll>F(eta_star)){
        F(eta_star)=ll;
        eta_prime=it1->first;                                
        cp[eta_star]=cp[eta_prime];
        cp[eta_star].push_back(eta_prime);        
      }
      it1++;  
    }
    if (ptr.size()>1 && PeltOrNot==true){
      bool idx(false);
      auto it2(ptr.begin());
      ++it2;  
      while (it2!=ptr.end()){
        if ((it2->second + beta) <= F(eta_star)){
          idx=true;
          it2=ptr.erase(it2); // minore o uguale perghé io sto selezinando i punti da levare!
        }
        else ++it2;
      }
    }
  }
  Rcpp::NumericVector bp(cp[eta_star-1]);
  bp.erase(0,1);
  bp.push_back(U);
  
  //--Lambdas
  
  int D(bp.length()-1);
  arma::cube Lbd(K,K,D);
  for (int d(0); d<D; d++){
    double delta(ptn(bp[d+1]-1)-ptn(bp[d]));
    for (int k(0); k<K; k++){
      for (int g(k); g<K; g++){ 
        pair<int,int> tmp(k,g);
        if (k!=g){
          if (T[tmp]==0) Lbd(g,k,d)=Lbd(k,g,d)=0;
          else Lbd(g,k,d)=Lbd(k,g,d)=Str[tmp]->at(bp[d],bp[d+1]-1)/(delta*T[tmp]);
        }
        else{
          if (T[tmp]==0) Lbd(k,g,d)=0;
          else Lbd(k,g,d)=Str[tmp]->at(bp[d],bp[d+1]-1)/(delta*T[tmp]);
        }
      }
    }
  }
  //--- step 2: estimation of pi's
  
  Rcpp::NumericVector Pi;
  for (int k(0); k<K; k++) Pi.push_back(mean(tau.col(k)));  
  
  return Rcpp::List::create(Rcpp::Named("ptn")=bp,
                            Rcpp::Named("Pi")=Pi,
                            Rcpp::Named("Lbd")=Lbd
  );
  
  
}




//--- In the next part we compute the lower bound for the likelihood 

// Entropia + pi
double Htp(const arma::mat &tau, const arma::vec& pi,  const int k){
  double out(0);
  const int N(tau.n_rows);
  for (int i(0); i<N; i++){
    //if (tau.at(i,k)!=0) out+=tau.at(i,k)*log(pi.at(k)/tau.at(i,k));
    if (tau.at(i,k)>1e-200) out+=tau.at(i,k)*log(pi.at(k)/tau.at(i,k));
    //cout<< "i "<<i<<"k "<<k<<" out "<<out<<endl;  
  }
  return out;
}

//[[Rcpp::export]]
double LowerBound_MP(Rcpp::NumericVector Etau,  // responsibilities output from VE (by ref.)
                     Rcpp::NumericVector EGnu,  // adjacency tensor (by ref.)
                     Rcpp::NumericVector ptn,   // bp output by VM
                     arma::vec Pi,              //  
                     Rcpp::NumericVector Lbd    // (by ref.) 
){
  // D and segments lengths
  int K(Pi.size());
  int D(ptn.size()-1);

  // Lambda
  Rcpp::IntegerVector dimLbd(Lbd.attr("dim"));
  arma::cube Lambda(Lbd.begin(), dimLbd[0], dimLbd[1], dimLbd[2], false); 
  // tau
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  int N(tau.n_rows);
  // Gnu
  Rcpp::IntegerVector dimEGnu(EGnu.attr("dim"));
  arma::mat Gnu(EGnu.begin(), dimEGnu[0], dimEGnu[1], false);
  int U(Gnu.n_cols);
  // breaks
  Rcpp::NumericVector breaks;
  breaks.push_back(0);
  auto it(ptn.begin());
  it++;
  while (it!=ptn.end()){
    breaks.push_back(Gnu(0, (*it)-1));
    it++;
  }
  //
  Rcpp::NumericVector delta_ptn(Rcpp::diff(breaks));

  // some useful data structures: for each pair (k,g) a matrix (à la LEBARBIER) containing at position (i,j) Y_{k,g}^{]i,j]} for j \geq than i
  // T: a map linking the pair (k,g) to its variational size
  std::map<pair<int, int>, shared_ptr<arma::mat> > Str;
  std::map<pair<int, int>, double> T;
  for (int k(0); k<K; k++){
    for (int g(k); g<K; g++){
      Rcpp::NumericVector Ykg(U);
      ToYkg(Ykg,Gnu,tau,k,g);
      pair<int, int> tmp(k,g);
      //arma::mat tmp_M(arma::zeros(U,U));
      shared_ptr<arma::mat> tmp_M(new arma::mat(arma::zeros(U,U)));
      ToMkg(*(tmp_M), Ykg);
      //ToMkg(tmp_M, Ykg);
      //shared_ptr<arma::sp_mat> p_tmp_M(new arma::sp_mat(tmp_M));
      Str.insert(pair<pair<int,int>, shared_ptr<arma::mat> >(tmp, tmp_M));
      T.insert(pair< pair<int, int>, double >(tmp, VS(tau,k,g)));
    }
  }
  //-- Lower bound
  double out(0); 
  for (int k(0); k<K; k++){
    for (int g(k); g<K; g++){ 
      double tmp1(0);
      double tmp2(0);
      for (int d(0); d<D; d++){
        //if (Lambda.at(k,g,d)!=0){ 
        //cout<<" ci passa "<<endl;
        tmp1+=Lambda.at(k,g,d)*delta_ptn(d); 
        tmp2+=log(Lambda.at(k,g,d))*Str[pair<int,int>(k,g)]->at(ptn(d),ptn(d+1)-1);
        //}
        //if (Lambda.at(k,g,d)==0 && Str[pair<int,int>(k,g)].at(ptn(d),ptn(d+1)-1)!=0) cout<<" warning: the lower bound is "<< log(Lambda.at(k,g,d))*Str[pair<int,int>(k,g)].at(ptn(d),ptn(d+1)-1)<<endl; 
      }
      out+= -tmp1*VS(tau,k,g) + tmp2;
    }
    out+=Htp(tau,Pi,k);
  }
  double alpha(.25*K*(K+1)*D*log((N*(N-1)/2)*U)), beta(.5*(K-1)*log(N));
  out-=(beta+alpha);
  return out;
}


// Expectation


///------ EXPECTATION -------------

// managing log-odds

void MLO(const arma::mat& log_p, arma::mat& tau, int idc, const int idr){
  //if (idc==1) cout<<" passo da qui !"<<idc<<endl;
  tau.at(idr,idc)=0;
  if (tau.n_cols==(idc+2)) tau(idr,idc+1)=1;
  else{
    double den(0);
    idc++;
    for (int ctr(idc); ctr<tau.n_cols; ctr++) den+=exp(log_p(idr,ctr)-log_p(idr,idc));
    if (isfinite(den)){
      for (int ctr(idc); ctr<tau.n_cols; ctr++) tau(idr,ctr)=exp(log_p(idr,ctr)-log_p(idr,idc))/den;
    }
    else MLO(log_p,tau,idc, idr); 
  }
}

// -- Weighted degree

double WD(const arma::mat& tau, int i_0, int d, int g, const arma::cube& X){// ci sono troppe copie qui dentro, da rivedere
  arma::mat tmp1(X.slice(d));
  //cout<<" WD 1"<<endl;
  arma::rowvec tmp2(tmp1.row(i_0));
  //cout<<" WD 2"<<endl;
  arma::colvec tmp3(tau.col(g));
  //cout<<" WD 3"<<endl;
  double out(arma::as_scalar(tmp2*tmp3)-tau(i_0,g)*X(i_0,i_0,d)); // en vrai on n'a pas besoin d'enlever tau*X(i_0,i_0,d) car il est nul..
  return(out);
}

// -- Counts aggregations over time segments (the break point lies in the previous segment)

/*vector<map<int, vector<int> > > SegmentDegree(const arma::mat& Gnu,
                                              const Rcpp::NumericVector& ptn
){
  int idx(0);
  for (int d(0); d<ptn.size(); d++){
   map<int, vector<int> > tmp;
   while (idx<ptn(d)){
     auto it1(tmp.find(Gnu(1,idx)));
     auto it2(tmp.find(Gnu(2,idx)));
     if (it1!=tmp.end() && it2!=tmp.end()){
       it1->second.push_back(it2->first);
       it2->second.push_back(it1->first); 
     }
     if (it1!=tmp.end() && it2==tmp.end()){
       it1->second.push_back(Gnu(2,idx));
       vector<int> tmp_vec;
       tmp_vec.push_back(it1->first);
       tmp.insert(pair<int, vector<int> >(Gnu(2,idx), tmp_vec));
     }
     if (it1==tmp.end() && it2!=tmp.end()){
       vector<int> tmp_vec;
       tmp_vec.push_back(it2->first);
       tmp.insert(pair<int, vector<int> >(Gnu(1,idx), tmp_vec));
       it2->second.push_back(Gnu(1,idx));
     }
     if (it1==tmp.end() && it2==tmp.end()){
       vector<int> tmp_vecA, tmp_vecB;
       tmp_vecA.push_back(Gnu(1, idx));
       tmp_vecB.push_back(Gnu(2, idx));
       tmp.insert(pair<int, vector<int> >(Gnu(1, idx), tmp_vecB));
       tmp.insert(pair<int, vector<int> >(Gnu(2, idx), tmp_vecA));
       
     };
   }  
  }
  
} */
  


arma::cube Aggreg(const arma::cube& X, const Rcpp::NumericVector& ptn){
  int d(0);
  arma::cube out(X.n_rows, X.n_cols, ptn.size());
  arma::mat tmp(arma::zeros(X.n_rows, X.n_cols));
  for (int ctr(0); ctr<ptn[d]; ctr++) tmp+=X.slice(ctr);
  out.slice(d)=tmp;
  for (d=1; d<ptn.size(); d++){
    tmp.zeros();
    for (int ctr(ptn[d-1]); ctr<ptn[d]; ctr++) tmp+=X.slice(ctr);
    out.slice(d)=tmp;
  }
  return out;
}

arma::mat VE_MP_(Rcpp::NumericVector Etau,  //(by val)
                 Rcpp::NumericVector ptn,   // by val
                 Rcpp::NumericVector Pi,    // by val
                 Rcpp::NumericVector Lbd,   // by re
                 Rcpp::NumericVector EGnu
){
  // D and segments lengths
  int K(Pi.size());
  int D(ptn.size());
  Rcpp::IntegerVector dimEGnu(EGnu.attr("dim"));
  arma::mat Gnu(EGnu.begin(), dimEGnu[0], dimEGnu[1], dimEGnu[2], false); // the final false serves to take the objects by reference from R!!!
  Rcpp::NumericVector breaks;
  for (auto it(ptn.begin()); it!=ptn.end(); ++it) breaks.push_back(Gnu(0, (*it)-1));
  //
  Rcpp::NumericVector delta_ptn(Rcpp::diff(breaks));
  delta_ptn.push_front(breaks.at(0));
  // Lambda
  Rcpp::IntegerVector dimLbd(Lbd.attr("dim"));
  arma::cube Lambda(Lbd.begin(), dimLbd[0], dimLbd[1], dimLbd[2], false); // the final false serves to take the objects by reference from R!!!
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  //
  int M(Gnu.n_cols);
  vector<arma::mat> tensor(M);
    for (int ctr(0); ctr<M; ctr++){
      arma::umat locations;
      locations << Gnu(1,ctr)-1 << Gnu(2,ctr)-1 << arma::endr
                << Gnu(2,ctr)-1 << Gnu(1,ctr)-1 << arma::endr;
      arma::vec values;
      values << 1 << 1 << arma::endr;
    tensor.at(ctr)=arma::sp_mat(locations,  values);
    }
   return tau;  
}

// variational expectation
// [[Rcpp::export]]
arma::mat VE_MP(Rcpp::NumericVector Etau,  //(by val)
                Rcpp::NumericVector ptn,   // by val
                Rcpp::NumericVector Pi,    // by val
                Rcpp::NumericVector Lbd,   // by ref
                Rcpp::NumericVector EX,
                Rcpp::NumericVector EGnu
                ){   
  // D and segments lengths
  int K(Pi.size());
  int D(ptn.size());
  Rcpp::IntegerVector dimEGnu(EGnu.attr("dim"));
  arma::mat Gnu(EGnu.begin(), dimEGnu[0], dimEGnu[1], dimEGnu[2], false); // the final false serves to take the objects by reference from R!!!
  Rcpp::NumericVector breaks;
  for (auto it(ptn.begin()); it!=ptn.end(); ++it) breaks.push_back(Gnu(0, (*it)-1));
  //
  Rcpp::NumericVector delta_ptn(Rcpp::diff(breaks));
  delta_ptn.push_front(breaks.at(0));
  // Lambda
  Rcpp::IntegerVector dimLbd(Lbd.attr("dim"));
  arma::cube Lambda(Lbd.begin(), dimLbd[0], dimLbd[1], dimLbd[2], false); // the final false serves to take the objects by reference from R!!!
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  // X: the original tesor (U slices)
  Rcpp::IntegerVector dimEX(EX.attr("dim"));
  arma::cube OX(EX.begin(), dimEX[0], dimEX[1], dimEX[2], false);  
  // X: aggregating over segments.
  arma::cube X(Aggreg(OX,ptn));
  // initialization of taus
  int N=X.slice(0).n_cols;
  arma::mat old_tau(N,K);
  old_tau.fill(1.0/tau.n_cols);
  double test(arma::accu(abs(tau-old_tau))); 
  int counter(0);
  int HowManyMLO(0);
  bool InMLO(false);
  while (test>0.000000000001 && counter<5){
    if (InMLO) {
      InMLO=false;
      HowManyMLO++;
    }
    counter++;
    //cout<<"test "<<test<<"   ";
    old_tau=tau;
    // preliminar data structures
    arma::mat VC(N,K);     // variational cardinality: at position (i_0, g) \sum_{i!=i_0} \tau^{(i,g)}
    std::vector<arma::mat> S; // at position i_0 there is a matrix containing at position (g,d) sum_{j!=i_0}\tau^{j,g} \DeltaN_{i_0,j}^d
    for (int i_0(0); i_0<N; i_0++){
      arma::mat tmp(arma::zeros(K,D));
      for (int g(0); g<K; g++){
        VC(i_0,g)=sum(tau.col(g))-tau(i_0,g);
        //if (VC(i_0, g)==0) cout<<" Warning 1! "<<endl;
        for (int d(0); d<D; d++) tmp(g,d)=WD(tau,i_0,d,g,X);
      }
      S.push_back(tmp);
    }
    arma::mat log_tau(arma::zeros(N,K));
    for (int i_0(0); i_0<N; i_0++){
      for (int k_0(0); k_0<K; k_0++){  
        double T1(0),T2(0);
        for(int g(0); g<K; g++){
          for (int d(0); d<D; d++){
            T1+=Lambda(k_0,g,d)*VC(i_0, g)*delta_ptn(d);
            T2+=log(Lambda(k_0,g,d))*S.at(i_0).at(g,d);
            if (Lambda(k_0,g,d)==0 && S.at(i_0).at(g,d)!=0) Rcpp::Rcout<<" so cazzi "<<endl;
          }
        }
        log_tau(i_0, k_0)=log(Pi(k_0)) - T1 + T2;
      }
    }
    // normalization
    for (int idr(0); idr<tau.n_rows; idr++){
      double den(0);
      for (int idc(0); idc<tau.n_cols; idc++) den+=exp(log_tau(idr,idc)-log_tau(idr,0));
      if (isfinite(den)){
        for (int idc(0); idc<tau.n_cols; idc++) tau(idr,idc)=exp(log_tau(idr, idc)-log_tau(idr,0))/den;
      }
      else {
        InMLO=true;
        //cout<<" MLO "<<endl;
        int idc(0);
        MLO(log_tau,tau,idc,idr);
      }
    }
    test=arma::accu(abs(tau-old_tau)); 
  }
  //cout<<" How many IMLO "<<HowManyMLO<<endl;
  return tau;
}  


/// Initialization

/*arma::mat GetDistance2(Rcpp::NumericVector EX){
    Rcpp::IntegerVector dimEX(EX.attr("dim"));
    arma::mat X(EX.begin(), dimEX[0], dimEX[1],false);
    const int N(X.n_rows), M(X.n_cols);
    arma::mat out(arma::zeros(N,N));
    for (int i(0); i<N; i++){
      for (int j(i+1); j<N; j++){
        for (int k(0); k<M; k++){
          if (k!=i && k!=j) out(j,i)+=pow((X(i,k)-X(j,k)),2);
        }
        out(j,i)=sqrt(out(j,i));
      }
    }
    return out;
  }*/

    ///------ SIMULATION functions -------

// stepwise function
double lbd_C(const arma::vec& breaks, const arma::vec lambdas, double t){
  int N(lambdas.size());
  int ctr(0);
  while (t>breaks(ctr) && t<breaks(N-1)) ctr++;
  return lambdas(ctr);
}

//NHPP increments
Rcpp::NumericVector& BuildCounts_C(Rcpp::NumericVector& out, const vector<double>& nu, const arma::vec& ptn){
  const int U(ptn.size());
  if (nu.empty()) out.fill(0);
  else{
    auto it=nu.begin();
    auto old_it=it;
    for (int i(0); i<U;i++){
      if ((*it)>ptn(i)) out[i]=0;
      else{
        while((*it)<ptn(i) && it<nu.end()) it++;
        out[i]=it-old_it;  
      }
      old_it=it;
    }
  }
  return out;
}

// One trajectory C++
vector<double>& One_Traj_C(vector<double>& nu,
                          const arma::vec& breaks, 
                          const arma::vec& lambdas, 
                          const int i, 
                          const int j){
  double U(breaks.at(breaks.size()-1));
  double lbar(lambdas.max());
  double T(0);
  while (T<U){
    double tmp(R::rexp(1/lbar));
    T+=tmp;
    if (R::runif(0,1) < lbd_C(breaks,lambdas, T)/lbar) nu.push_back(T);
  }
  return nu;
}



//[[Rcpp::export]]
arma::cube  GetX(const int U, 
                 const arma::vec ptn,
                 const arma::vec z,
                 const arma::vec breaks,
                 const arma::vec lIN,
                 const arma::vec lOut){
  int N(z.size());
  arma::cube X(N,N,U);
  X.fill(0);
  for (int i(0); i<(N-1); i++){
    for (int j(i+1); j<N; j++){
      vector<double> nu;
      if (z(i)==z(j)) nu=One_Traj_C(nu,breaks, lIN,i,j);
      else nu=One_Traj_C(nu,breaks, lOut,i,j);
      Rcpp::NumericVector tmp1(U);
      BuildCounts_C(tmp1,nu,ptn);
      for (int u(0); u<U; u++) X(j,i,u)=X(i,j,u)=tmp1(u);
    }
  }
  return X;
}


///------ EXPECTATION -------------

// managing log-odds

/*void MLO(const arma::mat& log_p, arma::mat& tau, int idc, const int idr){
  //if (idc==1) cout<<" passo da qui !"<<idc<<endl;
  tau.at(idr,idc)=0;
  if (tau.n_cols==(idc+2)) tau(idr,idc+1)=1;
  else{
    double den(0);
    idc++;
    for (int ctr(idc); ctr<tau.n_cols; ctr++) den+=exp(log_p(idr,ctr)-log_p(idr,idc));
    if (isfinite(den)){
      for (int ctr(idc); ctr<tau.n_cols; ctr++) tau(idr,ctr)=exp(log_p(idr,ctr)-log_p(idr,idc))/den;
    }
    else MLO(log_p,tau,idc, idr); 
  }
}*/

// -- Weighted degree

/*double WD(const arma::mat& tau, int i_0, int d, int g, const arma::cube& X){// ci sono troppe copie qui dentro, da rivedere
  arma::mat tmp1(X.slice(d));
  arma::rowvec tmp2(tmp1.row(i_0));
  //cout<<" WD 2"<<endl;
  arma::colvec tmp3(tau.col(g));
  //cout<<" WD 3"<<endl;
  double out(arma::as_scalar(tmp2*tmp3)-tau(i_0,g)*X(i_0,i_0,d)); // en vrai on n'a pas besoin d'enlever tau*X(i_0,i_0,d) car il est nul..
  return(out);
}*/

// -- Counts aggregations over time segments (the break point lies in the previous segment)

arma::cube& Aggreg(arma::cube& out, const arma::cube& X, const Rcpp::NumericVector& ptn){
  int d(0);
  //arma::cube out(X.n_rows, X.n_cols, ptn.size());
  arma::mat tmp(arma::zeros(X.n_rows, X.n_cols));
  //for(arma::cube::const_iterator it=X.begin(); it<(X.begin()+ptn[d]); ++it) tmp+=(*it);
  for (int ctr(0); ctr<ptn[d]; ctr++) tmp+=X.slice(ctr);
  out.slice(d)=tmp;
  for (d=1; d<ptn.size(); d++){
    tmp.zeros();
    for (int ctr(ptn[d-1]); ctr<ptn[d]; ctr++) tmp+=X.slice(ctr);
    out.slice(d)=tmp;
  }
return out;
}

// variational expectation
// [[Rcpp::export]]
arma::mat VE_(Rcpp::NumericVector Etau,  //(by val)
              Rcpp::NumericVector ptn,   // by val
              Rcpp::NumericVector ValuesAtPtn,
              Rcpp::NumericVector Pi,    // by val
              Rcpp::NumericVector Lbd,   // by ref
              Rcpp::NumericVector EX){   // by ref
  // D and segments lengths
  int K(Pi.size());
  int D(ptn.size());
  //Lambda
  Rcpp::IntegerVector dimLbd(Lbd.attr("dim"));
  arma::cube Lambda(Lbd.begin(), dimLbd[0], dimLbd[1], dimLbd[2], false); // the final false serves to take the objects by reference from R!!!
  // tau
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  // X: the original tesor (U slices)
  Rcpp::IntegerVector dimEX(EX.attr("dim"));
  arma::cube OX(EX.begin(), dimEX[0], dimEX[1], dimEX[2], false);  
  // X: aggregating over segments.
  arma::cube X(OX.n_rows, OX.n_cols, ptn.size());
  Aggreg(X, OX, ptn);
  //
  //ptn.push_front(0);
  ValuesAtPtn.push_front(0);
  arma::colvec delta_ptn(Rcpp::diff(ValuesAtPtn)); 
  //
  int N=X.slice(0).n_cols;
  int counter(0);
  int HowManyMLO(0);
  bool InMLO(false);
  while (/*test>0.00000000000000000000001 &&*/ counter<5){
    if (InMLO) {
      InMLO=false;
      HowManyMLO++;
    }
    counter++;
    arma::mat log_tau(arma::zeros(N,K));
    arma::rowvec ST(arma::sum(tau,0)); 
    int i_0;
    #pragma omp parallel for num_threads(20) //parallel for private(k_0)
    for (i_0=0; i_0<N; i_0++){
      int k_0;
       for (k_0=0; k_0<K; k_0++){
        double T1(0),T2(0), test_T2(0);
        arma::mat tmp1(Lambda.subcube(k_0,0,0,k_0,Lambda.n_cols-1, Lambda.n_slices-1));
        arma::mat tmp2(X.subcube(i_0,0,0,i_0,X.n_cols-1, X.n_slices-1));
        int i(0);
        if (D==1){ // se D=1, per qualche ragione (?) inverte le dimensioni delle matrici qui sopra 
          tmp1=tmp1.t();
          tmp2=tmp2.t();
        }
        arma::mat ttau(tau.t());
        /*while(i<N){ 
          if (i==i_0) i++;
          else{
           //cout<<" 11 "<<endl;
           T2+=arma::as_scalar(ttau.col(i).t()*log(tmp1)*(tmp2_.col(i)));
           //cout<<" 22 "<<endl;
           i++;
          }
        }*/
        T2 = arma::accu(log(tmp1) % (ttau*tmp2 - ttau.col(i_0)*tmp2.row(i_0)));
        T1 = arma::as_scalar((ST-tau.row(i_0))*tmp1*delta_ptn);
        log_tau(i_0, k_0)=log(Pi(k_0)) - T1 + T2;
      }
    }
    arma::colvec ml(arma::max(log_tau,  1));
    tau=ml*arma::ones(1, K);
    tau=log_tau-tau;
    tau=arma::sp_mat(exp(arma::mat(tau)));
    tau=tau/(arma::sum(tau, 1)*arma::ones(1,K));
  }
  return arma::mat(tau);
}  


//------ MAXIMIZATION -------

// this fucntion creates a vector of sparse Adjacency matrices from the data frame Gnu. 
// Hopefully useful in case the minimal partition is used with big networks..
vector<arma::sp_mat>& GnuToVecSPMat(vector<arma::sp_mat>& out,
                                    const arma::mat& Gnu,
                                    int NrowAndcols
                                    ){
  for (auto it=Gnu.begin(); it!=Gnu.end(); it+=3){
    // batch insertion of two values 
    arma::umat locations;
    locations << *(it+1)-1 << *(it+2)-1 << arma::endr
              << *(it+2)-1 << *(it+1)-1 << arma::endr;
    //
    arma::vec values;
    values << 1 << 1 <<arma::endr;
    //
    out.push_back(arma::sp_mat(locations, values,NrowAndcols, NrowAndcols));
  }
  return out;
}


// function to move (for each (k,g)) from the X tensor (N*N*U) to a vector Y_{kg}^u=\sum_{j>i}^N (\tau^{i,k}\tau^{j,g} + \tau^{i,g}\tau^{j,k})X[i,j,u]
arma::vec& ToYkg(arma::vec& out, 
                 const arma::cube& X,
                 const arma::mat& tau, 
                 const int k, 
                 const int g){
  if (k==g){
    for (int u(0); u<X.n_slices; u++){
      double tmp(0);
      for (int i(0); i<X.n_rows; i++){
        for (int j(i+1); j<X.n_cols; j++) tmp+=tau(i,k)*tau(j,k)*X(i,j,u);
      }
      out(u)=tmp;
    }
  }
  else if (g>k){
    for (int u(0); u<X.n_slices; u++){
      double tmp(0);
      for (int i(0); i<X.n_rows; i++){
        for (int j(i+1); j<X.n_cols; j++) tmp+=(tau(i,k)*tau(j,g) + tau(i,g)*tau(j,k))*X(i,j,u);
      }
      out(u)=tmp;
    }
  }
 else Rcpp::Rcout<<" Wrong value of g, greater than k in function ToYkg"<<endl;  
 return out;
}


// function to create, for each pair (k,g) a matrix (à la LEBARBIER) containing at position (i,j) Y_{k,g}^{]i,j]} for j \geq than i
void ToMkg(shared_ptr<arma::mat> out, const arma::vec& Ykg){
  int U(Ykg.size());
  for (int i(0); i<U; i++){
    for (int j(i); j<U; j++) out->at(i,j)=sum(Ykg(arma::span(i,j)));
  }
}

// function to aggregate adjacencies over all possible time segments
void ToM(map<pair<int, int>, arma::sp_mat>& out, const vector<arma::sp_mat>& inv){
  int U(inv.size());
  int i(0);
  pair<int, int> tmp(0,0);
  out.insert(pair< pair<int, int>, arma::sp_mat>(tmp, inv[0]));
  for (int j(1); j<U; j++){
    tmp=make_pair(i,j);
    out.insert(pair< pair<int, int>, arma::sp_mat>(tmp,(out[make_pair(i, j-1)]+inv[j])));
  }
 /* for (int i(1); i<U; i++){
    for (int j(i); j<U; j++){
      pair<int,int> tmp(i,j);
      out.insert(pair< pair<int, int>, arma::sp_mat >(tmp, out[make_pair(i-1,j)]-out[make_pair(i-1, i-1)]));      
    }
 }*/   
}


// function that return for a fixed pair (k,g) the "variational size" \sum_i \sum_{j>i} [\tau^{i,k}\tau^{j,g} + \tau^{i,g}\tau^{j,k}

/*double VS(const arma::mat& tau, const int k, const int g){
  double out(0);
  if (g>k){
    for (int i(0); i<tau.n_rows; i++){
      for (int j(i+1); j<tau.n_rows;j++) out+=tau(i,k)*tau(j,g)+ tau(i,g)*tau(j,k);
    }
  }
  else if (g==k){
    for (int i(0); i<tau.n_rows; i++){
      for (int j(i+1); j<tau.n_rows;j++) out+=tau(i,k)*tau(j,k);
    }
  }
  else cout<<"Wrong value of g, greater than k in function ToYkg"<<endl; 
  if (out==0) cout<<" warning! Null variational size "<<k<<" "<<g<<endl;
  return out;
}*/


// Maximization (old, not matrix form)

//[[Rcpp::export]]
Rcpp::List VM(Rcpp::NumericVector Etau, // responsibilities output from VE
              arma::vec ptn,            //a priori (1:U)
              Rcpp::NumericVector EX,    //adjacency tensor
              bool PeltOrNot=true  
              ){

      // tau
      Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
      arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
      // X
      Rcpp::IntegerVector dimEX(EX.attr("dim"));
      arma::cube X(EX.begin(), dimEX[0], dimEX[1], dimEX[2], false);
      //
      int K(tau.n_cols), N(tau.n_rows), U(X.n_slices);
      //double step(ptn(1)-ptn(0));
      
      ///---- step 1: estimation of break points
      
      // some useful data structures: for each pair (k,g) a matrix (à la LEBARBIER) containing at position (i,j) Y_{k,g}^{]i,j]} for j \geq than i
      // T: a map linking the pair (k,g) to its variational size
      std::map<pair<int, int>, shared_ptr<arma::mat> > Str;
      std::map<pair<int, int>, double> T;
        for (int k(0); k<K; k++){
        for (int g(k); g<K; g++){
          arma::vec Ykg(X.n_slices);
          ToYkg(Ykg,X,tau,k,g);
          pair<int, int> tmp(k,g);
          shared_ptr<arma::mat> tmp_M(new arma::mat(arma::zeros(U,U)));
          ToMkg(tmp_M, Ykg);
          Str.insert(pair<pair<int,int>, shared_ptr<arma::mat> >(tmp, tmp_M));
          T.insert(pair< pair<int, int>, double >(tmp, VS(tau,k,g)));
        }
      }
      // opt --
      double beta(.25*K*(K+1)*log((N*(N-1)/2)*U));
        arma::vec F(U+1);
        // change points
        map<int, Rcpp::NumericVector> cp;
        cp[0].push_back(0);
        list<pair<int, double> > ptr;
        ptr.push_back(pair<int, double>(0,0.01));
        int eta_star;
        for (eta_star=1; eta_star<=U; eta_star++){
          if (ptr.empty()){ 
            Rcpp::Rcout<< " Warning!!!!!!!!!!!!!!!!!!!!!" <<endl;
            break;
          }
          double eta(0);
          double delta(ptn(eta_star-1)-eta);
          F(eta_star)=0;
          for (int k(0); k<K; k++){
            for (int g(k); g<K; g++){
              pair<int, int> tmp(k,g);  // to have access to the map
              if ((Str[tmp]->at(0,eta_star-1))!=0.0 && T[tmp]!=0){ // normally this should NOT arrive...
                double lbd_kgd(Str[tmp]->at(0,eta_star-1)/(delta*T[tmp]));
                F(eta_star)+=(-lbd_kgd*delta*T[tmp] + Str[tmp]->at(0,eta_star-1)*log(lbd_kgd));
              }
            }
          }
          double eta_prime(eta);  
          cp[eta_star]=cp[eta_prime];
          cp[eta_star].push_back(eta_prime);
          auto it1(ptr.begin());
          it1->second=F(eta_star);
          if ((eta_star-1)!=0) ptr.push_back(pair<int, double>((eta_star-1), 0.01));
          it1++;
          while (it1!=ptr.end()){
            delta=ptn(eta_star-1)-ptn((it1->first) - 1);
            double ll(0);
            for (int k(0); k<K; k++){
              for (int g(k); g<K; g++){
                pair<int, int> tmp(k,g);  // to have access to the map
                if (Str[tmp]->at(it1->first,eta_star-1)==0 || T[tmp]==0) {ll+=0;}
                else{
                  double lbd_kgd(Str[tmp]->at(it1->first,eta_star-1)/(delta*T[tmp]));
                  ll+=- lbd_kgd*delta*T[tmp] + Str[tmp]->at(it1->first,eta_star-1)*log(lbd_kgd);
                }
              }
            }
            ll+=F(it1->first)-beta;
            it1->second=ll;  
            if (ll>F(eta_star)){
              F(eta_star)=ll;
              eta_prime=it1->first;                                
              cp[eta_star]=cp[eta_prime];
              cp[eta_star].push_back(eta_prime);        
            }
            it1++;  
          }
          if (ptr.size()>1 && PeltOrNot==true){
            bool idx(false);
            auto it2(ptr.begin());
            ++it2;  
            while (it2!=ptr.end()){
              if ((it2->second + beta) <= F(eta_star)){
                idx=true;
                it2=ptr.erase(it2); // minore o uguale perghé io sto selezinando i punti da levare!
              }
              else ++it2;
            }
          }
        }
        Rcpp::NumericVector bp(cp[eta_star-1]);
        bp.erase(0,1);
        bp.push_back(U);
/*        map<int, Rcpp::NumericVector> cp;
        cp[0].push_back(0);
        for (int eta_star(1); eta_star<=U; eta_star++){
          double eta(0);
          double delta(ptn(eta_star-1)-eta);
          F(eta_star)=0;
          for (int k(0); k<K; k++){
            for (int g(k); g<K; g++){
              pair<int, int> tmp(k,g);  // to have access to the map
              if ((Str[tmp]->at(0,eta_star-1))!=0.0 && T[tmp]!=0){ 
               double lbd_kgd(Str[tmp]->at(0,eta_star-1)/(delta*T[tmp]));
               F(eta_star)+=(-lbd_kgd*delta*T[tmp] + Str[tmp]->at(0,eta_star-1)*log(lbd_kgd));
              }
            }
          }
          //cout<<" occhio "<<F(eta_star)<<endl;
          double eta_prime(eta);  
          cp[eta_star]=cp[eta_prime];
          cp[eta_star].push_back(eta_prime);
          while (eta<eta_star-1){
            eta++;
            delta=ptn(eta_star-1)-ptn(eta-1);
            double ll(0);
            for (int k(0); k<K; k++){
              for (int g(k); g<K; g++){
                pair<int, int> tmp(k,g);  // to have access to the map
                if (Str[tmp]->at(eta,eta_star-1)==0 || T[tmp]==0) {ll+=0; //F(eta)-beta;
                }
                else{
                 double lbd_kgd(Str[tmp]->at(eta,eta_star-1)/(delta*T[tmp]));
                 ll+=- lbd_kgd*delta*T[tmp] + Str[tmp]->at(eta,eta_star-1)*log(lbd_kgd);
                }
              }
            }
            ll+=F(eta)-beta;
            if (ll>F(eta_star)){
              F(eta_star)=ll;
              eta_prime=eta;
              cp[eta_star]=cp[eta_prime];
              cp[eta_star].push_back(eta_prime);        
            }
          }
      }
      Rcpp::NumericVector bp(cp[U]);
      bp.erase(0,1);
      bp.push_back(U);*/

      //--Lambdas
      
      int D(bp.length()-1);
      arma::cube Lbd(K,K,D);
      for (int d(0); d<D; d++){
        double delta(bp[d+1]-bp[d]);
        for (int k(0); k<K; k++){
          for (int g(k); g<K; g++){ 
            pair<int,int> tmp(k,g);
            if (k!=g){
              if (T[tmp]==0) Lbd(g,k,d)=Lbd(k,g,d)=0;
              else Lbd(g,k,d)=Lbd(k,g,d)=Str[tmp]->at(bp[d],bp[d+1]-1)/(delta*T[tmp]);
            }
            else{
              if (T[tmp]==0) Lbd(k,g,d)=0;
              else Lbd(k,g,d)=Str[tmp]->at(bp[d],bp[d+1]-1)/(delta*T[tmp]);
            }
          }
        }
      }
      //--- step 2: estimation of pi's
      
      Rcpp::NumericVector Pi;
      for (int k(0); k<K; k++) Pi.push_back(mean(tau.col(k)));  
      
      return Rcpp::List::create(Rcpp::Named("ptn")=bp,
                                Rcpp::Named("Pi")=Pi,
                                Rcpp::Named("Lbd")=Lbd
                                );
    

}


/*
//[[Rcpp::export]]
Rcpp::List VM_(Rcpp::NumericVector Etau, // responsibilities output from VE
               arma::vec ptn,            //a priori (1:U)
               Rcpp::NumericVector EX,    //adjacency tensor
               bool PeltOrNot=true,
               bool MinimalPartition=true
){
  
  // tau
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  // X
  Rcpp::IntegerVector dimEX(EX.attr("dim"));
  arma::cube IX(EX.begin(), dimEX[0], dimEX[1], dimEX[2], false);
  //
  int K(tau.n_cols), N(tau.n_rows), U(IX.n_slices);
  
  // a more efficient structure for X:
  vector<arma::sp_mat> X;
  for (int u(0); u<U; u++)  X.push_back(arma::sp_mat(IX.slice(u)-arma::trimatl(IX.slice(u))));
  
  //cout<<"check 1"<<endl;
  
  // aggregating adjacencies in a map 
  map<pair<int, int>, arma::sp_mat> AgA;
  ToM(AgA, X);
  
  //cout<<" check 2"<<endl;
  
  // We calculate the denominator of L^{d}: \tau^{T} Ones(N) \tau + [\tau^{T} Ones(N,N) \tau]'
  arma::mat UU(arma::ones(N,N));
  UU=arma::mat(UU-arma::diagmat(UU));
  arma::mat DN(tau.t()*UU*tau);
  //DN+=DN.t();
  
  // opt --
  double beta(.25*K*(K+1)*log((N*(N-1)/2)*U));
  arma::vec F(U+1);
  // change points
  map<int, Rcpp::NumericVector> cp;
  cp[0].push_back(0);
  list<pair<int, double> > ptr;
  ptr.push_back(pair<int, double>(0,0.01));
  int eta_star;
  for (eta_star=1; eta_star<=U; eta_star++){
    if (ptr.empty()){ 
      cout<< " Warning!!!!!!!!!!!!!!!!!!!!!" <<endl;
      break;
    }
    double eta(0);
    double delta(ptn(eta_star-1)-eta);
    F(eta_star)=0;
    // lambda star
    arma::mat tmp(tau.t()*AgA[make_pair(eta, eta_star-1)]*tau);
    arma::mat NM(tmp + tmp.t());
    arma::mat L(NM/(delta*DN));
    // gain function
    F(eta_star)+=arma::accu(tmp % (log(L)-arma::ones(K,K)));
    double eta_prime(eta);  
    cp[eta_star]=cp[eta_prime];
    cp[eta_star].push_back(eta_prime);
    auto it1(ptr.begin());
    it1->second=F(eta_star);
    if ((eta_star-1)!=0) ptr.push_back(pair<int, double>((eta_star-1), 0.01));
    it1++;
    while (it1!=ptr.end()){
      delta=ptn(eta_star-1)-ptn((it1->first) - 1);
      //Lambda
      tmp=tau.t()*AgA[make_pair(it1->first, eta_star-1)]*tau;
      NM=tmp + tmp.t();
      L=NM/(delta*DN);
      // gain function
      double ll(arma::accu(tmp%(log(L)-arma::ones(K,K))));
      ll+=F(it1->first)-beta;
      it1->second=ll;  
      if (ll>F(eta_star)){
        F(eta_star)=ll;
        eta_prime=it1->first;                                
        cp[eta_star]=cp[eta_prime];
        cp[eta_star].push_back(eta_prime);        
      }
      it1++;  
    }
    // implementing pelt
    if (ptr.size()>1 && PeltOrNot==true){
      bool idx(false);
      auto it2(ptr.begin());
      ++it2;  
      while (it2!=ptr.end()){
        if ((it2->second + beta) <= F(eta_star)){
          idx=true;
          it2=ptr.erase(it2); // minore o uguale perghé io sto selezinando i punti da levare!
        }
        else ++it2;
      }
    }
  }//end for
  Rcpp::NumericVector bp(cp[eta_star-1]);
  bp.erase(0,1);
  bp.push_back(U);  

  //--Lambdas

  int D(bp.length()-1);
  arma::cube Lbd(K,K,D);
  for (int d(0); d<D; d++){
    double delta(0);
    if (MinimalPartition==false) delta=bp[d+1]-bp[d]; //  - ******************** !!!!
    else {
      if (d==0) delta=ptn(bp[d+1]-1);
      else delta=ptn(bp[d+1]-1) - ptn(bp[d]-1);
       
    }
    // L^{d}
    arma::mat NM(tau.t()*AgA[make_pair(bp[d], bp[d+1]-1)]*tau);
    NM+=NM.t();
    Lbd.slice(d)=NM/(delta*DN);
  }
  //--- step 2: estimation of pi's
  
  Rcpp::NumericVector Pi;
  for (int k(0); k<K; k++) Pi.push_back(mean(tau.col(k)));  
  
  return Rcpp::List::create(Rcpp::Named("ptn")=bp,
                            Rcpp::Named("Pi")=Pi,
                            Rcpp::Named("Lbd")=Lbd
  );
}
*/

// Uses a slowrer but slighter data strucutre
//[[Rcpp::export]]
Rcpp::List VM_(Rcpp::NumericVector Etau, // responsibilities output from VE
               arma::vec ptn,            //a priori (1:U)
               Rcpp::NumericVector EX,
               bool PeltOrNot=true,
               bool MinimalPartition=true
){
  
  // tau
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  // X
  Rcpp::IntegerVector dimEX(EX.attr("dim"));
  arma::cube IX(EX.begin(), dimEX[0], dimEX[1], dimEX[2], false);
  
  //Rcpp::IntegerVector dimGnu(Gnu.attr("dim"));
  //arma::mat Gnu_(Gnu.begin(), dimGnu[0], dimGnu[1], false);
  
  //
  int K(tau.n_cols), N(tau.n_rows), U(IX.n_slices);
  
  // a more efficient structure for X:
  vector<arma::sp_mat> X;
  //GnuToVecSPMat(X, Gnu, N);
  
  //cout<<"new code ok!!"<<endl;
  
  for (int u(0); u<U; u++)  X.push_back(arma::sp_mat(IX.slice(u)-arma::trimatl(IX.slice(u))));
  
  //cout<<"check 1"<<endl;
  
  // aggregating adjacencies in a map 
  map<pair<int, int>, arma::sp_mat> AgA;
  ToM(AgA, X);
  //cout<<" phase 1: end"<<endl;
  //cout<<" check 2"<<endl;
  
  // We calculate the denominator of L^{d}: \tau^{T} Ones(N) \tau + [\tau^{T} Ones(N,N) \tau]'
  arma::mat UU(arma::ones(N,N));
  UU=arma::mat(UU-arma::diagmat(UU));
  arma::mat DN(tau.t()*UU*tau);
  //DN+=DN.t();
  
  // opt --
  double beta(.25*K*(K+1)*log((N*(N-1)/2)*U));
  arma::vec F(U+1);
  // change points
  map<int, Rcpp::NumericVector> cp;
  cp[0].push_back(0);
  list<pair<int, double> > ptr;
  ptr.push_back(pair<int, double>(0,0.01));
  int eta_star;
  for (eta_star=1; eta_star<=U; eta_star++){
    //cout<<" "<<eta_star;
    if (ptr.empty()){ 
      Rcpp::Rcout<< " Warning!!!!!!!!!!!!!!!!!!!!!" <<endl;
      break;
    }
    double eta(0);
    double delta(ptn(eta_star-1)-eta);
    F(eta_star)=0;
    // lambda star
    arma::mat tmp(tau.t()*AgA[make_pair(eta, eta_star-1)]*tau);
    arma::mat NM(tmp + tmp.t());
    arma::mat L(NM/(delta*DN));
    // gain function
    F(eta_star)+=arma::accu(tmp % (log(L)-arma::ones(K,K)));
    double eta_prime(eta);  
    cp[eta_star]=cp[eta_prime];
    cp[eta_star].push_back(eta_prime);
    auto it1(ptr.begin());
    it1->second=F(eta_star);
    if ((eta_star-1)!=0) ptr.push_back(pair<int, double>((eta_star-1), 0.01));
    it1++;
    while (it1!=ptr.end()){
      delta=ptn(eta_star-1)-ptn((it1->first) - 1);
      //Lambda  /******************/
      arma::sp_mat TM;
      if (it1->first !=0) TM=(AgA[make_pair(0, eta_star-1)] - AgA[make_pair(0, it1->first-1)]);
      else TM=AgA[make_pair(it1->first, eta_star-1)];
      tmp=tau.t()*TM*tau;
      NM=tmp + tmp.t();
      L=NM/(delta*DN);
      // gain function
      double ll(arma::accu(tmp%(log(L)-arma::ones(K,K))));
      ll+=F(it1->first)-beta;
      it1->second=ll;  
      if (ll>F(eta_star)){
        F(eta_star)=ll;
        eta_prime=it1->first;                                
        cp[eta_star]=cp[eta_prime];
        cp[eta_star].push_back(eta_prime);        
      }
      it1++;  
    }
    // implementing pelt
    if (ptr.size()>1 && PeltOrNot==true){
      bool idx(false);
      auto it2(ptr.begin());
      ++it2;  
      while (it2!=ptr.end()){
        if ((it2->second + beta) <= F(eta_star)){
          idx=true;
          it2=ptr.erase(it2); // minore o uguale perghé io sto selezinando i punti da levare!
        }
        else ++it2;
      }
    }
  }//end for
  Rcpp::NumericVector bp(cp[eta_star-1]);
  bp.erase(0,1);
  bp.push_back(U);  
  
  //--Lambdas
  
  int D(bp.length()-1);
  arma::cube Lbd(K,K,D);
  for (int d(0); d<D; d++){
    double delta(0);
    if (MinimalPartition==false) delta=bp[d+1]-bp[d]; 
    else {
      if (d==0) delta=ptn(bp[d+1]-1);
      else delta=ptn(bp[d+1]-1) - ptn(bp[d]-1);
      
    }
    // L^{d}
    arma::sp_mat TM;
    if (bp[d]!=0) TM=AgA[make_pair(0, bp[d+1]-1)]-AgA[make_pair(0, bp[d]-1)];
    else TM=AgA[make_pair(bp[d], bp[d+1]-1)];
    arma::mat NM(tau.t()*TM*tau);
    NM+=NM.t();
    Lbd.slice(d)=NM/(delta*DN);
  }
  //--- step 2: estimation of pi's
  
  Rcpp::NumericVector Pi;
  for (int k(0); k<K; k++) Pi.push_back(mean(tau.col(k)));  
  
  return Rcpp::List::create(Rcpp::Named("ptn")=bp,
                            Rcpp::Named("Pi")=Pi,
                            Rcpp::Named("Lbd")=Lbd
  );
}

//--- In the next part we compute the lower bound for the likelihood 

// Entropia + pi
/*double Htp(const arma::mat &tau, const arma::vec& pi,  const int k){
  double out(0);
  const int N(tau.n_rows);
  for (int i(0); i<N; i++){
    //if (tau.at(i,k)!=0) out+=tau.at(i,k)*log(pi.at(k)/tau.at(i,k));
    if (tau.at(i,k)>1e-200) out+=tau.at(i,k)*log(pi.at(k)/tau.at(i,k));
    //cout<< "i "<<i<<"k "<<k<<" out "<<out<<endl;  
  }
  return out;
}*/

//[[Rcpp::export]]
double LowerBound(Rcpp::NumericVector Etau,  // responsibilities output from VE (by ref.)
                  Rcpp::NumericVector EX,    // adjacency tensor (by ref.)
                  Rcpp::NumericVector ptn,   // bp output by VM
                  arma::vec Pi,              //  
                  Rcpp::NumericVector Lbd    // (by ref.) 
                  ){
  // D and segments lengths
  int K(Pi.size());
  int D(ptn.size()-1);
  Rcpp::NumericVector delta_ptn(Rcpp::diff(ptn));
  // Lambda
  Rcpp::IntegerVector dimLbd(Lbd.attr("dim"));
  arma::cube Lambda(Lbd.begin(), dimLbd[0], dimLbd[1], dimLbd[2], false); 
  // tau
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat tau(Etau.begin(), dimEtau[0], dimEtau[1], true);
  int N(tau.n_rows);
  // X
  Rcpp::IntegerVector dimEX(EX.attr("dim"));
  arma::cube X(EX.begin(), dimEX[0], dimEX[1], dimEX[2], false);
  int U(X.n_slices);
  // some useful data structures: for each pair (k,g) a matrix (à la LEBARBIER) containing at position (i,j) Y_{k,g}^{]i,j]} for j \geq than i
  // T: a map linking the pair (k,g) to its variational size
  std::map<pair<int, int>, shared_ptr<arma::mat> >Str;
  std::map<pair<int, int>, double> T;
  for (int k(0); k<K; k++){
    for (int g(k); g<K; g++){
      arma::vec Ykg(U);
      ToYkg(Ykg, X,tau,k,g);
      pair<int, int> tmp(k,g);
      shared_ptr<arma::mat> tmp_M(new arma::mat(arma::zeros(U,U)));
      ToMkg(tmp_M, Ykg);
      Str.insert(pair<pair<int,int>, shared_ptr<arma::mat> >(tmp, tmp_M));
      T.insert(pair< pair<int, int>, double >(tmp, VS(tau,k,g)));
    }
  }
  //-- Lower bound
  double out(0); 
  for (int k(0); k<K; k++){
    for (int g(k); g<K; g++){ 
      double tmp1(0);
      double tmp2(0);
      for (int d(0); d<D; d++){
        //if (Lambda.at(k,g,d)!=0){ 
          //cout<<" ci passa "<<endl;
          tmp1+=Lambda.at(k,g,d)*delta_ptn(d);
          tmp2+=log(Lambda.at(k,g,d))*Str[pair<int,int>(k,g)]->at(ptn(d),ptn(d+1)-1);
        //}
        //if (Lambda.at(k,g,d)==0 && Str[pair<int,int>(k,g)].at(ptn(d),ptn(d+1)-1)!=0) cout<<" warning: the lower bound is "<< log(Lambda.at(k,g,d))*Str[pair<int,int>(k,g)].at(ptn(d),ptn(d+1)-1)<<endl; 
      }
      out+=-tmp1*VS(tau,k,g) +  tmp2;
    }
    //cout<<"T1: "<<out<<endl;
    out+=Htp(tau,Pi,k);
    //cout<<" Htp: "<<out<<endl;
  }
  double alpha(.25*K*(K+1)*D*log((N*(N-1)/2)*U)), beta(.5*(K-1)*log(N));
  out-=(alpha+beta);
  return out;
}




//[[Rcpp::export]]
double LowerBound_(Rcpp::NumericVector Etau,   // responsibilities output from VE (by ref.)
                   Rcpp::NumericVector EX,     // adjacency tensor (by ref.)
                   Rcpp::NumericVector ptn,    // bp output by VM
                   arma::vec Pi,               //  
                   Rcpp::NumericVector Lbd,    // (by ref.) 
                   Rcpp::NumericVector LLbd,   // (by ref.)
                   Rcpp::NumericVector ValuesAtPtn // =Rcpp::NumericVector::create()  
){   
  // D and segments lengths
  int K(Pi.size());
  int D(ptn.size());
  //Rcpp::NumericVector delta_ptn(Rcpp::diff(ptn));
  Rcpp::NumericVector delta_ptn;
  /*if (MinimalPartition==false){
    delta_ptn=Rcpp::diff(ptn);
    delta_ptn.push_front(*ptn.begin());
  }*/
  ValuesAtPtn.push_front(0);
  delta_ptn=Rcpp::diff(ValuesAtPtn);
  //delta_ptn.push_front(*ValuesAtPtn.begin());
  // Lambda
  Rcpp::IntegerVector dimLbd(Lbd.attr("dim"));
  arma::cube Lambda(Lbd.begin(), dimLbd[0], dimLbd[1], dimLbd[2], false);
  // LL
  Rcpp::IntegerVector dimLLbd(Lbd.attr("dim"));
  arma::cube LL(LLbd.begin(), dimLLbd[0], dimLLbd[1], dimLLbd[2], false);
  // tau
  Rcpp::IntegerVector dimEtau(Etau.attr("dim"));
  arma::mat Ttau(Etau.begin(), dimEtau[0], dimEtau[1], true); 
  int N(Ttau.n_rows);
  arma::sp_mat tau(Ttau);
  // X
  Rcpp::IntegerVector dimEX(EX.attr("dim"));
  arma::cube BX(EX.begin(), dimEX[0], dimEX[1], dimEX[2], false);
  int U(BX.n_slices);
  // aggregating X
  arma::cube X(N,N,D);
  Aggreg(X,BX,ptn);
  // T1
  arma::sp_mat L(arma::mat(K,K, arma::fill::zeros)) ;
  for (int d(0); d<D; d++) L+=delta_ptn[d]*arma::sp_mat(Lambda.slice(d));
  //double tmp(0);
  arma::sp_mat ntau(tau.t());
  arma::sp_mat tmp(tau*L*ntau);
  tmp-=arma::sp_mat(arma::diagmat(arma::mat(tmp)));
  double T1(arma::accu(tmp));
  T1=T1/2;
  //T2
  double T2(0);
  for (int d(0); d<D; d++){
    arma::sp_mat tmp(ntau*arma::sp_mat(X.slice(d)-arma::trimatl(X.slice(d)))*tau);
    tmp=LL.slice(d)%tmp;
    T2+=arma::accu(tmp);
  } 
  // Hentropy
  double T3(0);
  for (arma::sp_mat::iterator it=tau.begin(); it!=tau.end(); ++it) T3+=(*it)*log(Pi[it.col()]/(*it));
  double alpha((.25*K*(K+1)*D + .5*(K-1))*log((N*(N-1)/2)*U));// beta(.5*(K-1)*log(N));
  return -T1+T2+T3-alpha;  
}


// -- creating a distance matrix

//[[Rcpp::export]]
arma::mat GetDistance(Rcpp::NumericVector EX){
  // X
  Rcpp::IntegerVector dimEX(EX.attr("dim"));
  arma::cube X(EX.begin(), dimEX[0], dimEX[1], dimEX[2], false);
  //
  const int N(X.n_rows), U(X.n_slices);
  arma::mat out(arma::zeros(N, N));
  for (int i(0); i<N; i++){
    for (int j(i+1); j<N; j++ ){
      double d1(0);
      for (int u(0); u<U; u++){
        double d2(0);
        for (int k(0); k<N; k++){
           if (k!=i && k!=j) d2+=(abs((X(i,k,u) - X(j,k,u)))); 
        }
      d1+=d2;
      }
    out(j,i)=d1;  
    }
  }
  return out;
}



//[[Rcpp::export]]
arma::mat GetDistance2(Rcpp::NumericVector EX){
  Rcpp::IntegerVector dimEX(EX.attr("dim"));
  arma::mat X(EX.begin(), dimEX[0], dimEX[1],false);
  const int N(X.n_rows), M(X.n_cols);
  arma::mat out(arma::zeros(N,N));
  for (int i(0); i<N; i++){
    for (int j(i+1); j<N; j++){
      for (int k(0); k<M; k++){
        if (k!=i && k!=j) out(j,i)+=pow((X(i,k)-X(j,k)),2);
      }
      out(j,i)=sqrt(out(j,i));
    }
  }
  return out;
}






  
