// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


double compute_logLike(const arma::colvec& Y, const arma::mat& X, const arma::colvec& beta, const std::vector<arma::mat>& vc, const arma::colvec& theta, bool REML = true) {

  int i;
  arma::uword n = X.n_rows, nVC = vc.size(); // k = X.n_cols, 
  arma::mat Omega(n,n);
  double logDet = 0, logDet2 = 0;
  double sign;
  double logLike;

  Omega.zeros();
  for (i=0; i<nVC; i++) {
    Omega += theta(i)*vc[i];
  }

  // Invert the estimated variance matrix
  // and compute the determinant
  arma::mat IOmega = inv_sympd(Omega);
  log_det(logDet, sign, Omega);

  if (REML) {
      arma::mat IOmegaX = (IOmega*X);
      arma::mat xOx = arma::trans(X) * (IOmegaX);
      arma::mat IxOmegax = arma::inv(xOx);
      log_det(logDet2, sign, xOx);
      arma::mat P = IOmega - (IOmegaX*IxOmegax)*arma::trans(IOmegaX);   
      arma::mat IOmega2 = P*P;
      
      arma::mat PY = P*Y;
      logLike = - 0.5*(n*log(2*arma::datum::pi) + logDet + logDet2 + arma::as_scalar(arma::trans(PY)*Y));
      logLike -= - 0.5*beta.size()*log(2*arma::datum::pi);
	
  } else {
    arma::colvec resid = Y - X*beta;
    logLike = - 0.5*(n*log(2*arma::datum::pi) + logDet + arma::as_scalar(arma::trans(resid)*IOmega*(resid)));
  }
  

  return logLike;
}


// Currently missing in the algorithm in the code below:
// OK. Step halving
// 2. transform VC coefficients to have non-negative values - waste of time
// 3. Choose between ML / REML
// 4. Add convergence tolerance
// 5. Clustered input


//' Maximize linear mixed model with user-specified covariance matrices
//'
//' @description Fast computation of simple regression slopes for each predictor represented by a column in a matrix
//' @param y A vector of outcomes.
//' @param x A design matrix of regressor variables. Must have the same number of rows as the length of y.
//' @param vc A list of fixed variance-covariance matrices
//' @param maxiter Maximum number of iterations
//' @param REML A logical that determines if restricted maximum likelihood (REML - the default) or maximum likelihood (ML) 
//' @param tolerance The Maximum number of iterations
//' @return A data frame with two variables: coefficients and stderr that gives the slope estimate and corresponding standard error for each column in x.
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @export
// [[Rcpp::export]]
List lmm_Maximize_cpp(NumericVector y,
		      NumericMatrix x,
		      List vc,
		      int maxiter = 25,
		      bool REML = true,
		      double tolerance = 0.000001,
		      bool reparam = false,
		      bool scale = false
		      ) {
  arma::uword n = x.nrow(), k = x.ncol(), nVC = vc.size();

  // Should do sanity checks
  
  arma::mat X(x.begin(), n, k, false);
  arma::colvec Y = Rcpp::as<arma::colvec>(y); // Making a full copy on purpose since I want to scale it later. Otherwise use Y(y.begin(), y.size(), false);
  arma::colvec theta = arma::zeros(nVC+1);
  
  // Convert List/VC to at list of arma matrices
  std::vector<arma::mat> VC(vc.size()+1);

  for (int i=0; i < nVC; i++) {
    Rcpp::NumericMatrix tmpcv = vc[i];
    VC[i] = Rcpp::as<arma::mat>(tmpcv);
  }
  VC[nVC] = arma::eye<arma::mat>(n,n);
  
  
  // Initial estimate for beta (from independence model)

  double scalingfactor = stddev(Y);
  if (scale) {
    Y /= scalingfactor;
  }

  arma::colvec beta = arma::solve(X, Y);
  arma::colvec deltaBeta;
  
  arma::colvec resid = Y - X*beta;
  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  double logLike = 0, oldLogLike = 0;
  double logDet = 0, logDet2 = 0;
  double sign;

  // Initialize thetas
  theta += 0.005;
  theta(nVC) = sig2;
  
  if (reparam)
    theta = log(theta);
  
  int errorcode = 25;

  arma::mat Omega(n,n);
  arma::mat IOmega, P, IxOmegax, xOx, IOmegaX, IOmega2, IFisher;
  arma::colvec mu = X * beta;
  arma::mat Fisher(nVC+1, nVC+1), InvFisher(nVC+1, nVC+1) ;
  arma::colvec Deriv(nVC+1);
  arma::colvec WorkingTheta(nVC+1), PY(n), WorkingBeta(k);

  int i, j, iter;

  // Main iteration loop
  for (iter=0; iter<maxiter; iter++) {

    // printf("Iter %d\n", iter);

    // Compute the Inverse variance matrix
    Omega.zeros();
    
    for (i=0; i<=nVC; i++) {
      if (reparam)
	Omega += exp(theta(i))*VC[i];
      else 
	Omega += theta(i)*VC[i];  
    }
    // Invert the estimated variance matrix
    // and compute the determinant
    IOmega = inv_sympd(Omega);
    log_det(logDet, sign, Omega);

    Deriv.zeros();
    Fisher.zeros();
    
    // X^t Omega X
    IOmegaX = (IOmega*X);
    xOx = arma::trans(X) * (IOmegaX);
    IxOmegax = arma::inv(xOx);    

    if (REML) {
      log_det(logDet2, sign , xOx);
      P = IOmega - (IOmegaX*IxOmegax)*arma::trans(IOmegaX);   
      IOmega2 = P*P;
      
      PY = P*Y;
      for (i = 0; i < nVC+1; i++) {
	if (reparam) {
	  Deriv(i) = (-arma::trace(P*VC[i]) + arma::as_scalar((arma::trans(PY))*VC[i]*(PY)))*exp(theta(i));
	}
	else {
	  Deriv(i) = - arma::trace(P*VC[i]) + arma::as_scalar((arma::trans(PY))*VC[i]*(PY));
	}
	
	for (j = i; j < nVC+1; j++) {
	  Fisher(i,j) = arma::trace(IOmega2*VC[i]*VC[j]);
	  // Just do the VC VC multiplications once and store them for speed?	
	  
	  if (reparam) {
	    if (i==j) {
	      Fisher(i,j) = Fisher(i,j)*exp(2*theta(i)) + Deriv(i);
	    } else {
	      Fisher(i,j) *= exp(theta(i)+theta(j));
	    }
	  }	
	  Fisher(j,i) = Fisher(i,j);
	}
      }
    } else {   // ML
      resid =  (Y-X*beta);
      deltaBeta =  arma::trans(IOmegaX) * resid;

      IOmega2 = IOmega*IOmega;

      resid = IOmega*resid; // Note this change!!
      for (i = 0; i < nVC+1; i++) {
	Deriv(i) = - arma::trace(VC[i]*IOmega) + arma::as_scalar(arma::trans(resid)*VC[i]*resid);// Note the change in definition of resid
	
	for (j = i; j < nVC+1; j++) {
	  Fisher(i, j) = arma::trace(VC[j]*VC[i]*IOmega2);
	  Fisher(j, i) = Fisher(i,j);
	}
      }
    }
    
    // Remember to scale the 1st and 2nd derivatives by 1/2
    Deriv  *= .5;
    Fisher *= .5;

    // Inverts the matrix
    IFisher = inv(Fisher);

    // Calculate the change in delta
    arma::colvec DeltaTheta = IFisher*Deriv;
    // Possibly do step-halving
    double step = 2.0;
    double newll = 0;
    int stephalvingcount =0;
    do {
      step /= 2;
      stephalvingcount += 1;
      WorkingTheta = theta + DeltaTheta*step;
      // WorkingBeta = beta + IxOmegax*deltaBeta;

      if (!reparam) {
	for (i = 0; i<WorkingTheta.size(); i++) {
	  if (WorkingTheta(i)<=0)
	    WorkingTheta(i)=.0001;
	}
      }      
      newll = compute_logLike(Y, X, WorkingBeta, VC, WorkingTheta, REML);
    }
    while (newll<oldLogLike & stephalvingcount < 20);

    theta = WorkingTheta;

    if (!REML) 
      beta += IxOmegax*deltaBeta;
    
    // logLike = - 0.5*(n*log(2*arma::datum::pi) + logDet + logDet2 + arma::as_scalar(arma::trans(PY)*Y));
    // logLike -= - 0.5*beta.size()*log(2*arma::datum::pi);
    logLike = newll;

    /*
    Rcout << theta << arma::endl;
    Rcout << "New log.like   " << logLike << arma::endl;
    Rcout << "Log likelihood " << oldLogLike << arma::endl;
    Rcout << "XXX likelihood " << compute_logLike(Y, X, beta, VC, theta) << arma::endl;
    Rcout << "Difference     " << 10000*(logLike - oldLogLike) << arma::endl;
    */

    if (iter>0 & std::abs(logLike-oldLogLike)<tolerance) {
      errorcode=0;
      break;
    }

    // Prepare for next iteration
    oldLogLike = logLike;
    
  }

  if (reparam)
    theta = exp(theta);

  if (scale) {
    beta *= scalingfactor;
    Y *= scalingfactor;
    theta *= (scalingfactor)*scalingfactor;
    logLike = compute_logLike(Y, X, beta, VC, theta, REML);
  }
  
  // create a new data frame and return it
  return List::create(Rcpp::Named("coefficients")=as<NumericVector>(wrap(beta)),
		      Rcpp::Named("sigmas")=theta,
		      Rcpp::Named("logLik")=logLike,
		      Rcpp::Named("convergence")=errorcode,
		      Rcpp::Named("VarMatrix")=IxOmegax,
		      Rcpp::Named("InvFisher")=IFisher,		      		      
		      Rcpp::Named("REML")=REML,
		      Rcpp::Named("iterations")=iter
		      );
  
}

