// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fast marginal simple regresion analyses
//'
//' @description Fast computation of simple regression slopes for each predictor represented by a column in a matrix
//' @param y A vector of outcomes.
//' @param x A matrix of regressor variables. Must have the same number of rows as the length of y.
//' @param addintercept A logical that determines if the intercept should be included in all analyses (TRUE) or not (FALSE)
//' @return A data frame with two variables: coefficients and stderr that gives the slope estimate and corresponding standard error for each column in x.
//' @author Claus Ekstrøm <claus@@rprimer.dk>
//' @export
// [[Rcpp::export]]
DataFrame lmm_maximize_cpp(NumericVector y, NumericMatrix x, int addintercept) {
  arma::uword n = x.nrow(), k = x.ncol();
  int df = n-1;

  // Sanity checks
  if (y.size() != n) {
    stop("The length of y and the number of rows in x must match");
  }
  
  arma::mat X(x.begin(), n, k, false);
  arma::colvec Y(y.begin(), y.size(), false);
  arma::mat newX;
  arma::mat x0=arma::ones<arma::mat>(n,1);

  if (addintercept) {
    df = n-2;
  }

  arma::colvec rescoef = arma::zeros(k), resse = arma::zeros(k);
  arma::colvec coef, resid, stderrest;
  double sig2;
  
  for (arma::uword i=0; i<k; i++) {
    if (addintercept) {
      newX = join_rows(X.cols(i,i), x0);
    } else {
      newX = X.cols(i,i);
    }

    coef = arma::solve(newX, Y);
    rescoef(i) = coef(0);
    resid = Y - newX*coef;
    sig2 = arma::as_scalar(arma::trans(resid)*resid/df);
    stderrest = arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(newX)*newX)) );
    resse(i) = stderrest(0);
  }










  
  // create a new data frame and return it
  return DataFrame::create(Rcpp::Named("coefficients")=rescoef,
			   Rcpp::Named("stderr")=resse);
}




//' @export
// [[Rcpp::export]]
List lmm_Maximize_cpp(NumericVector y, NumericMatrix x, List vc, int maxiter) {
  arma::uword n = x.nrow(), k = x.ncol(), nVC = vc.size();

  // Should do sanity checks
  
  arma::mat X(x.begin(), n, k, false);
  arma::colvec Y(y.begin(), y.size(), false);
  arma::colvec theta = arma::zeros(nVC+1);

  // Convert List/VC to at list of arma matrices
  std::vector<arma::mat> VC(vc.size()+1);

  for (int i=0; i < nVC; i++) {
    Rcpp::NumericMatrix tmpcv = vc[i];
    VC[i] = Rcpp::as<arma::mat>(tmpcv);
    // VC[i](tmpcv.begin(), n, n, false);
  }
  VC[nVC].set_size(n,n);
  VC[nVC].eye();
  
  

  // Initial estimate for beta (from independence model)
  arma::colvec beta = arma::solve(X, Y);
  
  arma::colvec resid = Y - X*beta;
  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  theta(nVC) = sig2;

  arma::mat Omega = sig2*arma::eye<arma::mat>(n,n);
  arma::mat IOmega, P, xOmegax, IOmegaX, IOmega2;
  arma::colvec mu = X * beta;
  arma::mat Fisher(nVC+1, nVC+1), InvFisher(nVC+1, nVC+1);
  arma::colvec Deriv(nVC+1);
  arma::colvec WorkingTheta(k);

  int i, j;

  // Main iteration
  for (int iter=0; iter<maxiter; iter++) {

    printf("Iter %d\n", iter);

    // Compute the Inverse variance matrix
    Omega.zeros();
    
    for (i=0; i<=nVC; i++) {
      Omega += theta(i)*VC[i];  
    }

    IOmega = inv_sympd(Omega);

    
    // X^t Omega X
    P = IOmega;
    IOmegaX = (IOmega*X);

    xOmegax = arma::inv(arma::trans(X) * (IOmegaX));
    P = IOmega - (IOmegaX*xOmegax)*arma::trans(IOmegaX);   

    IOmega2 = P*P;
    
    mu = X * beta;

    Deriv.zeros();
    Fisher.zeros();
    
    for (i = 0; i < nVC+1; i++) {
      Deriv(i) += - arma::trace(P*VC[i]) + arma::as_scalar((arma::trans(Y)*P)*VC[i]*(P*(Y)));
      
      for (j = i; j < nVC+1; j++) {
	Fisher(i,j) += arma::trace(IOmega2*VC[i]*VC[j]);
	Fisher(j,i)  = Fisher(i,j);
      }
    }

    // Remember to scale the 1st and 2nd derivatives by 1/2
    Deriv  *= .5;
    Fisher *= .5;
    /*
    arma::colvec OldBeta = beta;
*/
    /*
    if (Method == METHOD_ML) {
      beta += Inverse(xOmegax)*DeltaBeta;
      DeltaBeta = Inverse(xOmegax)*DeltaBeta;
    }
    */


        
    // Inverts the matrix
    arma::mat IFisher = inv(Fisher);

    // Calculate the change in delta
    arma::colvec DeltaTheta = IFisher*Deriv;
    WorkingTheta = DeltaTheta;

    /*
      dLogDet = 0;
      for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	IOmega[nFam] = 0;
	for (i = 0; i < nVC; i++) {	
	  IOmega[nFam] += VC[nFam+nPedigrees*i]*theta(i+1,1);
	}
	IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	dLogDet += dLogDet2;
      }
      // Calculate P
      InvOmega2 = P*P;

      //	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
      //	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));

      for (i = 0; i < nVC; i++) {
	Deriv(i+1,1) += - Trace(P*NewVC[i]) 
	  + ((Transpose(NewY)*P)*NewVC[i]*(P*(NewY)))(1,1);
	  
	for (j = i; j < nVC; j++) {
	  Fisher(i+1,j+1) += Trace(InvOmega2*NewVC[i]*NewVC[j]);
	  Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	}
      }

      LogLike = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
      LogLike -= - 0.5*nMeanParam*log(2*PI);
    */

    
    
  }

  
  
    // create a new data frame and return it
  return DataFrame::create(Rcpp::Named("coefficients")=beta,
			   Rcpp::Named("theta")=theta,
			   Rcpp::Named("Omega")=Omega,
			   Rcpp::Named("newtheta")=WorkingTheta
			   );
  
}



/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * This file is part of the statistical program Pedipet.
 * Functions to do ML and REML of mixed model
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "maximize.h"

// #define debug
#define constrained       // Do constrained maximization (i.e. sigma >= 0)
#define TOLERANCE         0.00001
#define CLOSETOEDGE       0.001       // When is a parameter close to the edge
#define PCLOSETOEDGE      0.001       // When is the mixing parameter close to the edge
#define SCALEEPSILON      0.1
#define MINIMUMSTEPSIZE   -30         // På log-skala
#ifndef PI                // Definition of PI when lacking in compiler
#define PI                3.141592654 
#endif

#ifndef SQR
static double sqrarg;
#define  SQR(a)                ((sqrarg=(a))==0.0 ? 0.0 : sqrarg*sqrarg)
#endif

#define REML

#ifndef min
#define min(a,b)          ((a) < (b) ? (a) : (b))
#define max(a,b)          ((a) > (b) ? (a) : (b))
#endif

#define constrainmu       .05
#define output

//#define newthings

//#define improveconv
#define fix_error
#define MIXTURECONVLAG    8
//#define REALCALC          // If a parameter is close to 0 than CLOSETOEDGE, then set it to 0 for likelihood calculations
#define GEM               // Use general EM
#define obsdatalikelihood
#define diffbeta
#define edgewalk          // If true then walk along edges without trying to get close to 0

/*
char *ConvergenceText[5] = {
  "No convergence",
  "No change in likelihoods",
  "No change in mean parameters",
  "Stepsize too small",
  "Max iterations reached",
};
 

/*

  MaximizeMixedModel maximizes a MM under the constraints, that the VC should be non-negative.
  No missing data should be present in the input data

  Input    : nPedigrees    number of pedigrees
             nPedSize      array describing the size of the nPedigrees pedigrees
             y             array of matrices (each of size nPedSize[i]*1) of responses
             x             array of matrices (each of size nPedSize[i]*p) of design matrices
	     nVC           number of variance components
             VC            array of variance components matrices. 
                           Indexed by families (the first nPedigrees matrices are alle VC_1, 
                           the next nPedigrees are VC_2 etc.)
	     start         a vector of starting values for the variance components
             beta          a vector of starting values for the mean parameters
	     Method        Used  maximization method (see below). 1 = ML
             Constrain     Should we use constrained optimization? 0 = no, otherwise yes
                           Isn't implemented. Allways constrained VC to be non-negative.
         


  Returns the negative log likelihood value at the minimum. 
  Start holds the estimated variance components paraemters

  REML doesn't work, as you can't work on families but have to work with the complete
  variance matrix (P isn't block-diagonal)
             
  (C) Claus Ekstrøm 1999--2000

*/


/*
double MaximizeMixedModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, int Method, int Constrain, int PrintInfo)
{
  int i, j, nIter, nFam, nCode, nTotal=0;
  int nMeanParam, convergence;
  double dLogLike, dLogDet, dLogDet2, LogLike, downscale, change;
  double step, newll;

  convergence = 0;

  // Generel matrices
  MATRIX xT;
  MATRIX Omega, InvOmega, P, InvOmega2;
  MATRIX Fisher(nVC, nVC), InvFisher(nVC, nVC);
  MATRIX Deriv(nVC, 1);
  MATRIX theta(nVC, 1), DeltaTheta(nVC, 1), WorkingTheta(nVC, 1);;  
  MATRIX DeltaBeta = beta;
  MATRIX OldBeta = beta;
  MATRIX xOmegax(beta.Rows(), beta.Rows());

  // Matrices for REML estimation
  MATRIX *IOmega = NULL;
  MATRIX NewY, NewX, *NewVC = NULL;


  nMeanParam = x[0].Cols();

#ifdef REML
  if (Method != METHOD_ML && Method != METHOD_REML ) {
    printf("ERROR: Trying to maximize likelihood using non-(RE)ML methods.\nThis is not implemented yet.\n");
    exit(1);
  }
#else
  if (Method != METHOD_ML) {
    printf("ERROR: Trying to maximize likelihood using non-ML methods.\nThis is not implemented yet.\n");
    exit(1);
  }
#endif

  // Should combine all the matrices to one large set
  if (Method == METHOD_REML) {

#ifdef debug
  printf("Preparing REML code\n");
#endif

    NewVC = new MATRIX[nVC];

    nTotal=0;
    for (i = 0; i< nPedigrees; i++) {
      nTotal += nPedSize[i];
    }

    for (i = 0; i< nVC; i++) {
      NewVC[i].Resize(nTotal,nTotal);
      NewVC[i] = CombineMatrices(&VC[i*nPedigrees], nPedigrees);
    }

    NewY.Resize(nTotal,1);
    NewX.Resize(nTotal,beta.Rows());

    NewY = AppendMatrices(y, nPedigrees);
    NewX = AppendMatrices(x, nPedigrees);

    xT.Resize(beta.Rows(), nTotal);
    xT = Transpose(NewX);

    MATRIX tempmat(nTotal, nTotal);
    tempmat = 0;
    for (i = 1; i<=nTotal; i++) {
      tempmat(i,i) = 1;
    }
    NewY = (tempmat - NewX*Inverse(xT*NewX)*xT)*NewY;

    // Holds the pedigree variances below
    IOmega = new MATRIX[nPedigrees];
    for (nFam = 0; nFam < nPedigrees; nFam++){      
      IOmega[nFam].Resize(nPedSize[nFam], nPedSize[nFam]);
    }
    P.Resize(nTotal, nTotal);
    InvOmega.Resize(nTotal, nTotal);
    InvOmega2.Resize(nTotal, nTotal);


  }

#ifdef debug
  printf("Start maximizing\n");
#endif

  // Start checking that everythink looks ok

  theta = start;  

  // Starts iterating
  for (nIter = 1; nIter <= MIXEDMODELMAXITER; nIter++)  {

#ifdef debug
    printf("Iteration: %3d\n", nIter);
#endif
    Fisher = 0.0;
    Deriv  = 0.0;
    DeltaBeta = 0.0;

    dLogLike = 0.0;
    LogLike = 0.0;

    xOmegax = 0.0;
    

    if (Method == METHOD_ML) {
      // This part calculates the first derivative and the Fisher information
      // matrix using Maximum likelihood methods

#ifdef debug
    theta.Print();
#endif


      // Go through each family/pedigree
      // The block diagonal structure is kept when doing ML
      for (nFam = 0; nFam < nPedigrees; nFam++)	{
	// Calculate the current variance for the family
	Omega.Resize(nPedSize[nFam], nPedSize[nFam]);
	InvOmega.Resize(nPedSize[nFam], nPedSize[nFam]);
	P.Resize(nPedSize[nFam], nPedSize[nFam]);
	xT.Resize(x[nFam].Cols(), x[nFam].Rows());
	xT = Transpose(x[nFam]);
	
	Omega = 0.0;

	nCode = nFam*nPedigrees;

	for (i = 0; i < nVC; i++)
	  Omega += VC[nFam+nPedigrees*i]*theta(i+1,1);

	// Inverting Omega_i
	InvOmega = Inverse(Omega, &dLogDet);

	// Add to the log likelihood
	MATRIX mu = (x[nFam] * beta);
	
	xOmegax   += (xT*InvOmega)*x[nFam];
	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));
     
	InvOmega2.Resize(nPedSize[nFam], nPedSize[nFam]);
	MATRIX UsedMatrix(nPedSize[nFam], nPedSize[nFam]);
	UsedMatrix = InvOmega;

	InvOmega2 = UsedMatrix*UsedMatrix;
	// Calculating the derivatives and Fisher scoring matrix
	for (i = 0; i < nVC; i++) {
	  Deriv(i+1,1) += - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
	    + ((Transpose(y[nFam] - mu)*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-mu)))(1,1);

	  for (j = i; j < nVC; j++) {
	    Fisher(i+1,j+1) += Trace(InvOmega2*VC[nFam + nPedigrees*i]*VC[nFam + nPedigrees*j]);
	    Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	  }
	}
      }
    }
    else if (Method == METHOD_REML) {
      // This part calculates the first derivative and the Fisher information
      // matrix using Restricted maximum likelihood methods

      // Calculate the current inverse variance matrix
      dLogDet = 0;
      for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	IOmega[nFam] = 0;
	for (i = 0; i < nVC; i++) {	
	  IOmega[nFam] += VC[nFam+nPedigrees*i]*theta(i+1,1);
	}
	IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	dLogDet += dLogDet2;
      }
      InvOmega = CombineMatrices(IOmega, nPedigrees);

      // Calculate P
      xOmegax = Inverse((xT*InvOmega)*NewX, &dLogDet2);
      P = InvOmega - ((InvOmega*NewX)*xOmegax)*(xT*InvOmega);   
      InvOmega2 = P*P;

      MATRIX mu = (NewX * beta);

      //	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
      //	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));

      for (i = 0; i < nVC; i++) {
	Deriv(i+1,1) += - Trace(P*NewVC[i]) 
	  + ((Transpose(NewY)*P)*NewVC[i]*(P*(NewY)))(1,1);
	  
	for (j = i; j < nVC; j++) {
	  Fisher(i+1,j+1) += Trace(InvOmega2*NewVC[i]*NewVC[j]);
	  Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	}
      }

      LogLike = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
      LogLike -= - 0.5*nMeanParam*log(2*PI);
    }


    // Remember to scale the 1st and 2nd derivatives by 1/2
    Deriv  *= .5;
    Fisher *= .5;

    OldBeta = beta;

    if (Method == METHOD_ML) {
      beta += Inverse(xOmegax)*DeltaBeta;
      DeltaBeta = Inverse(xOmegax)*DeltaBeta;
    }

#ifdef constrained
    // Modifies the diagonal of the Fisher Matrix
    for (i=1; i<=nVC; i++)
      Fisher(i,i) += constrainmu/theta(i,1);

#endif

    // Inverts the matrix
    InvFisher = Inverse(Fisher);

    // Calculates the change in delta
    DeltaTheta = InvFisher*Deriv;
    WorkingTheta = DeltaTheta;

    // This change should be constrained
#ifdef constrained

    downscale = 1.0;  // Maximum steplength to stay positive

    // Walk along edge for parameters close to the edge
    for (i=1; i<=nVC; i++) {
      if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	WorkingTheta(i,1)=0;
      // Now finds the maximum allowed steplength 
      if (WorkingTheta(i,1)<0)
	downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));
    }
    // Scales the step down
    DeltaTheta = DeltaTheta*downscale;
    WorkingTheta = WorkingTheta*downscale;
#endif

    step = 2.0;
    do  {
      step *= .5;
      WorkingTheta = DeltaTheta*step;

      for (i=1; i<=nVC; i++) {
	if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	  WorkingTheta(i,1)=0;
      }

      newll = 0.0;
      // Calculate the new Omega

      if (Method == METHOD_ML) {
	for (nFam = 0; nFam < nPedigrees; nFam++)	{
	  // Calculate the current variance for the family
	  MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX xT = Transpose(x[nFam]);
	  
	  Omega = 0.0;
	  
	  nCode = nFam*nPedigrees;
	  
	  for (i = 0; i < nVC; i++)
	    Omega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + WorkingTheta(i+1,1));
	  
	  // Inverting Omega_i
	  InvOmega = Inverse(Omega, &dLogDet);
	  
	  // Add to the log likelihood
#ifdef fix_error
	  MATRIX mu = (x[nFam] * beta);
#else
	  MATRIX mu = (x[nFam] * OldBeta);
#endif
	  newll += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));
	}
	  
      } else if (Method == METHOD_REML) {
	// Calculate the current inverse variance matrix
	dLogDet = 0;
	for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	  IOmega[nFam] = 0;
	  for (i = 0; i < nVC; i++) {	
	    IOmega[nFam] += VC[nFam+nPedigrees*i]*(theta(i+1,1) + WorkingTheta(i+1,1));
	  }
	  IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	  dLogDet += dLogDet2;
	}
	InvOmega = CombineMatrices(IOmega, nPedigrees);
	xOmegax = Inverse((xT*InvOmega)*NewX, &dLogDet2);
	P = InvOmega - ((InvOmega*NewX)*xOmegax)*(xT*InvOmega);   

	newll = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
    	newll -= -0.5*nMeanParam*log(2*PI);	



      }    
    }
    while (newll < LogLike && fabs(newll-LogLike)>TOLERANCE*downscale*step);

    // Modifies Theta
    theta += WorkingTheta;

    // XXX
    // Should also check the change in beta for ML 
    // Exiting due to small increas in LL or theta
    change = 0.0;

    for (i=1; i<= nVC; i++)
      change += fabs(WorkingTheta(i,1));

    // if ML and no change in mean parameters
    if (Method == METHOD_ML) {
      for (i=1; i<= beta.Rows(); i++)
	change += fabs(DeltaBeta(i,1));
    }


    // Should perhaps multiply by step also?
    // If yes, then same criteria should be used above in the while loop
    if (fabs(newll-LogLike)<TOLERANCE*downscale*step || change<TOLERANCE*downscale*step) {
      if (fabs(newll-LogLike)<TOLERANCE*downscale*step) {
	convergence = 1;
      }
      else {
	convergence = 2;
      }

      // Calculating the mean parameter estimates based on the variances
      if (Method == METHOD_REML) {
	NewY = AppendMatrices(y, nPedigrees);
	beta = xOmegax*(xT*InvOmega)*NewY;

	// Fixes xOmegax
	// Above, when using REML xOmegax is the INVERSE of xomega x, but for
	// ML is is xomegax
	// Therefore, change it back sp the right variance of the
	// mean is printed

	xOmegax = Inverse(xOmegax);
	
      }

      RES1 = theta(1,1);
      if (theta.Rows()>=2)       
	RES2 = theta(2,1);
      if (theta.Rows()>=3) 
	RES3 = theta(3,1);
      break;
    }
  }
#ifdef output
  if (PrintInfo != 0) {
    printf("Convergence reached after %d iterations (%s)\n",nIter, ConvergenceText[convergence]);
    printf("Estimating using %s.\n", Method==METHOD_ML ? "ML" : "REML");
    printf("Estimated variance components:\n");
    theta.Print();
    printf("Estimated mean parameters:\n");
    beta.Print();
    printf("Log likelihood at maximum: %6.3f\n", newll);
#ifdef showvariances	
    printf("Inverse Fisher information matrix for VC:\n");
    InvFisher.Print();
    printf("Variances of the mean parameters:\n");
    Inverse(xOmegax).Print();
#endif
  }
    
#endif

  fflush(stdout);

  if (Method == METHOD_REML) {
    delete[] IOmega;
    delete[] NewVC;
  }


  return (newll);
}


*/

