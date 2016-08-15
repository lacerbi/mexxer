#include "mex.h"
#include "math.h"

/*
 * VestBMS_likec1qtrapz.c
 *
 * multiplies two matrices elementwise and executes trapezoidal 
 * integration along first dimensions.
 *
 * This is a MEX-file for MATLAB.
 */

void calculate_d(double *d, int M, double sigma, int nS, int Nold, int N, double *S, int Srows, double* X)
{
    mwSize i,j,k;
    double J,Jss,tmp;
    double *S0,*X0,*d0,sum,SUM;
        
    /* store initial position */
    S0 = S;     
    d0 = d;
        
    /* Compute d = M/2*log(1+1/sigma^2) + 0.5*sum(SNew.^2,2) */

    J = 1/(sigma*sigma) + 1;
    
    /* this summation is awkward because MATLAB is column-major */
    for (i=0; i<Srows; i++) {
        sum = 0.;
        *d = 0.5 * M * log(J);
        S = S0 + i;   /* initialize pointer at current row */
        for (j=0; j<M; j++) {
            sum += (*S)*(*S);
            S += Srows; /* move pointer to next column */
        }
        *(d++) += 0.5 * sum;
    }
    
    
/* Compute d = d + log(squeeze(mean(exp( ...
 *        -0.5*J*sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) ...
 *                    - repmat(X/J/sigma^2,[1,1,Nnew*nS])).^2,2) ...
 *        )))) */
    
    Jss = (1/J)/(sigma*sigma);
    
    X0 = X;
    d = d0;
    
    for (k=0; k < Srows; k++) {
        S = S0 + k;
        SUM = 0.;
        
        for (i=0; i < Nold; i++) {
            X = X0 + i;
            
            sum = 0.;            
            for (j=0; j < M; j++) {
                tmp = S[j*Srows] - (X[j*Nold] * Jss);
                sum += tmp * tmp;
            }
            SUM += exp(-0.5 * J * sum);
        }    
        *(d++) += log( SUM / (double) Nold );
        
    }
    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /* const mwSize *dims_a,*dims_b; */
  double J;
  int M,nS,Nnew,Nold;
  double *SNew,*SOld,*X,*d_new,*d_old;
  double sigma;
  /* size_t n1,n2,n3; */
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=8)
    mexErrMsgIdAndTxt( "MATLAB:calculate_d:invalidNumInputs",
            "Eight inputs required.");
  if(nlhs!=2) 
    mexErrMsgIdAndTxt( "MATLAB:calculate_d:invalidNumOutputs",
            "Two outputs required.");
  
  /* check to make sure the first input argument is a scalar */
  /* if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
    mexErrMsgIdAndTxt( "MATLAB:calculate_d:xNotScalar",
            "Input M must be a scalar.");
  } */
  /* I don't do any check on the inputs for the sake of speed */
  
  /* Get first input (M, scalar integer) */
  M = (int) mxGetScalar(prhs[0]);
  
  /* Get second input (SIGMA, scalar double) */
  sigma = (double) mxGetScalar(prhs[1]);

  /* Get third input (NS, scalar integer) */
  nS = (int) mxGetScalar(prhs[2]);
  
  /* Get fourth input (NNEW, scalar integer) */
  Nnew = (int) mxGetScalar(prhs[3]);
  
  /* Get fifth input (NOLD, scalar integer) */
  Nold = (int) mxGetScalar(prhs[4]);
  
  /* Get pointer to sixth input (SNEW, Nnew*nS x M double) */
  SNew = mxGetPr(prhs[5]);
  
  /* Get pointer to seventh input (SOLD, Nold*nS x M double) */
  SOld = mxGetPr(prhs[6]);
  
  /* Get pointer to eight input (X, Nold x M double) */
  X = mxGetPr(prhs[7]);
  
  /* dims_a = mxGetDimensions(prhs[0]); */   /*  get the dimensions */
      
  /* Check matrix sizes */
  /* N = dims_a[0];
  K = dims_a[2];
  if ( dims_b[0] != N || dims_b[1] != K ) {
      mexErrMsgIdAndTxt( "MATLAB:VestBMS_likec1qtrapzc:dimensionMismatch",
            "Second matrix dimensions do not match the first input.");      
  }  
  /* printf("%d %d %d\n",dims_a[0],dims_a[1],dims_a[2]); */

  /*  output pointer to the first output matrix (D_NEW, Nnew*nS x 1 double) */
  plhs[0] = mxCreateDoubleMatrix((mwSize) (Nnew * nS), 1, mxREAL);
  
  /*  output pointer to the second output matrix (D_OLD, Nold*nS x 1 double) */
  plhs[1] = mxCreateDoubleMatrix((mwSize) (Nold * nS), 1, mxREAL);
  
  /*  create C pointers to copies of the output matrices */
  d_new = mxGetPr(plhs[0]);
  d_old = mxGetPr(plhs[1]);
    
  /*  call the C subroutine */
  calculate_d(d_new, M, sigma, nS, Nold, Nnew, SNew, Nnew*nS, X);
  calculate_d(d_old, M, sigma, nS, Nold, Nnew, SOld, Nold*nS, X);
  
}
