#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <thread>

#define NTHREADS 16

/*
 * example2_mex.c
 *
 *EXAMPLE2 Compute log odds for Aspen Yoo's word recognition memory model.
 *
 * ================ INPUT VARIABLES ====================
 * M: number of features. [scalar] (integer)
 * SIGMA: memory noise. [scalar] (double)
 * NS: number of samples of SNew and SOld. [scalar] (integer)
 * NNEW: number of new words. [scalar] (integer)
 * NOLD: number of old words. [scalar] (integer)
 * SNEW: new words across S simulations. [Nnew*nS,M] (double)
 * X: noisy memories. [Nold,M] (double)
 *
 * ================ OUTPUT VARIABLES ==================
 * D_NEW: log odds of new trials. [Nnew*nS,1] (double)
 *
 * This is a MEX-file for MATLAB.
 */

void example2(double* St, double* Xt, double* d, const int k, const int Srows, const int Nold, const int M, const double J){
  double SUM,sum,tmp;
  for(int l=k;l<Srows;l+=NTHREADS){
    d[l]=0.5*M*log(J);
    for(int j=l*M;j<l*M+M;j++)
      d[l] += 0.5*St[j]*St[j];
    SUM=0.0;
    for (int i=0; i < Nold; i++) {
      sum = 0.0;
      for (int j=l*M, r=i*M; j < l*M+M; j++,r++)
        sum -= (St[j] - Xt[r]) * (St[j] - Xt[r]);
      SUM += exp(0.5 * J * sum);
    }
    d[l] += log(SUM)-log(Nold);
  }
}

void example2_threads( double *d, int M, double sigma, int nS, int Nold, int N, double *S, int Srows, double* X){
    std::thread t[NTHREADS];
    double J=1.0/(sigma*sigma) + 1.0;
    double Jss=(1.0/J)/(sigma*sigma);
    double *St=new double[Srows*nS];
    double *Xt=new double[Srows*Nold];
    for(int i=0;i<Srows;i++)
      for(int j=0;j<nS;j++)
        St[i*nS+j]=S[j*Srows+i];
    for(int i=0;i<Nold;i++)
      for(int j=0;j<M;j++)
        Xt[i*M+j]=Jss*X[j*Nold+i];
    for(int k=0;k<NTHREADS;k++)
      t[k]=std::thread(example2,S,X,d,k,Srows,Nold,M,J);
    for(int k=0;k<NTHREADS;k++)
      t[k].join();
    delete[] St;
    delete[] Xt;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *d_new, sigma, *SNew, *X;
	int M, nS, Nnew, Nold;

	/* Get first input (M, scalar int) */
	M = (int) mxGetScalar(prhs[0]);

	/* Get second input (SIGMA, scalar double) */
	sigma = (double) mxGetScalar(prhs[1]);

	/* Get third input (NS, scalar int) */
	nS = (int) mxGetScalar(prhs[2]);

	/* Get fourth input (NNEW, scalar int) */
	Nnew = (int) mxGetScalar(prhs[3]);

	/* Get fifth input (NOLD, scalar int) */
	Nold = (int) mxGetScalar(prhs[4]);

	/* Get sixth input (SNEW, Nnew*nS-by-M double) */
	SNew = (double*) mxGetPr(prhs[5]);

	/* Get seventh input (X, Nold-by-M double) */
	X = (double*) mxGetPr(prhs[6]);

	/* Pointer to first output (D_NEW, Nnew*nS-by-1 double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (Nnew*nS), (mwSize) (1), mxREAL);
	d_new = mxGetPr(plhs[0]);

	/* Call the C subroutine */
  example2_threads(d_new, M, sigma, nS, Nnew, Nold, SNew, Nnew*nS, X);

}
