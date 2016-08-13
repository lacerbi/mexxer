#include "mex.h"
#include "math.h"
#include "matrix.h"

/*
 * testmex_mex.c
 *
 *% calculate log odds
 * function used as a basis for Luigi to code up C code!
 *
 * ================ INPUT VARIABLES ====================
 * SIGMA: memory noise. [scalar] (double)
 * NS: number of samples of SNew and SOld. [scalar] (integer)
 * NNEW: number of new words. [scalar] (integer)
 * NOLD: number of old words. [scalar] (integer)
 * SNEW: new words across S simulations. [Nnew*nS,1,M] (double)
 * SOLD: old words across S simulations. [Nold*nS,TEMP] (double)
 * X: noisy memories. [Nold,M] (double)
 * 
 * ================ OUTPUT VARIABLES ==================
 * D_NEW: log odds of new trials. [Nnew*nS,1,nS] (double)
 * D_OLD: log odds of old trials. [Nold*nS,1] (double)
 *
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 13-Aug-2016 with MEXXER v0.1 
 * (https://github.com/lacerbi/mexxer).
*/

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void testmex( double *d_new, double *d_old, double sigma, int nS, int Nnew, int Nold, double *SNew, double *SOld, double *X, mwSize M, mwSize TEMP )
{
	
	/* Write your main calculations here... */
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *d_new, *d_old, sigma, *SNew, *SOld, *X;
	int nS, Nnew, Nold;
	mwSize *dims_SNew, *dims_SOld, *dims_X, dims_d_new[3];
	mwSize M, TEMP;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<7 || nrhs>7 )
		mexErrMsgIdAndTxt( "MATLAB:testmex:invalidNumInputs",
			"Seven inputs required.");
	if ( nlhs<2 || nlhs>2 )
		mexErrMsgIdAndTxt( "MATLAB:testmex:invalidNumOutputs",
			"Two outputs required.");

	/* Get first input (SIGMA, scalar double) */
	sigma = (double) mxGetScalar(prhs[0]);

	/* Get second input (NS, scalar int) */
	nS = (int) mxGetScalar(prhs[1]);

	/* Get third input (NNEW, scalar int) */
	Nnew = (int) mxGetScalar(prhs[2]);

	/* Get fourth input (NOLD, scalar int) */
	Nold = (int) mxGetScalar(prhs[3]);

	/* Get fifth input (SNEW, Nnew*nS-by-1-by-M double) */
	SNew = (double*) mxGetPr(prhs[4]);
	dims_SNew = (mwSize*) mxGetDimensions(prhs[4]);
	M = dims_SNew[2];

	/* Get sixth input (SOLD, Nold*nS-by-TEMP double) */
	SOld = (double*) mxGetPr(prhs[5]);
	dims_SOld = (mwSize*) mxGetDimensions(prhs[5]);
	TEMP = dims_SOld[1];

	/* Get seventh input (X, Nold-by-M double) */
	X = (double*) mxGetPr(prhs[6]);
	dims_X = (mwSize*) mxGetDimensions(prhs[6]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
	if ( ARGSCHECK ) {
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetN(prhs[0])*mxGetM(prhs[0])!=1) )
			mexErrMsgIdAndTxt("MATLAB:testmex:sigmaNotScalar", "Input SIGMA must be a scalar.");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetN(prhs[1])*mxGetM(prhs[1])!=1) )
			mexErrMsgIdAndTxt("MATLAB:testmex:nSNotScalar", "Input NS must be a scalar.");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetN(prhs[2])*mxGetM(prhs[2])!=1) )
			mexErrMsgIdAndTxt("MATLAB:testmex:NnewNotScalar", "Input NNEW must be a scalar.");

		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || (mxGetN(prhs[3])*mxGetM(prhs[3])!=1) )
			mexErrMsgIdAndTxt("MATLAB:testmex:NoldNotScalar", "Input NOLD must be a scalar.");

		if ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) )
				mexErrMsgIdAndTxt("MATLAB:testmex:SNewNotReal", "Input SNEW must be real.");
		if ( dims_SNew[0] != ((mwSize) (Nnew*nS)) )
			mexErrMsgIdAndTxt("MATLAB:testmex:SNewWrongSize", "The first dimension of input SNEW has the wrong size (should be Nnew*nS).");
		if ( dims_SNew[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:testmex:SNewWrongSize", "The second dimension of input SNEW has the wrong size (should be 1).");
		if ( dims_SNew[2] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:testmex:SNewWrongSize", "The third dimension of input SNEW has the wrong size (should be M).");

		if ( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) )
				mexErrMsgIdAndTxt("MATLAB:testmex:SOldNotReal", "Input SOLD must be real.");
		if ( dims_SOld[0] != ((mwSize) (Nold*nS)) )
			mexErrMsgIdAndTxt("MATLAB:testmex:SOldWrongSize", "The first dimension of input SOLD has the wrong size (should be Nold*nS).");
		if ( dims_SOld[1] != ((mwSize) (TEMP)) )
			mexErrMsgIdAndTxt("MATLAB:testmex:SOldWrongSize", "The second dimension of input SOLD has the wrong size (should be TEMP).");

		if ( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) )
				mexErrMsgIdAndTxt("MATLAB:testmex:XNotReal", "Input X must be real.");
		if ( dims_X[0] != ((mwSize) (Nold)) )
			mexErrMsgIdAndTxt("MATLAB:testmex:XWrongSize", "The first dimension of input X has the wrong size (should be Nold).");
		if ( dims_X[1] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:testmex:XWrongSize", "The second dimension of input X has the wrong size (should be M).");
	}

	/* Pointer to first output (D_NEW, Nnew*nS-by-1-by-nS double) */
	dims_d_new[0] = (mwSize) (Nnew*nS);
	dims_d_new[1] = (mwSize) (1);
	dims_d_new[2] = (mwSize) (nS);
	plhs[0] = mxCreateNumericArray(3, dims_d_new, mxDOUBLE_CLASS, mxREAL);
	d_new = mxGetPr(plhs[0]);

	/* Pointer to second output (D_OLD, Nold*nS-by-1 double) */
	plhs[1] = mxCreateDoubleMatrix((mwSize) (Nold*nS), (mwSize) (1), mxREAL);
	d_old = mxGetPr(plhs[1]);

	/* Call the C subroutine */
	testmex(d_new, d_old, sigma, nS, Nnew, Nold, SNew, SOld, X, M, TEMP);

}
