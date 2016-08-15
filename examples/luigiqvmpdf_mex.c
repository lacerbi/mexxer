#include "mex.h"
#include "math.h"
#include "matrix.h"

/*
 * luigiqvmpdf_mex.c
 *
 * calculate von mises probability 
 * (with many bsxfuns so you don't have to repmat input if you want to 
 * calculate many probabilities for many different distributions)
 *
 * INPUT ===
 *   x : the angle(s) for the pdf to be evaluated [M-by-N] (double)
 *   mu: the mean(s) of the Von Mises distribution [M-by-N] (double)
 *   k : the concentration(s) of the Von Mises distribution [M-by-1] (double)
 *   B : besseli(0, k, 1) [M-by-1] (double)
 *
 * OUTPUT ===
 *   p : the probability of the given angle(s) [M-by-N] (double)
 *
 * Thanks for giving me this a year ago Luigi!
 *
 * Also maybe this is just slow because the thing to do was not to stick
 * bsxfuns everywhere? I'm don't really know anything about how bsxfun works
 * 
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 13-Aug-2016 with MEXXER v0.1 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void luigiqvmpdf( double *p, double *x, double *mu, double *k, double *B, mwSize M, mwSize N )
{
	mwSize i,j;
    double C1,C2,*k0,*B0,*tmp,*tmp0;
    C1 = 3.14159265358979323846264338327950288 / 180.0;
    C2 = log(360.0);
    
    tmp = malloc(M * sizeof(double) );
        
    k0 = k;
    B0 = B;
    tmp0 = tmp;
    
    for (i=0;i<M;i++) {
        *(tmp++) = -(C2 + log(*(B++)) + *(k++));
    }    
    
    for (j=0;j<N;j++) {
        k = k0;
        tmp = tmp0;
        for (i=0;i<M;i++,k++,tmp++) {
            *(p++) = exp(*k * cos( C1 * (*(x++) - *(mu++)) ) + *tmp);
        }
    }
    
    /* p = exp( bsxfun(@plus, bsxfun(@times, k, cos((pi/180).*bsxfun(@minus, x,mu))),...
    -1*(log(360)+log(B) + k)  )); */

    /* free(tmp); */
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *p, *x, *mu, *k, *B;
	mwSize *dims_x, *dims_mu, *dims_k, *dims_B;
	mwSize M, N;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<4 || nrhs>4 )
		mexErrMsgIdAndTxt( "MATLAB:luigiqvmpdf:invalidNumInputs",
			"Four inputs required.");
	if ( nlhs<1 || nlhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:luigiqvmpdf:invalidNumOutputs",
			"One outputs required.");

	/* Get first input (X, M-by-N double) */
	x = (double*) mxGetPr(prhs[0]);
	dims_x = (mwSize*) mxGetDimensions(prhs[0]);
	M = dims_x[0];
	N = dims_x[1];

	/* Get second input (MU, M-by-N double) */
	mu = (double*) mxGetPr(prhs[1]);
	dims_mu = (mwSize*) mxGetDimensions(prhs[1]);

	/* Get third input (K, M-by-1 double) */
	k = (double*) mxGetPr(prhs[2]);
	dims_k = (mwSize*) mxGetDimensions(prhs[2]);

	/* Get fourth input (B, M-by-1 double) */
	B = (double*) mxGetPr(prhs[3]);
	dims_B = (mwSize*) mxGetDimensions(prhs[3]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
	if ( ARGSCHECK ) {
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:xNotReal", "Input X must be real.");
		if ( dims_x[0] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:xWrongSize", "The first dimension of input X has the wrong size (should be M).");
		if ( dims_x[1] != ((mwSize) (N)) )
			mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:xWrongSize", "The second dimension of input X has the wrong size (should be N).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:muNotReal", "Input MU must be real.");
		if ( dims_mu[0] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:muWrongSize", "The first dimension of input MU has the wrong size (should be M).");
		if ( dims_mu[1] != ((mwSize) (N)) )
			mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:muWrongSize", "The second dimension of input MU has the wrong size (should be N).");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:kNotReal", "Input K must be real.");
		if ( dims_k[0] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:kWrongSize", "The first dimension of input K has the wrong size (should be M).");
		if ( dims_k[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:kWrongSize", "The second dimension of input K has the wrong size (should be 1).");

		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
				mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:BNotReal", "Input B must be real.");
		if ( dims_B[0] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:BWrongSize", "The first dimension of input B has the wrong size (should be M).");
		if ( dims_B[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:luigiqvmpdf:BWrongSize", "The second dimension of input B has the wrong size (should be 1).");
	}

	/* Pointer to first output (P, M-by-N double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (M), (mwSize) (N), mxREAL);
	p = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	luigiqvmpdf(p, x, mu, k, B, M, N);

}
