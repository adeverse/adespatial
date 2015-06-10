/*
  Whittaker and Robinson periodogram; Whittaker and Robinson (1924),
  Legendre & Legendre (1998, 2012, Section 12.4.1).
  x : a vector of quantitative data
  T1: first period to analyse
  T2: last period to analyse (T2 <= n/2)
  nperm: number of permutations for tests of significance
  mult : methods for correction for multiple testing; "sidak" or "bonferroni"

  License: GPL-2 
  Authors:: Pierre Legendre, September 2012; Guillaume Guenard, March 2014,
  C function
*/

#include"WRperiodogram.h"

void BBCMVAR(double* x, int* nx, int* T1, int* T2, double* out, double* cmacc, int* cmden)
{
  int i, j, k;                  // i: index of period, j: index of column, k: index of column elements
  int gmden, varden;            // Denominators of the great column means and column variance
  double gmacc, varacc, subuf;  // Accumulators for the great column means and column variance, and subtraction buffer
  for(i = 0; i <= *T2 - *T1; i++)    // Period offset starting from T1
    {
      varacc = gmacc = 0.0;        // Resetting accumulators
      varden = gmden = 0;          // Resetting denominators
      for(j = 0; j < *T1 + i; j++) // Resurse columns
	{
	  cmacc[j] = 0.0;          // Resetting column means accumulators
	  cmden[j] = 0;            // Resetting column means denominators
	  for(k = j; k < *nx; k += *T1 + i)  // Apply offset T1 + i to recurse through rows
	    {
	      if(!ISNA(x[k]))          // Accumulate only when not NA.
		{
		  cmacc[j] += x[k];
		  cmden[j]++;
		}
	    }
	  if(cmden[j])                 // If anything was accumulated:
	    {
	      cmacc[j] /= cmden[j];    // Compute column j mean value
	      gmacc += cmacc[j];       // Accumulate the value to the great column mean
	      gmden++;
	    }
	}
      if(gmden)                        // If anything was accumulated for the great mean:
	{
	  gmacc /= gmden;
	  for(j = 0; j < *T1 + i; j++)  // Calculate variance (mean square deviation of the
	    {                           // individual column means about the great column mean.
	      if(cmden[j]) {
		subuf = cmacc[j] - gmacc;
		varacc += subuf * subuf;
		varden++;
	      }
	    }
	  out[i] = varacc / varden;
	}
      else
	out[i] = NA_REAL;               // If nothing was accumulated no variance calculation
    }
  return;
}

void WRperiodogram(double* x, int* nx, int* T1, int* T2, double* out, int* nperm, int* pidx, int* npidx, int* permout)
{
  // Allocate memory for the calculation of column means
  double rnb;
  int i, j;
  double* cmacc = (double*)Calloc(*T2, double);  // Accumulator
  int* cmden = (int*)Calloc(*T2, int);           // Denominator
  if(cmacc == NULL || cmden == NULL)
    error("Dynamic memory allocation failure in C function BBCMVAR");
  BBCMVAR(x, nx, T1, T2, out, cmacc, cmden);
  if(*nperm > 0)
    {
      double buffer;
      int idx;
      double* outperm = (double*)Calloc(*T2 - *T1 + 1, double);
      if(permout == NULL)
	error("Dynamic memory allocation failure in C function BBCMVAR");
      // Perform permulation test.
      GetRNGstate();
      // The permutation loop.
      for(i = 0; i < *nperm; i++)
	{
	  // I sould add some break point.
	  for(j = 0; j < *npidx; j++)
	    {
	      do
		{
		  rnb = unif_rand();
		} while (rnb == 1.0);
	      idx = (int)(rnb * *npidx);
	      buffer = x[pidx[idx]];
	      x[pidx[idx]] = x[pidx[j]];
	      x[pidx[j]] = buffer;
	    }
	  BBCMVAR(x, nx, T1, T2, outperm, cmacc, cmden);
	  for(j = 0; j < *T2 - *T1 + 1; j++)
	    {
	      if(outperm[j] >= out[j])
		permout[j]++;
	    }
	}
      // End of the permutation loop
      PutRNGstate();
      Free(outperm);
    }
  Free(cmden);
  Free(cmacc);
  return;
}
