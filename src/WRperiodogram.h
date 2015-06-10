/*
  Whittaker and Robinson periodogram; Whittaker and Robinson (1924),
  Legendre & Legendre (1998, 2012, Section 12.4.1).
  x : a vector of quantitative data
  T1: first period to analyse
  T2: last period to analyse (T2 <= n/2)
  nperm: number of permutations for tests of significance
  mult : methods for correction for multiple testing; "sidak" or "bonferroni"

  License: GPL-2 
  Authors:: Pierre Legendre, September 2012, Guillaume Guenard, March 2014
  C header
*/

#ifndef WRPERIODOGRAM_H
#define WRPERIODOGRAM_H
#include<R.h>

void BBCMVAR(double* x, int* nx, int* T1, int* T2, double* out, double* cmacc, int* cmden);
/* array out should be pre-allocated to size T2-T1+1, i.e., the number of periods calculated.
   array cm and cmden should be pre-allocated to size T2, i.e. the maximum period.
   they were not internally allocated to spare system calls (and therefore save time) */

void WRperiodogram(double* x, int* nx, int* T1, int* T2, double* out, int* nperm, int* pidx, int* npidx, int* permout);
/* array out and permout should be pre-allocated to size T2-T1+1 when calling from R
*/

#endif  // WRPERIODOGRAM_H
