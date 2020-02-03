/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
|                                                            |
|  CONSTRAINED HIERARCHICAL CLUSTERING                       |
|  using (user-specified) criterion                          |
|                                                            |
|  C implementation of the Lance and Williams (1967) method  |
|  for hierarchical clustering with or without the spatial   |
|  contiguity constraint.                                    |
|                                                            |
|  Guillaume Guénard, Université de Montréal, Québec, Canada |
|  August 2018 - September 2019                              |
|                                                            |
|  The present implementation has been vastly inspirered by  |
|  the original HCLUST Fortran code by                       |
|  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       |
|  with modifications for R: Ross Ihaka, Dec 1996            |
|                            Fritz Leisch, Jun 2000          |
|                                                            |
\+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// C header

// Definitions

// #define testing      // Enables testing functions (useful during development)
// #define show_links
// Shows links before and after the clustering function to check link updating
#define true -1
#define false 0
#define INTMOD 16      // Interval to check for user interruption request

// Accessory functions (for Lance and Williams)
void setLWUpdate(unsigned int n, int method, double* diss0,
                 void (**update)(int,int*,int*,double*,double*,unsigned int,unsigned int));
void      getmin(unsigned int n, int* flag, double* diss0, unsigned int* nn_i,
                 unsigned int* nn_j, double* nn_dist);
void  getminlink(unsigned int n, double* diss0, unsigned int nl, int* linkl,
                 unsigned int* nn_i, unsigned int* nn_j, double* nn_dist);
void   mergelink(unsigned int nl, int* linkl, unsigned int i2, unsigned int j2);

/* 
 * A suite of update functions to recalculate the distances from the newly
 * aggregated clusters to the other entities (Lance and Williams)
 */

void     lw_Ward(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2);
void   lw_single(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2);
void lw_complete(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2);
void    lw_UPGMA(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2);
void    lw_WPGMA(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2);
void    lw_UPGMC(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2);
void    lw_WPGMC(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2);
void lw_flexible(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2); 

/*
 * Note on adding new classificatory strategies:
 * 
 * Additional classificatory strategies can be implemented by adding update
 * functions, add a case for it function setLWUpdate in order for it to be able
 * to provide its address to the clustering functions clust and constClust, and
 * add the new method to R language binding function constr.hclust in order for
 * the new type code to be given to the clustering functions.
 */

/* 
 * Lance and Williams algorithm: clustering functions (without or with the
 * spatial contiguity constraint) */
void       clust(int* n, int* membr, int* flag, int* ia, int* ib, double* crit,
                 double* diss0, int* method, double* par);
void  constClust(int* n, int* membr, int* flag, int* ia, int* ib, double* crit,
                 double* diss0, int* nl, int* linkl, int* method, double* par);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
|  C version of the hcass2 Fortran routine by                    |
|  F. Murtagh, ESA/ESO/STECF, Garching, June 1991                |
|                                                                |
|  Original description:                                         |
|  Given a HIERARCHIC CLUSTERING, described as a sequence of     |
|  agglomerations, prepare the seq. of aggloms. and "horiz."     |
|  order of objects for plotting the dendrogram using S routine  |
|  'plclust'.                                                    |
\+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void hcass2(int* n, int* ia, int* ib, int* iorder, int* iia, int* iib);

// R wrapper function (called from R using the .C() interface)
void cclust(int* n, int* merge, double* height, int* order, double* diss0,
            int* nl, int* linkl, int* method, double* par, int* type);

// Testing functions (called from R using the .C() interface during development)
// Not to be compiled for regular execution purposes
#ifdef testing

void R_getminlink(int* n, double* diss0, int* nl, int* linkl, int* nn_i,
                  int* nn_j, double* nn_dist);
void     R_getmin(int* n, int* flag, double* diss0, int* nn_i, int* nn_j,
                  double* nn_dist);

#endif

#ifdef show_links

void R_printlink(unsigned int nl, int* linkl);

#endif
