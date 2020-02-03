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

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
|                                                            |
|  Notes                                                     |
|  1. To access a value i, j, where i < j, the index k in    |
|     the distance data vector is calculated as follows:     |
|     k = j + i*n - (i+1)*(i+2)/2                            |
|     and where i > j, as follows:                           |
|     k = i + j+n - (j+1)((j+2)/2                            |
|                                                            |
\+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// C functions declarations

#include<R.h>
#include<Rmath.h>
#include"constr.hclust.h"

/*
 * Setup the update function (for ward.D2, method == 2, the distances are
 * squared and the criterion need to be square rooted at the end) */
void setLWUpdate(unsigned int n, int method, double* diss0,
                 void (**update)(int,int*,int*,double*,double*,unsigned int,unsigned int)) {
  unsigned int i;
  switch(method) {
  case 1:
    *update = &lw_Ward;
    break;
  case 2:
    *update = &lw_Ward;
    unsigned int len = n*(n-1)/2;
    for(i = 0; i < len; i++)
      diss0[i] *= diss0[i];
    break;
  case 3:
    *update = &lw_single;
    break;
  case 4:
    *update = &lw_complete;
    break;
  case 5:
    *update = &lw_UPGMA;
    break;
  case 6:
    *update = &lw_WPGMA;
    break;
  case 7:
    *update = &lw_UPGMC;
    break;
  case 8:
    *update = &lw_WPGMC;
    break;
  case 9:
    *update = &lw_flexible;
    break;
  default:
    error("Bad method number %d",method);
  }
  return;
}

// Finds the minimum of diss0 for unmasked locations
void getmin(unsigned int n, int* flag, double* diss0, unsigned int* nn_i,
            unsigned int* nn_j, double* nn_dist) {
  unsigned int i, j, k = 0;
  for(i = 0; i < (n - 1); i++) {
    if(flag[i]) {
      for(j = (i + 1); j < n; j++, k++) {
        if(flag[j]) {
          if(diss0[k] < *nn_dist) {
            *nn_dist = diss0[k];
            *nn_i = i;
            *nn_j = j;
          }
        }
      }
    } else
      k += (n - i - 1);
  }
  return;
}

/* 
 * Finds the minimum of diss0 for unmasked and connected locations using a link
 * list */
void getminlink(unsigned int n, double* diss0, unsigned int nl, int* linkl,
                unsigned int* nn_i, unsigned int* nn_j, double* nn_dist) {
  unsigned int i, ii, j, jj, k;
  for(i = 0, j = nl; i < nl; i++, j++) {
    if(linkl[i]!=linkl[j]) {
      ii = linkl[i];
      jj = linkl[j];
      if(ii<jj) {
        k = jj+ii*n-(ii+1)*(ii+2)/2;
        if(diss0[k] < *nn_dist) {
          *nn_dist = diss0[k];
          *nn_i = ii;
          *nn_j = jj;
        }
      } else {
        k = ii+jj*n-(jj+1)*(jj+2)/2;
        if(diss0[k] < *nn_dist) {
          *nn_dist = diss0[k];
          *nn_i = jj;
          *nn_j = ii;
        }
      }
    }
  }
  return;
}

// Merge the links shared between i2 and j2
void mergelink(unsigned int nl, int* linkl, unsigned int i2, unsigned int j2) {
  unsigned int i, j;
  for(i = 0, j = nl; i < nl; i++, j++) {
    if(linkl[i]==j2) linkl[i] = i2;
    if(linkl[j]==j2) linkl[j] = i2;
  }
  return;
}

/* 
 * A suite of update functions to recalculate the distances from the newly
 * aggregated clusters to the other elements */
void lw_Ward(int n, int* flag, int* membr, double* diss0, double* par,
             unsigned int i2, unsigned int j2) {
  unsigned int i, ind1, ind2, ind3;
  double temp;
  for(i = 0; i < n; i++) {
    if(flag[i] && (i != i2)) {
      if(i2 < i)
        ind1 = i + i2*n - (i2+1)*(i2+2)/2;
      else
        ind1 = i2 + i*n - (i+1)*(i+2)/2;
      if(j2 < i)
        ind2 = i + j2*n - (j2+1)*(j2+2)/2;
      else
        ind2 = j2 + i*n - (i+1)*(i+2)/2;
      ind3 = j2 + i2*n - (i2+1)*(i2+2)/2;  // D12 = diss0[ind3];
      temp = (membr[i2] + membr[i])*diss0[ind1] +
        (membr[j2] + membr[i])*diss0[ind2] - membr[i]*diss0[ind3];
      diss0[ind1] = temp / (membr[i2] + membr[j2] + membr[i]);
    }
  }
  return;
}

void lw_single(int n, int* flag, int* membr, double* diss0, double* par,
               unsigned int i2, unsigned int j2) {
  unsigned int i, ind1, ind2;
  for(i = 0; i < n; i++) {
    if(flag[i] && (i != i2)) {
      if(i2 < i)
        ind1 = i + i2*n - (i2+1)*(i2+2)/2;
      else
        ind1 = i2 + i*n - (i+1)*(i+2)/2;
      if(j2 < i)
        ind2 = i + j2*n - (j2+1)*(j2+2)/2;
      else
        ind2 = j2 + i*n - (i+1)*(i+2)/2;
      if(diss0[ind1] > diss0[ind2])
        diss0[ind1] = diss0[ind2];
    }
  }
  return;
}

void lw_complete(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2) {
  unsigned int i, ind1, ind2;
  for(i = 0; i < n; i++) {
    if(flag[i] && (i != i2)) {
      if(i2 < i)
        ind1 = i + i2*n - (i2+1)*(i2+2)/2;
      else
        ind1 = i2 + i*n - (i+1)*(i+2)/2;
      if(j2 < i)
        ind2 = i + j2*n - (j2+1)*(j2+2)/2;
      else
        ind2 = j2 + i*n - (i+1)*(i+2)/2;
      if(diss0[ind1] < diss0[ind2])
        diss0[ind1] = diss0[ind2];
    }
  }
  return;
}

void lw_UPGMA(int n, int* flag, int* membr, double* diss0, double* par,
              unsigned int i2, unsigned int j2) {
  unsigned int i, ind1, ind2;
  for(i = 0; i < n; i++) {
    if(flag[i] && (i != i2)) {
      if(i2 < i)
        ind1 = i + i2*n - (i2+1)*(i2+2)/2;
      else
        ind1 = i2 + i*n - (i+1)*(i+2)/2;
      if(j2 < i)
        ind2 = i + j2*n - (j2+1)*(j2+2)/2;
      else
        ind2 = j2 + i*n - (i+1)*(i+2)/2;
      diss0[ind1] = (membr[i2]*diss0[ind1] + membr[j2]*diss0[ind2]) /
        (membr[i2] + membr[j2]);
    }
  }
  return;
}

void lw_WPGMA(int n, int* flag, int* membr, double* diss0, double* par,
              unsigned int i2, unsigned int j2) {
  unsigned int i, ind1, ind2;
  for(i = 0; i < n; i++) {
    if(flag[i] && (i != i2)) {
      if(i2 < i)
        ind1 = i + i2*n - (i2+1)*(i2+2)/2;
      else
        ind1 = i2 + i*n - (i+1)*(i+2)/2;
      if(j2 < i)
        ind2 = i + j2*n - (j2+1)*(j2+2)/2;
      else
        ind2 = j2 + i*n - (i+1)*(i+2)/2;
      diss0[ind1] = 0.5*diss0[ind1] + 0.5*membr[j2]*diss0[ind2];
    }
  }
  return;
}

void lw_UPGMC(int n, int* flag, int* membr, double* diss0, double* par,
              unsigned int i2, unsigned int j2) {
  unsigned int i, ind1, ind2, ind3;
  double temp;
  for(i = 0; i < n; i++) {
    if(flag[i] && (i != i2)) {
      if(i2 < i)
        ind1 = i + i2*n - (i2+1)*(i2+2)/2;
      else
        ind1 = i2 + i*n - (i+1)*(i+2)/2;
      if(j2 < i)
        ind2 = i + j2*n - (j2+1)*(j2+2)/2;
      else
        ind2 = j2 + i*n - (i+1)*(i+2)/2;
      ind3 = j2 + i2*n - (i2+1)*(i2+2)/2;
      temp = membr[i2]*diss0[ind1] + membr[j2]*diss0[ind2] -
        membr[i2]*membr[j2]*diss0[ind3]/(membr[i2] + membr[j2]);
      diss0[ind1] = temp / (membr[i2] + membr[j2]);
    }
  }
  return;
}

void lw_WPGMC(int n, int* flag, int* membr, double* diss0, double* par,
              unsigned int i2, unsigned int j2) {
  unsigned int i, ind1, ind2, ind3;
  for(i = 0; i < n; i++) {
    if(flag[i] && (i != i2)) {
      if(i2 < i)
        ind1 = i + i2*n - (i2+1)*(i2+2)/2;
      else
        ind1 = i2 + i*n - (i+1)*(i+2)/2;
      if(j2 < i)
        ind2 = i + j2*n - (j2+1)*(j2+2)/2;
      else
        ind2 = j2 + i*n - (i+1)*(i+2)/2;
      ind3 = j2 + i2*n - (i2+1)*(i2+2)/2;
      diss0[ind1] = 0.5*diss0[ind1] + 0.5*diss0[ind2] - 0.25*diss0[ind3];
    }
  }
  return;
}

void lw_flexible(int n, int* flag, int* membr, double* diss0, double* par,
                 unsigned int i2, unsigned int j2) {
  // alpha: par[0] ; beta: par[1]
  unsigned int i, ind1, ind2, ind3;
  for(i = 0; i < n; i++) {
    if(flag[i] && (i != i2)) {
      if(i2 < i)
        ind1 = i + i2*n - (i2+1)*(i2+2)/2;
      else
        ind1 = i2 + i*n - (i+1)*(i+2)/2;
      if(j2 < i)
        ind2 = i + j2*n - (j2+1)*(j2+2)/2;
      else
        ind2 = j2 + i*n - (i+1)*(i+2)/2;
      ind3 = j2 + i2*n - (i2+1)*(i2+2)/2;
      diss0[ind1] = par[0]*diss0[ind1] + par[0]*diss0[ind2] - par[1]*diss0[ind3];
    }
  }
  return;
}

/* 
 * Plain Lance and Williams clustering function (without spatial contiguity
 * constraints) */
void clust(int* n, int* membr, int* flag, int* ia, int* ib, double* crit,
           double* diss0, int* method, double* par) {
  unsigned int i, nn_i, nn_j, i2, j2, ncl;
  double nn_dist;
  void (*update)(int,int*,int*,double*,double*,unsigned int,unsigned int);
  setLWUpdate((unsigned int)(*n),*method,diss0,&update);
  for(i = 0; i < *n; i++) {
    membr[i] = 1;
    flag[i] = 1;
  }
  for(ncl = *n; ncl > 1;) {
    if(!(ncl%INTMOD))
      R_CheckUserInterrupt();
    nn_dist = R_PosInf;
    getmin((unsigned int)(*n),flag,diss0,&nn_i,&nn_j,&nn_dist);
    ncl--;
    i2 = (nn_i<nn_j) ? nn_i : nn_j;
    j2 = (nn_i>nn_j) ? nn_i : nn_j;
    ia[*n-ncl-1] = i2;
    ib[*n-ncl-1] = j2;
    crit[*n-ncl-1] = nn_dist;
    flag[j2] = false;
    update(*n,flag,membr,diss0,par,i2,j2);
    membr[i2] += membr[j2];
  }
  for(i = 0; i < (*n - 1); i++) {
    ia[i]++;
    ib[i]++;
  }
  if(*method == 2)
    for(i = 0; i < (*n - 1); i++)
      crit[i] = sqrt(crit[i]);
  return;
}

// Lance and Williams clustering function with the spatial contiguity constraint
void constClust(int* n, int* membr, int* flag, int* ia, int* ib, double* crit,
                double* diss0, int* nl, int* linkl, int* method, double* par) {
  void (*update)(int,int*,int*,double*,double*,unsigned int,unsigned int);
  unsigned int i, ncl, nn_i, nn_j, i2, j2;
  double nn_dist;
  setLWUpdate((unsigned int)(*n),*method,diss0,&update);
  for(i = 0; i < *n; i++) {
    membr[i] = 1;
    flag[i] = 1;
  }
  for(ncl = *n; ncl > 1;) {
    if(!(ncl%INTMOD))
      R_CheckUserInterrupt();
    nn_dist = R_PosInf;
    getminlink((unsigned int)(*n),diss0,(unsigned int)(*nl),linkl,&nn_i,&nn_j,
               &nn_dist);
    ncl--;
    i2 = (nn_i<nn_j) ? nn_i : nn_j;
    j2 = (nn_i>nn_j) ? nn_i : nn_j;
    ia[*n-ncl-1] = i2;
    ib[*n-ncl-1] = j2;
    crit[*n-ncl-1] = nn_dist;
    flag[j2] = false;
    update(*n,flag,membr,diss0,par,i2,j2);
    membr[i2] += membr[j2];
    mergelink((unsigned int)(*nl),linkl,i2,j2);
  }
  for(i = 0; i < (*n - 1); i++) {
    ia[i]++;
    ib[i]++;
  }
  /* For ward.D2, distances were squared and the criterion thus has to be
   * squared rooted */
  if(*method == 2)
    for(i = 0; i < (*n - 1); i++)
      crit[i] = sqrt(crit[i]);
  return;
}

// C version of the hcass2 Fortran routine
void hcass2(int *n, int *ia, int *ib, int *iorder, int *iia, int *iib) {
  int i, j, k, k1, k2, loc, dobreak;
  for(i = 0; i < (*n - 1); i++) {
    iia[i] = ia[i];
    iib[i] = ib[i];
  }
  for(i = 0; i < (*n - 2); i++) {
    k = (ia[i] < ib[i]) ? ia[i] : ib[i];
    for(j = i + 1; j < (*n - 1); j++) {
      if(ia[j] == k)
        iia[j] = -(i+1);
      if(ib[j] == k)
        iib[j] = -(i+1);
    }
  }
  for(i = 0; i < (*n - 1); i++) {
    iia[i] = -iia[i];
    iib[i] = -iib[i];
  }
  for(i = 0; i < (*n - 1); i++) {
    if((iia[i] > 0) && (iib[i] < 0)) {
      k = iia[i];
      iia[i] = iib[i];
      iib[i] = k;
    }
    if((iia[i] > 0) && (iib[i] > 0)) {
      k1 = (iia[i]<iib[i]) ? iia[i] : iib[i];
      k2 = (iia[i]>iib[i]) ? iia[i] : iib[i];
      iia[i] = k1;
      iib[i] = k2;
    }
  }
  iorder[0] = iia[*n-2];
  iorder[1] = iib[*n-2];
  loc = 2;
  for(i = *n - 3; i >= 0; i--) {
    for(j = 0; j < loc; j++) {
      dobreak = false;
      if(iorder[j] == (i+1)) {
        iorder[j] = iia[i];
        if(j == (loc-1)) {
          loc++;
          iorder[loc-1] = iib[i];
          dobreak = true;
        }
        loc++;
        for(k = (loc-2); k > (j+1); k--)
          iorder[k] = iorder[k-1];
        iorder[j+1] = iib[i];
        dobreak = true;
      }
      if(dobreak)
        break;
    }
  }
  for(i = 0; i < *n; i++)
    iorder[i] = -iorder[i];
  return;
}

// R wrapper function
void cclust(int* n, int* merge, double* height, int* order, double* diss0,
            int* nl, int* linkl, int* method, double* par, int* type) {
  int* flag = (int*)R_alloc(*n,sizeof(int));
  int* membr = (int*)R_alloc(*n,sizeof(int));
  int* ia = (int*)R_alloc(*n-1,sizeof(int));
  int* ib = (int*)R_alloc(*n-1,sizeof(int));
  int* linkc;
  switch(*type) {
  case 0:
    clust(n,membr,flag,ia,ib,height,diss0,method,par);
    break;
  case 1:
    for(int i = 0, j = *nl; i < *nl; i++, j++) {
      linkl[i]--;
      linkl[j]--;
    }
#ifdef show_links
    printf("\nBefore:\n");
    R_printlink((unsigned int)(*nl), linkl);
#endif
    constClust(n,membr,flag,ia,ib,height,diss0,nl,linkl,method,par);
#ifdef show_links
    printf("\nAfter:\n");
    R_printlink((unsigned int)(*nl), linkl);
#endif
    break;
  case 2:
    *nl = *n - 1;
    linkc = (int*)R_alloc(2*(*nl),sizeof(int));
    for(int i = 0, j = *nl; i < *nl; i++, j++) {
      linkc[i] = i;
      linkc[j] = i+1;
    }
#ifdef show_links
    printf("Before:\n");
    R_printlink((unsigned int)(*nl), linkc);
#endif
    constClust(n,membr,flag,ia,ib,height,diss0,nl,linkc,method,par);
#ifdef show_links
    printf("After:\n");
    R_printlink((unsigned int)(*nl), linkc);
#endif
    break;
  default:
    error("Bad method number %d",*type);
  }
  hcass2(n,ia,ib,order,merge,&merge[*n-1]);
  // Free(ib);
  // Free(ia);
  // Free(flag);
  // Free(membr);
  return;
}

// Testing functions (called from R using the .C() interface during development)
#ifdef testing

void R_getminlink(int* n, double* diss0, int* nl, int* linkl, int* nn_i,
                  int* nn_j, double* nn_dist) {
  *nn_dist = R_PosInf;
  *nn_i = 0;
  *nn_j = 0;
  getminlink((unsigned int)(*n),diss0,(unsigned int)(*nl),linkl,
             (unsigned int*)nn_i,(unsigned int*)nn_j,nn_dist);
  (*nn_i)++;
  (*nn_j)++;
  return;
}

void R_getmin(int* n, int* flag, double* diss0, int* nn_i, int* nn_j,
              double* nn_dist) {
  *nn_dist = R_PosInf;
  *nn_i = 0;
  *nn_j = 0;
  getmin((unsigned int)(*n),flag,diss0,(unsigned int*)nn_i,
         (unsigned int*)nn_j,nn_dist);
  (*nn_i)++;
  (*nn_j)++;
  return;
}

#endif

#ifdef show_links

void R_printlink(unsigned int nl, int* linkl) {
  int i, j;
  for(i = 0, j = nl; i < nl; i++, j++)
    printf("%d -> %d\n",linkl[i],linkl[j]);
  return;
}

#endif
