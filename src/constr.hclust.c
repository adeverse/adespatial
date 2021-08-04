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
|  August 2018 - February 2020                               |
|                                                            |
|  The present implementation has been vastly inspirered by  |
|  the original HCLUST Fortran code by                       |
|  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       |
|  with modifications for R: Ross Ihaka, Dec 1996            |
|                            Fritz Leisch, Jun 2000          |
|                                                            |
|  Also contains:                                            |
|  Least squares (non distance-based) agglomerative          |
|  clustering with or without the spatial contiguity         |
|  constraint.                                               |
|                                                            |
|  Performance was enhanced by carefully buffering the       |
|  identity (unconstrained case) and the values of nearest   |
|  neighbour dissimilarities (both unconstrained and         |
|  constrained cases) for each observations (unconstrained   |
|  case) or links (constrained case). This type of           |
|  performance enhancement is present in recent hclust       |
|  implementations and has been originally proposed by       |
|  Langfelder and Horvath (2012, Journal of Statistical      |
|  Software, 46(11), 1-17) for the unconstrained case in     |
|  package flashClust.                                       |
|                                                            |
\+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
|                                                            |
|  Notes                                                     |
|  1. To access a value i, j, where i < j, the index k in    |
|     the distance data vector is calculated as follows:     |
|     k = j + i*n - (i+1)*(i+2)/2                            |
|     and where i > j, as follows:                           |
|     k = i + j+n - (j+1)*(j+2)/2                            |
|                                                            |
\+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// C functions declarations

#include<R.h>
#include<Rmath.h>
#include"constr.hclust.h"

/*
 * Set the update function (for ward.D2, method == 2, the distances are squared
 * and the criterion need to be square rooted at the end) */
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

void initNNlist(unsigned int n, double* diss0, unsigned int* nn_idx,
                double* nn_diss, unsigned int* min_idx) {
  unsigned int i, j, k, nn;
  double min_diss, global_min = R_PosInf;
  for(i = 0, k = 0; i < (n - 1); i++) {
    min_diss = R_PosInf;
    for(j = (i + 1); j < n; j++, k++) {
      if(diss0[k] < min_diss) {
        nn = j;
        min_diss = diss0[k];
      }
    }
    nn_idx[i] = nn;
    nn_diss[i] = min_diss;
    if(min_diss < global_min) {
      *min_idx = i;
      global_min = min_diss;
    }
  }
  return;
}

void updateNNlist(unsigned int n, int* flag, double* diss0,
                  unsigned int* nn_idx, double* nn_diss, unsigned int i) {
  unsigned int j, k, nn;
  double min_diss = R_PosInf;
  for(j = (i + 1), k = j + i*n - (i + 1)*(i + 2)/2; j < n; j++, k++) {
    if(flag[j])
      if(diss0[k] < min_diss) {
        nn = j;
        min_diss = diss0[k];
      }
  }
  nn_idx[i] = nn;
  nn_diss[i] = min_diss;
  return;
}

void fixNNlist(unsigned int n, double* diss0, unsigned int* nn_idx,
               double* nn_diss, unsigned int i, unsigned int j) {
  unsigned int k = j + i*n - (i + 1)*(i + 2)/2;
  double min_diss;
  min_diss = diss0[k];
  if(min_diss < nn_diss[i]) {
    nn_idx[i] = j;
    nn_diss[i] = min_diss;
  }
  return;
}

void initNNlink(unsigned int n, double* diss0, unsigned int nl, int* linkl,
                unsigned int* nn_i, unsigned int* nn_j, double* nn_diss,
                double* min_diss) {
  unsigned int i, j, ii, jj, k;
  *min_diss = R_PosInf;
  for(i = 0, j = nl; i < nl; i++, j++) {
    if(linkl[i] > linkl[j]) {
      jj = linkl[i];
      ii = linkl[j];
      linkl[i] = ii;
      linkl[j] = jj;
    } else {
      ii = linkl[i];
      jj = linkl[j];
    }
    k = jj + ii*n - (ii + 1)*(ii + 2)/2;
    nn_diss[i] = diss0[k];
    if(nn_diss[i] < *min_diss) {
      *min_diss = nn_diss[i];
      *nn_i = ii;
      *nn_j = jj;
    }
  }
  return;
}

void updateNNlink(unsigned int n, double* diss0, unsigned int nl, int* linkl,
                  unsigned int* nn_i, unsigned int* nn_j, double* nn_diss,
                  double* min_diss, unsigned int i2, unsigned int j2) {
  unsigned int i, ii, j, jj, k, changed;
  *min_diss = R_PosInf;
  for(i = 0, j = nl; i < nl; i++, j++) {
    if(linkl[i] != linkl[j]) {
      changed = false;
      if(linkl[i] == j2) {
        linkl[i] = i2;
        changed = true;
      }
      if(linkl[j] == j2) {
        linkl[j] = i2;
        changed = true;
      }
      if(linkl[i] != linkl[j]) {
        if(changed) {
          if(linkl[i] > linkl[j]) {
            ii = linkl[j];
            linkl[j] = linkl[i];
            linkl[i] = ii;
            jj = linkl[j];
          } else {
            ii = linkl[i];
            jj = linkl[j];
          }
          k = jj + ii*n - (ii + 1)*(ii + 2)/2;
          nn_diss[i] = diss0[k];
        } else if((linkl[i] == i2) || (linkl[j] == i2)) {
          ii = linkl[i];
          jj = linkl[j];
          k = jj + ii*n - (ii + 1)*(ii + 2)/2;
          nn_diss[i] = diss0[k];
        }
        if(nn_diss[i] < *min_diss) {
          *min_diss = nn_diss[i];
          *nn_i = linkl[i];
          *nn_j = linkl[j];
        }
      }
    }
  }
  return;
}

#ifdef with_LS
void initNNlistLS(unsigned int n, unsigned int m, double* x, double* xx,
                  unsigned int* nn_idx, double* nn_diss,
                  unsigned int* min_idx) {
  unsigned int i, j, k, ind_i, ind_j;
  double acc, tmp, min_diss = R_PosInf;
  for(i = 0; i < (n - 1); i++) {
    nn_diss[i] = R_PosInf;
    for(j = (i + 1); j < n; j++) {
      ind_i = i;
      ind_j = j;
      acc = 0.0;
      for(k = 0; k < m; k++, ind_i += n, ind_j += n) {
        tmp = x[ind_i] + x[ind_j];
        tmp *= tmp;
        acc += (xx[ind_i] + xx[ind_j]) - tmp/2.0;
      }
      if(acc < nn_diss[i]) {
        nn_idx[i] = j;
        nn_diss[i] = acc;
      }
    }
    if(nn_diss[i] < min_diss) {
      *min_idx = i;
      min_diss = nn_diss[i];
    }
  }
  return;
}

void updateNNlistLS(unsigned int n, int* membr, int* flag, unsigned int m,
                    double* x, double* xx, unsigned int* nn_idx,
                    double* nn_diss, unsigned int i) {
  unsigned int j, k, ind_i, ind_j;
  double acc, tmp;
  nn_diss[i] = R_PosInf;
  for(j = (i + 1); j < n; j++)
    if(flag[j]) {
      ind_i = i;
      ind_j = j;
      acc = 0.0;
      for(k = 0; k < m; k++, ind_i += n, ind_j += n) {
        tmp = x[ind_i] + x[ind_j];
        tmp *= tmp;
        acc += (xx[ind_i] + xx[ind_j]) - tmp/(membr[i] + membr[j]);
      }
      if(acc < nn_diss[i]) {
        nn_idx[i] = j;
        nn_diss[i] = acc;
      }
    }
  return;
}

void initNNlinkLS(unsigned int n, int* membr, unsigned int m, double* x,
                  double* xx, unsigned int nl, int* linkl, unsigned int* nn_i,
                  unsigned int* nn_j, double* nn_diss, double* min_diss) {
  unsigned int i, j, ii, jj, k, ind_i, ind_j;
  double acc, tmp;
  *min_diss = R_PosInf;
  for(i = 0, j = nl; i < nl; i++, j++) {
    if(linkl[i] > linkl[j]) {
      jj = linkl[i];
      ii = linkl[j];
      linkl[i] = ii;
      linkl[j] = jj;
    } else {
      ii = linkl[i];
      jj = linkl[j];
    }
    ind_i = ii;
    ind_j = jj;
    acc = 0.0;
    for(k = 0; k < m; k++, ind_i += n, ind_j += n) {
      tmp = x[ind_i] + x[ind_j];
      tmp *= tmp;
      acc += (xx[ind_i] + xx[ind_j]) - tmp/(membr[ii] + membr[jj]);
    }
    nn_diss[i] = acc;
    if(nn_diss[i] < *min_diss) {
      *min_diss = nn_diss[i];
      *nn_i = ii;
      *nn_j = jj;
    }
  }
  return;
}

void updateNNlinkLS(unsigned int n, int* membr, unsigned int m, double* x,
                    double* xx, unsigned int nl, int* linkl, unsigned int* nn_i,
                    unsigned int* nn_j, double* nn_diss, double* min_diss,
                    unsigned int i2, unsigned int j2) {
  unsigned int i, ii, j, jj, k, changed, ind_i, ind_j;
  double acc, tmp;
  *min_diss = R_PosInf;
  for(i = 0, j = nl; i < nl; i++, j++) {
    if(linkl[i] != linkl[j]) {
      changed = false;
      if(linkl[i] == j2) {
        linkl[i] = i2;
        changed = true;
      }
      if(linkl[j] == j2) {
        linkl[j] = i2;
        changed = true;
      }
      if(linkl[i] != linkl[j]) {
        if(changed) {
          if(linkl[i] > linkl[j]) {
            ii = linkl[j];
            linkl[j] = linkl[i];
            linkl[i] = ii;
            jj = linkl[j];
          } else {
            ii = linkl[i];
            jj = linkl[j];
          }
          ind_i = ii;
          ind_j = jj;
          acc = 0.0;
          for(k = 0; k < m; k++, ind_i += n, ind_j += n) {
            tmp = x[ind_i] + x[ind_j];
            tmp *= tmp;
            acc += (xx[ind_i] + xx[ind_j]) - tmp/(membr[ii] + membr[jj]);
          }
          nn_diss[i] = acc;
        } else if((linkl[i] == i2) || (linkl[j] == i2)) {
          ii = linkl[i];
          jj = linkl[j];
          ind_i = ii;
          ind_j = jj;
          acc = 0.0;
          for(k = 0; k < m; k++, ind_i += n, ind_j += n) {
            tmp = x[ind_i] + x[ind_j];
            tmp *= tmp;
            acc += (xx[ind_i] + xx[ind_j]) - tmp/(membr[ii] + membr[jj]);
          }
          nn_diss[i] = acc;
        }
        if(nn_diss[i] < *min_diss) {
          *min_diss = nn_diss[i];
          *nn_i = linkl[i];
          *nn_j = linkl[j];
        }
      }
    }
  }
  return;
}

void updateLS(unsigned int n, int* membr, unsigned int m, double* x, double* xx,
             unsigned int i2, unsigned int j2) {
  unsigned int k;
  for(k = 0; k < m; k++, i2 += n, j2 += n) {
    x[i2] += x[j2];
    xx[i2] += xx[j2];
  }
  return;
}

double getdistLS(unsigned int n, int* membr, unsigned int m, double* x,
                 unsigned int i2, unsigned int j2, unsigned int squared) {
  unsigned int k, ii2 = i2, jj2 = j2;
  double tmp, acc = 0.0;
  for(k = 0; k < m; k++, ii2 += n, jj2 += n) {
    tmp = x[ii2]/membr[i2] - x[jj2]/membr[j2];
    tmp *= tmp;
    acc += tmp;
  }
  return(squared ? acc : sqrt(acc));
}
#endif

/* A suite of update functions to recalculate the distances from the newly
 * aggregated clusters to the other elements.
 * New methods can be implemented by adding new functions in the C code, adding
 * the necessary code to dispatch them in setLWUpdate(), and update the R
 * language header to handle the new method properly. Notably, function's
 * argument beta may then need to pass >1 value. One may need to update the
 * function's formal arguments while making sure that it doesn't break backward
 * compatibility with the methods that are already implemented.
 */
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
      diss0[ind1] = 0.5*diss0[ind1] + 0.5*diss0[ind2];
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

/* Plain Lance and Williams clustering function (without spatial contiguity
 * constraints) flashClust version */
void clust(int* n, int* membr, int* flag, int* ia, int* ib, double* crit,
           double* diss0, unsigned int* nn_idx, double* nn_diss, int* method,
           double* par, int fastUpdate) {
  unsigned int i, k, nn_i, nn_j, i2, j2, ncl, min_idx;
  double min_diss;
  void (*update)(int,int*,int*,double*,double*,unsigned int,unsigned int);
  setLWUpdate((unsigned int)(*n),*method,diss0,&update);
  for(i = 0; i < *n; i++)
    flag[i] = true;
  initNNlist((unsigned int)(*n),diss0,nn_idx,nn_diss,&nn_i);
  nn_j = nn_idx[nn_i];
  min_diss = nn_diss[nn_i];
  for(ncl = *n; ncl > 1;) {
    if(!(ncl%INTMOD))
      R_CheckUserInterrupt();
    ncl--;
    i2 = nn_i;
    j2 = nn_j;
    ia[*n-ncl-1] = i2;
    ib[*n-ncl-1] = j2;
    crit[*n-ncl-1] = min_diss;
    flag[j2] = false;
    update(*n,flag,membr,diss0,par,i2,j2);
    membr[i2] += membr[j2];
    if(ncl > 1) {
      if(fastUpdate) {
        min_diss = R_PosInf;
        for(i = 0; i < (*n - 1); i++)
          if(flag[i]) {
            if((nn_idx[i] == i2) || (nn_idx[i] == j2))
              updateNNlist((unsigned int)(*n),flag,diss0,nn_idx,nn_diss,i);
            if(nn_diss[i] < min_diss) {
              nn_i = i;
              min_diss = nn_diss[i];
            }
          }
        nn_j = nn_idx[nn_i];
      } else {
        min_diss = R_PosInf;
        for(i = 0; i < (*n - 1); i++)
          if(flag[i]) {
            if((i == i2) || (nn_idx[i] == i2) || (nn_idx[i] == j2))
              updateNNlist((unsigned int)(*n),flag,diss0,nn_idx,nn_diss,i);
            else
              if(i < i2)
                fixNNlist((unsigned int)(*n),diss0,nn_idx,nn_diss,i,i2);
            if(nn_diss[i] < min_diss) {
              nn_i = i;
              min_diss = nn_diss[i];
            }
          }
        nn_j = nn_idx[nn_i];
      }
    }
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

/* Lance and Williams clustering function with the spatial contiguity
 * constraint */
void constClust(int* n, int* membr, int* flag, int* ia, int* ib, double* crit,
                double* diss0, double* nn_diss, int* nl, int* linkl,
                int* method, double* par) {
  void (*update)(int,int*,int*,double*,double*,unsigned int,unsigned int);
  unsigned int i, ncl, nn_i, nn_j, i2, j2;
  double min_diss;
  setLWUpdate((unsigned int)(*n),*method,diss0,&update);
  for(i = 0; i < *n; i++)
    flag[i] = true;
  initNNlink((unsigned int)(*n),diss0,(unsigned int)(*nl),linkl,&nn_i,&nn_j,
             nn_diss,&min_diss);
  for(ncl = *n; ncl > 1;) {
    if(!(ncl%INTMOD))
      R_CheckUserInterrupt();
    ncl--;
    //
    /* When all links have been depleted, merge the remaining disjoint clusters,
     * from that with smallest index > 0 to that with the largest index, to
     * cluster 0 while inserting NA as the merging dissimilarity.*/
    if(min_diss == R_PosInf) {
      for(i = 1; (i < *n) && !flag[i]; i++)
        ;
      nn_i = 0;
      nn_j = i;
      min_diss = NA_REAL;
    }
    //
    i2 = nn_i;
    j2 = nn_j;
    ia[*n-ncl-1] = i2;
    ib[*n-ncl-1] = j2;
    crit[*n-ncl-1] = min_diss;
    //
    flag[j2] = false;
    update(*n,flag,membr,diss0,par,i2,j2);
    membr[i2] += membr[j2];
    if(ncl > 1) 
      updateNNlink((unsigned int)(*n),diss0,(unsigned int)(*nl),linkl,
                   &nn_i,&nn_j,nn_diss,&min_diss,i2,j2);
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

#ifdef with_LS
/* Main routine for plain least squares agglomerative clustering without a
 * contiguity constraint. */
void clustLS(int* n, int* membr, int* flag, int* ia, int* ib, double* crit,
             int* m, double* x, double* xx, unsigned int* nn_idx,
             double* nn_diss, int* out) {
  unsigned int i, j, k, nn_i, nn_j, i2, j2, ncl;
  double min_diss;
  for(i = 0; i < *n; i++) {
    membr[i] = 1;
    flag[i] = true;
  }
  initNNlistLS((unsigned int)(*n),(unsigned int)(*m),x,xx,nn_idx,nn_diss,&nn_i);
  nn_j = nn_idx[nn_i];
  min_diss = nn_diss[nn_i];
  for(ncl = *n; ncl > 1;) {
    if(!(ncl%INTMOD))
      R_CheckUserInterrupt();
    ncl--;
    i2 = nn_i;
    j2 = nn_j;
    ia[*n-ncl-1] = i2;
    ib[*n-ncl-1] = j2;
    switch(*out) {
    case 1:
      crit[*n-ncl-1] = sqrt(min_diss);
      break;
    case 2:
      crit[*n-ncl-1] = min_diss;
      break;
    case 3:
      crit[*n-ncl-1] = getdistLS((unsigned int)(*n),membr,(unsigned int)(*m),x,
                                 i2,j2,false);
      break;
    case 4:
      crit[*n-ncl-1] = getdistLS((unsigned int)(*n),membr,(unsigned int)(*m),x,
                                 i2,j2,true);
      break;
    default:
      error("Unknown output type: %d\n",*out);
    }
    flag[j2] = false;
    updateLS((unsigned int)(*n),membr,(unsigned int)(*m),x,xx,i2,j2);
    membr[i2] += membr[j2];
    if(ncl > 1) {
      min_diss = R_PosInf;
      for(i = 0; i < (*n - 1); i++)
        if(flag[i]) {
          if((nn_idx[i] == i2) || (nn_idx[i] == j2))
            updateNNlistLS((unsigned int)(*n),membr,flag,(unsigned int)(*m),x,
                           xx,nn_idx,nn_diss,i);
          if(nn_diss[i] < min_diss) {
            nn_i = i;
            min_diss = nn_diss[i];
          }
        }
      nn_j = nn_idx[nn_i];
    }
  }
  for(i = 0; i < (*n - 1); i++) {
    ia[i]++;
    ib[i]++;
  }
  return;
}

/* Main routine for contiguity-constrained least squares agglomerative
 * clustering. */
void constClustLS(int* n, int* membr, int* ia, int* ib, double* crit, int* m,
                  double* x, double* xx, int* nl, int* linkl, double* nn_diss,
                  int* out) {
  unsigned int i, j, nn_i, nn_j, i2, j2, ncl;
  double min_diss;
  for(i = 0; i < *n; i++)
    membr[i] = 1;
  initNNlinkLS((unsigned int)(*n),membr,(unsigned int)(*m),x,xx,
               (unsigned int)(*nl),linkl,&nn_i,&nn_j,nn_diss,&min_diss);
  for(ncl = *n; ncl > 1;) {
    if(!(ncl%INTMOD))
      R_CheckUserInterrupt();
    ncl--;
    if(min_diss == R_PosInf) {
      for(i = 1; (i < *n) && !membr[i]; i++)
        ;
      nn_i = 0;
      nn_j = i;
    }
    i2 = nn_i;
    j2 = nn_j;
    ia[*n-ncl-1] = i2;
    ib[*n-ncl-1] = j2;
    if(min_diss == R_PosInf)
      crit[*n-ncl-1] = NA_REAL;
    else
      switch(*out) {
      case 1:
        crit[*n-ncl-1] = sqrt(min_diss);
        break;
      case 2:
        crit[*n-ncl-1] = min_diss;
        break;
      case 3:
        crit[*n-ncl-1] = getdistLS((unsigned int)(*n),membr,(unsigned int)(*m),
                                   x,i2,j2,false);
        break;
      case 4:
        crit[*n-ncl-1] = getdistLS((unsigned int)(*n),membr,(unsigned int)(*m),
                                   x,i2,j2,true);
        break;
      default:
        error("Unknown output type: %d\n",*out);
      }
    // This function has no flag (unlike all others).
    updateLS((unsigned int)(*n),membr,(unsigned int)(*m),x,xx,i2,j2);
    membr[i2] += membr[j2];
    membr[j2] = 0;  // Here, membr is used as flag
    updateNNlinkLS((unsigned int)(*n),membr,(unsigned int)(*m),x,xx,
                   (unsigned int)(*nl),linkl,&nn_i,&nn_j,nn_diss,&min_diss,
                   i2,j2);
  }
  for(i = 0; i < (*n - 1); i++) {
    ia[i]++;
    ib[i]++;
  }
  return;
}
#endif

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

// Service function to be called by R wrapper function: constr.hclust
// (flashClust version)
void cclust(int* n, int* merge, double* height, int* order, double* diss0,
            int* nl, int* linkl, int* method, double* par, int* type,
            int* membr) {
  int* flag = (int*)R_alloc(*n,sizeof(int));
  int* ia = (int*)R_alloc(*n-1,sizeof(int));
  int* ib = (int*)R_alloc(*n-1,sizeof(int));
  unsigned int* nn_idx;
  double* nn_diss;
  int* linkc;
  switch(*type) {
  case 0:
    nn_idx = (unsigned int*)R_alloc(*n-1,sizeof(unsigned int));
    nn_diss = (double*)R_alloc(*n-1,sizeof(double));
    clust(n,membr,flag,ia,ib,height,diss0,nn_idx,nn_diss,method,par,
          ((*method==7) || (*method==8)) ? false : true);
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
    nn_diss = (double*)R_alloc(*nl,sizeof(double));
    constClust(n,membr,flag,ia,ib,height,diss0,nn_diss,nl,linkl,method,par);
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
    nn_diss = (double*)R_alloc(*nl,sizeof(double));
#ifdef show_links
    printf("Before:\n");
    R_printlink((unsigned int)(*nl), linkc);
#endif
    constClust(n,membr,flag,ia,ib,height,diss0,nn_diss,nl,linkc,method,par);
#ifdef show_links
    printf("After:\n");
    R_printlink((unsigned int)(*nl), linkc);
#endif
    break;
  default:
    error("Bad method number %d",*type);
  }
  hcass2(n,ia,ib,order,merge,&merge[*n-1]);
  return;
}

#ifdef with_LS
// Service function to be called by R wrapper function: constr.lshclust
void cclustLS(int* n, int* merge, double* height, int* order, int* m,
              double* x, int* nl, int* linkl, int* type, int* out) {
  int* flag;
  int* membr = (int*)R_alloc(*n,sizeof(int));
  int* ia = (int*)R_alloc(*n-1,sizeof(int));
  int* ib = (int*)R_alloc(*n-1,sizeof(int));
  double* xx = (double*)R_alloc((*n)*(*m),sizeof(double));
  unsigned int* nn_idx;
  double* nn_diss;
  for(int i = 0; i < (*n)*(*m); i++)
    xx[i] = x[i]*x[i];
  int* linkc;
  switch(*type) {
  case 0:
    flag = (int*)R_alloc(*n,sizeof(int));
    nn_idx = (unsigned int*)R_alloc(*n-1,sizeof(unsigned int));
    nn_diss = (double*)R_alloc(*n-1,sizeof(double));
#ifdef show_links
    printf("\nNo links\n");
#endif
    clustLS(n,membr,flag,ia,ib,height,m,x,xx,nn_idx,nn_diss,out);
    break;
  case 1:
    for(int i = 0, j = *nl; i < *nl; i++, j++) {
      linkl[i]--;
      linkl[j]--;
    }
    nn_diss = (double*)R_alloc(*nl,sizeof(double));
#ifdef show_links
    printf("\nBefore:\n");
    R_printlink((unsigned int)(*nl),linkl);
#endif
    constClustLS(n,membr,ia,ib,height,m,x,xx,nl,linkl,nn_diss,out);
#ifdef show_links
    printf("\nAfter:\n");
    R_printlink((unsigned int)(*nl),linkl);
#endif
    break;
  case 2:
    *nl = *n - 1;
    linkc = (int*)R_alloc(2*(*nl),sizeof(int));
    for(int i = 0, j = *nl; i < *nl; i++, j++) {
      linkc[i] = i;
      linkc[j] = i+1;
    }
    nn_diss = (double*)R_alloc(*nl,sizeof(double));
#ifdef show_links
    printf("Before:\n");
    R_printlink((unsigned int)(*nl),linkc);
#endif
    constClustLS(n,membr,ia,ib,height,m,x,xx,nl,linkc,nn_diss,out);
#ifdef show_links
    printf("After:\n");
    R_printlink((unsigned int)(*nl),linkc);
#endif
    break;
  default:
    error("Bad method number %d",*type);
  }
  hcass2(n,ia,ib,order,merge,&merge[*n-1]);
  return;
}
#endif

// Testing functions (called from R using the .C() interface during development)
#ifdef testing

void R_printnninfo(unsigned int n, int* flag, unsigned int* nn_idx,
                   double* nn_diss) {
  unsigned int i, min_idx;
  double min_diss = R_PosInf;
  printf("List of the nearest neighbours\n--------\n");
  for(i = 0; i < n - 1; i++) {
    if(flag[i]) {
      printf("(%d,%d): %f\n", i, nn_idx[i], nn_diss[i]);
      if(nn_diss[i] < min_diss) {
        min_idx = i;
        min_diss = nn_diss[i];
      }
    }
  }
  printf("Minimum: diss(%d and %d) = %f\n", min_idx, nn_idx[min_idx], min_diss);
  return;
}

void R_printdiss(unsigned int n, int* flag, double* diss0) {
  unsigned int i, j, k;
  printf("Similarity matrix\n----------------\n\t");
  for(j = 1; j < n; j++)
    printf("%d\t",j);
  printf("\n");
  for(i = 0, k = 0; i < n - 1; i++) {
    if(flag[i]) {
      printf("%d\t",i);
      for(j = 1; j < n; j++, k++)
        if(flag[j] && (i < j))
          printf("%0.2f\t",diss0[k]);
        else
          printf("----\t");
      printf("\n");
    } else
      k += (n - i - 1);
  }
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
