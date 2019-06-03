#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "adesub.h"

void testglobal(double *eigenvec, double *eigenval, int *nlig, int *ncol, double *xR, int *nsim, double *sim);



/*--------------------------------------*/

void testglobal(double *eigenvec, double *eigenval, int *nlig, int *ncol, double *xR, int *nsim, double *sim){
  
  int k,nl,nc,i,j;
  double **Evec,*x,*xperm,*cor,*Eval;
  /* 1 I  : somme (lambda)*R2 */
  /* 2 I+ : somme (lambda posi)*R2 */
  /* 3 I- : somme (lambda nega)*R2 */
  nl=*nlig;
  nc=*ncol;
  taballoc(&Evec,nl,nc);
  vecalloc(&Eval,nc);
  vecalloc(&cor,nc);
  vecalloc(&x,nl);
  vecalloc(&xperm,nl);
  k = 0;
  for (i=1; i<=nl; i++) {
    for (j=1; j<=nc; j++) {
      Evec[i][j] = eigenvec[k];
      k = k + 1;
    }
  }
  
  
  for (i=1; i<=nl; i++) {
    x[i]=xR[i-1];
    xperm[i]=x[i];
  }
  
  for (i=1; i<=nc; i++) {
    Eval[i]=eigenval[i-1];
  }
  
  prodatBc(x, Evec,cor);

  for(i=1;i<=nc;i++){
    sim[0]=sim[0]+Eval[i]*cor[i]*cor[i];
    if(Eval[i]>0.0){
      sim[1]=sim[1]+Eval[i]*cor[i]*cor[i];
    }
    if(Eval[i]<0.0){
      sim[2]=sim[2]+Eval[i]*cor[i]*cor[i];
    }
  }
  

  
  for (i=1;i<=*nsim;i++){
    aleapermutvec (xperm);
    prodatBc(xperm, Evec,cor);
    for(k=1;k<=nc;k++){
      sim[3*i]=sim[3*i]+Eval[k]*cor[k]*cor[k];
      if(Eval[k]>0.0){
	sim[3*i+1]=sim[3*i+1]+Eval[k]*cor[k]*cor[k];
      }
      if(Eval[k]<0.0){
	sim[3*i+2]=sim[3*i+2]+Eval[k]*cor[k]*cor[k];
      }
      
    }


    
  }
  

  freevec(Eval);
  freetab(Evec);
  freevec(cor);
  freevec(x);
  freevec(xperm);
}
/*--------------------------------------*/