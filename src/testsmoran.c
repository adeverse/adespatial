#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>


void prodatBc(double *veca, double **matB,double *vecc); 
void testglobal(double *eigenvec, double *eigenval, int *nlig, int *ncol, double *xR, int *nsim, double *sim);
void vecalloc (double **vec, int n);




void freevec (double *vec);
void aleapermutvec (double *a);



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



void prodatBc(double *veca, double **matB,double *vecc) {
/*--------------------------------------------------
* Produit matriciel atB
--------------------------------------------------*/

    int j,  i, lig, col;
    double s;
    
    lig = matB[0][0];
    col = matB[1][0];
 
     for (j=1;j<=col;j++) {
        s = 0;
        for (i=1;i<=lig;i++) {
                 s = s + veca[i] * matB[i][j];
            }
        vecc[j] = s;
     }
}

/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
* Allocation de memoire pour un vecteur de longueur n
--------------------------------------------------*/
{
    if ( (*vec = (double *) calloc(n+1, sizeof(double))) != 0) {
        **vec = n;
        return;
    } else {
        return;
    }
}

/*************************/
void aleapermutvec (double *a)
{
    /* permute au hasard les ?l?ments du vecteur a
    Manly p. 42 Le vecteur est modifi?
    from Knuth 1981 p. 139*/
    int lig, i,j, k;
    double z;
    
    lig = a[0];
    for (i=1; i<=lig-1; i++) {
        j=lig-i+1;
        k = (int) (j*alea()+1);
        /*k = (int) (j*genrand()+1);*/
        if (k>j) k=j;
        z = a[j];
        a[j]=a[k];
        a[k] = z;
    }
}
/***********************************************************************/

/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
* liberation de memoire pour un vecteur
--------------------------------------------------*/
{
    free((char *) vec); 
}

/***********************************************************************/
