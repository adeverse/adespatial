#define USE_FC_LEN_T
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include "adesub.h"
#ifndef FCONE
# define FCONE
#endif
/* ============================= */
/*            Declaration        */
/* ============================= */
void forwardsel (double *tabXR, double *tabYR, int *nrowXR, int *ncolXR, int *ncolYR, double *pvalue, int *ordre,
    double *Fvalue, int *nperm, double *R2cum, double *adjR2cum, int *K, double *R2seuil, double *adjR2seuil,double *R2more, int *nbcovar, double *alpha, int *verbose);
void projX (double **tabX, double **projsurX);
void constnewX (double **tabX, double **tabnewX, int *vecsel);
double calcF (double R2X, double R2XZ, int q, int n, int p);
double calcR2 (double **tabY, double **tabYpred);
double calcR2adj (double R2, int n, int nvar);
void tabstandar (double **tab);
double testFreducedmodel(double **predX, double **projectX, double **projectXZ, double **tabres, double **tabXZ, double Fobs, int q, int n, int p, int nperm);
void dinvG(double **X, double **X_inv);


/* ============================= */
/*            Definition         */
/* ============================= */

void forwardsel (double *tabXR, double *tabYR, int *nrowXR, int *ncolXR, int *ncolYR, double *pvalue, int *ordre,
    double *Fvalue, int *nperm, double *R2cum, double *adjR2cum, int *K, double *R2seuil, double *adjR2seuil, double *R2more, int *nbcovar, double *alpha, int *verbose){
    /* CANOCO 4.5 p. 49 */
    /* nouvelle version 01/11/2004 */
    /* on ne stocke plus le projecteur au complet (qui est n x n)
     ce qui permet d'economiser de la RAM quand n est grand 
     modifie aussi testFreducedmodel */
    
    int i,j,k,l,nrowX,ncolX, ncolY, *vecrest, *vecsel, R2maxj=0;
    double **tabX, **tabY, **tabres, **tabpred1, **tabpred2, **tabXnew, R2max=0, adjR2max=0,R2j, R2X, R2XZ, Fobs,**projectX,**projectXZ, **provi;
    nrowX=*nrowXR;
    ncolX=*ncolXR;
    ncolY=*ncolYR;
    
    taballoc (&tabX, nrowX, ncolX);
    taballoc (&tabY, nrowX, ncolY);
    taballoc (&tabres, nrowX, ncolY);
    taballoc (&tabpred1, nrowX, ncolY);
    taballoc (&tabpred2, nrowX, ncolY);
    vecintalloc(&vecrest,ncolX);
    vecintalloc(&vecsel,ncolX);
    taballoc (&projectX, 1,nrowX);
    
    /* Passage des objets R en C */    
    k = 0;
    for (i=1; i<=nrowX; i++) {
        for (j=1; j<=ncolX; j++) {
            tabX[i][j] = tabXR[k];
            k = k + 1;
        }
    }
    
    k = 0;
    for (i=1; i<=nrowX; i++) {
        for (j=1; j<=ncolY; j++) {
            tabY[i][j] = tabYR[k];
            k = k + 1;
            tabpred1[i][j]=0;
        }
    }
    
    for (i=1; i<=ncolX; i++) vecrest[i]=i;
    R2X=0;
    
    for (i=1; i<=ncolX; i++){
        
        if ((R2max>*R2seuil) || (i>*K) || (i>(nrowX-1)) || (adjR2max>*adjR2seuil)) {
            if ((i>*K) )
                Rprintf("Procedure stopped (K criteria)\n");
            if ((i>(nrowX-1)))
                Rprintf("Procedure stopped: number of variables included equals number of rows minus one. All is explained ! Redo your analysis with other parameters.\n");
            break;
        }
        
        if(*verbose == 1){
            Rprintf("Testing variable %d\n",i);
        }
        R2max=0;
        adjR2max=0;
        R2maxj=0;
        taballoc(&tabXnew, nrowX, i);
        taballoc(&projectXZ, i, nrowX);    
        taballoc(&provi,i,ncolY);
        
        for (j=1; j<=ncolX; j++){
            if (vecrest[j]>0){ /* Selection de la variable j base sur le R2 */
                vecsel[i]=j;           
                constnewX(tabX,tabXnew,vecsel);
                projX(tabXnew,projectXZ);
                prodmatABC(projectXZ,tabY,provi);
                prodmatABC(tabXnew,provi,tabpred2);
                R2j=calcR2(tabY,tabpred2);
                if (R2j>R2max) {
                    R2max=R2j;
                    adjR2max=calcR2adj(R2max,nrowX,i);
                    R2maxj=j;
                }
            }
            
        } /* for (j=1; j<=ncolX; j++) */
    
    /* test de la variable j selectionne */
    if ((i>1) && (fabs(R2max-R2cum[i-2])<*R2more)){
        Rprintf("Procedure stopped (R2more criteria): variable %d explains only %f of the variance.\n",i,(fabs(R2max-R2cum[i-2])));
        freetab(tabXnew);
        break;
    }
    
    if ((R2max>*R2seuil)){
        Rprintf("Procedure stopped (R2thresh criteria) R2cum = %f with %d variables (> %f)\n",  R2max, i, *R2seuil);
        break;
        }
    if ((adjR2max>*adjR2seuil)){
        Rprintf("Procedure stopped (adjR2thresh criteria) adjR2cum = %f with %d variables (> %f)\n",  adjR2max, i, *adjR2seuil);
        break;
    }
    
    vecsel[i]=R2maxj;  
    R2cum[i-1]=R2max;
    adjR2cum[i-1]=adjR2max;
    ordre[i-1]=vecsel[i];
    vecrest[R2maxj]=0;
    constnewX(tabX,tabXnew,vecsel);
    projX(tabXnew,projectXZ);
    prodmatABC(projectXZ,tabY,provi);
    prodmatABC(tabXnew,provi,tabpred2);
    for (k=1;k<=nrowX;k++){
        for (l=1;l<=ncolY;l++){
            tabres[k][l]=tabY[k][l]-tabpred1[k][l];
        }
    }
    R2XZ=R2max;    
    Fobs=calcF(R2X,R2XZ,1,nrowX,i+(*nbcovar));
    Fvalue[i-1]=Fobs;
    pvalue[i-1]=testFreducedmodel(tabpred1,projectX,projectXZ,tabres, tabXnew, Fobs,1,nrowX,i+(*nbcovar),*nperm);
    for (k=1;k<=nrowX;k++){
        for (l=1;l<=ncolY;l++){
            tabpred1[k][l]=tabpred2[k][l];/* mettre dans pred1 les predictions du nouvo model */
    
        }
    }
    freetab(projectX);
    taballoc(&projectX, i, nrowX); 
    for (k=1;k<=i;k++){
        for (l=1;l<=nrowX;l++){
            projectX[k][l]=projectXZ[k][l];
            
        }
    }
    
    R2X=R2XZ;
    freetab(tabXnew);
    freetab(projectXZ);
    freetab(provi);
    
    if ((pvalue[i-1]>*alpha)) {
        Rprintf("Procedure stopped (alpha criteria): pvalue for variable %d is %f (> %f)\n",i,pvalue[i-1],*alpha);
        ordre[i-1]=0;
        break;
    }
    
    
    }/* for (i=1; i<=ncolX; i++) */
    
    freetab(projectX);
    freeintvec(vecsel);
    freeintvec(vecrest);
    freetab(tabX);
    freetab(tabY);
    freetab(tabres);
    freetab(tabpred1);
    freetab(tabpred2);
}

/* ============================= */
/* ============================= */

double testFreducedmodel(double **predX, double **projectX, double **projectXZ, double **tabres, double **tabXZ, double Fobs, int q, int n, int p, int nperm){
    int i,j,k,nrowY,ncolY,NGT=0, ncolX, ncolXZ;
    double **Yperm,R2X=0,R2XZ, Fi,**tabpredX,**tabpredXZ, **proviX, **proviXZ, **tabX;
    nrowY=predX[0][0];
    ncolY=predX[1][0];
    ncolXZ=tabXZ[1][0];
    ncolX=ncolXZ-1;
    taballoc(&Yperm,nrowY,ncolY);
    taballoc(&tabpredXZ,nrowY,ncolY);
    taballoc(&proviXZ,ncolXZ,ncolY);
    
    if(ncolX>0){
        taballoc(&tabpredX,nrowY,ncolY);
        taballoc(&proviX,ncolX,ncolY);
        taballoc(&tabX,nrowY, ncolX);
        for (i=1;i<=nrowY;i++){
            for (j=1;j<=ncolX;j++){
                tabX[i][j]=tabXZ[i][j];
            }
        }
    }
    
    for (i=1;i<=nperm;i++){
        aleapermutmat (tabres);
        for (k=1;k<=nrowY;k++){
            for (j=1;j<=ncolY;j++){
                Yperm[k][j]=predX[k][j]+tabres[k][j];
            }
        }
        
        prodmatABC(projectXZ,Yperm,proviXZ);
        prodmatABC(tabXZ,proviXZ,tabpredXZ);
        if(ncolX>0){
            prodmatABC(projectX,Yperm,proviX);
            prodmatABC(tabX,proviX,tabpredX);
        }
        if(ncolX>0){
            R2X=calcR2(Yperm,tabpredX);
        }   
        R2XZ=calcR2(Yperm,tabpredXZ);    
        Fi=calcF(R2X,R2XZ,q,n,p);
        if (Fi>=Fobs) NGT=NGT+1;
        
    }
    freetab(Yperm);
    freetab(tabpredXZ);
    freetab(proviXZ);
    if(ncolX>0){
        freetab(tabpredX);
        freetab(tabX);
        freetab(proviX);
    }
    return(((double)NGT+1)/((double)nperm+1));
    
    
}
/* ============================= */
/* ============================= */

double calcR2 (double **tabY, double **tabYpred){
    int i,j,nrowY,ncolY;
    double res=0;
    nrowY=tabY[0][0];
    ncolY=tabY[1][0];
    tabstandar(tabY);
    tabstandar(tabYpred);
    for (j=1;j<=ncolY;j++){
        for (i=1;i<=nrowY;i++){
            res=res+tabY[i][j]*tabYpred[i][j];
        }
    }
    res=pow(res/(double)(nrowY*ncolY),2);
    return(res);
}

/* ============================= */
/* ============================= */

double calcR2adj (double R2, int n, int nvar){
    double res=0;
    res =1-(1-R2)*( (double)(n - 1) )/( (double)(n - nvar -1));
    return(res);
}

/* ============================= */
/* ============================= */

void tabstandar (double **tab)
    /*--------------------------------------------------
     * tab est un tableau                               
     * la procedure retourne tab norme au total 
     * variance en n 
     --------------------------------------------------*/
{
    double      z, v2,x;
    int         j,i,l1,c1;
    double      moy=0, var=0;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    
    /*--------------------------------------------------
     * calcul du tableau centre/norme
     --------------------------------------------------*/
    moy=0;
    var=0;
    for (j=1;j<=c1;j++){
        
        for (i=1;i<=l1;i++) {
            moy = moy + tab[i][j];
        }
    }   
    moy=moy/(double)(l1*c1); 
    for (j=1;j<=c1;j++){
        for (i=1;i<=l1;i++) {
            x = tab[i][j] - moy;
            var = var +  x * x /(double)(l1*c1);
        }
    }
    
    
    v2 = var;
    if (v2<=0) v2 = 1;
    v2 = sqrt(v2);
    var = v2;
    
    for (j=1;j<=c1;j++){ 
        for (i=1;i<=l1;i++) {
            z = (tab[i][j] - moy)/var;
            tab[i][j] = z;
        }
    }
    
}

/*=========================================================================*/
/* ============================= */
/* ============================= */

void constnewX (double **tabX, double **tabnewX, int *vecsel){
    int nrowX, ncolnewX,i,j;
    nrowX=tabX[0][0];
    ncolnewX=tabnewX[1][0];
    for (i=1;i<=nrowX;i++){
        for (j=1;j<=ncolnewX;j++){
            tabnewX[i][j]=tabX[i][(vecsel[j])];
        }
    }
    
}

/* ============================= */
/* ============================= */
void projX (double **tabX, double **projsurX){
    int nlX,ncX,i,j;
    double **provi1, **provi2, **provi3;
    nlX=tabX[0][0];
    ncX=tabX[1][0];
    
    taballoc(&provi1,ncX,ncX);
    prodmatAtAB (tabX, provi1); /* provi1=XtX */
taballoc(&provi2,ncX,ncX);
dinvG(provi1,provi2); /* provi2=(XtX)-1 */
freetab(provi1);
/*taballoc(&provi1,nlX,ncX);
 prodmatABC(tabX,provi2,provi1);  provi1=X(XtX)-1 
 */

taballoc(&provi3,ncX,nlX);
for (i=1; i<=nlX;i++){
    for (j=1; j<=ncX;j++){
        provi3[j][i]=tabX[i][j]; /* provi3=Xt */
    }
}

prodmatABC(provi2,provi3,projsurX); /* projsurX=(XtX)-1Xt */
freetab(provi2);
freetab(provi3);

}


/* ============================= */
/* ============================= */

double calcF (double R2X, double R2XZ, int q, int n, int p){
    double res=0;
    res=(R2XZ-R2X)/(double)q;
    res=res/((1-R2XZ)/(double)(n-p-q));
    return(res);
    
}

/*============================== */
/*============================== */

/* inverse generalise par svd d'une matrice symetrique */
/*DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 $                   WORK, LWORK, INFO ) */
void dinvG(double **X, double **X_inv)
{
    int i,j, k, l,error,size,lwork;
    size=X[1][0];
    double *A = (double *)calloc((size_t)size*size, sizeof(double));/*doubleArray(size*size);*/
double *D = (double *)calloc((size_t)size, sizeof(double));/*doubleArray(size*1);*/
double *U = (double *)calloc((size_t)size*size, sizeof(double));
double *V = NULL,work1,*work;
double **XU, **XUred;
const char jobu='A',jobvt='N';

taballoc(&XU,size,size);
lwork=-1; 
for (i = 0, j = 1; j <= size; j++) {
    for (k = 1; k <= size; k++) {
        A[i] = X[k][j];
        i++;
    }
}
F77_CALL(dgesvd)(&jobu, &jobvt,&size, &size,A, &size, D,U,&size,V,&size,&work1, &lwork,&error FCONE FCONE);

lwork=(int)floor(work1);
if (work1-lwork>0.5) lwork++;
work=(double *)calloc((size_t)lwork,sizeof(double));
/* actual call */
F77_NAME(dgesvd)(&jobu, &jobvt,&size, &size,A, &size, D,U,&size,V,&size,work, &lwork,&error FCONE FCONE);
free(work);

if (error) {
    Rprintf("error in svd: %d\n", error);
}
i = 0;
l=0;
for ( j = 1; j <= size; j++) {
    for (k = 1; k <= size; k++) {
        XU[k][j] = U[i];
        i++;
    }
    
    if (D[j-1]>0.00000000001) l=l+1;
    
}

taballoc(&XUred,size,l);
for (i=1;i<=size;i++) {
    for (k=1;k<=l;k++) {
        XUred[i][k]=pow(D[k-1],-0.5)*XU[i][k];
    }
}     

prodmatAAtB (XUred, X_inv);

freetab(XUred);
free(A);
free(D);
free(U);
freetab(XU);
}



/* ============================= */
/* ============================= */

