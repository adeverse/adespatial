#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <R_ext/Applic.h> /* for dgemm */
#ifndef FCONE
# define FCONE
#endif
/***********sampleIntC : permutes a vector of integers**********************/

SEXP sampleIntC (SEXP x)
{
    x = PROTECT(coerceVector(x, INTSXP));
    
    SEXP res;
    R_len_t i, irdm;
    double rdm;
    int tmp;
    int I = length(x);
    
    PROTECT(res = allocVector(INTSXP,I));
    memset(INTEGER(res),0,I*sizeof(int));
    for(i=0;i<I;i++)  INTEGER(res)[i]=INTEGER(x)[i];
    
    GetRNGstate();
    for (i = 0 ; i < I ; i++)
    {    
        do {
            rdm = unif_rand();                 /* Sytem RNG  */
        } while (rdm == 1.0);                  /* pour ne pas avoir un indice irdm  = I */  
irdm = (R_len_t) (rdm * I); 
tmp = INTEGER(res)[irdm];
INTEGER(res)[irdm]  = INTEGER(res)[i];
INTEGER(res)[i] = tmp;  
    } 
    PutRNGstate();   
    
    UNPROTECT(2);
    return(res);
}
/*********** End  sampleIntC **********************/

/*********** restricted_perm permutes a vector of size n 
 in 3 different ways, nobs_bloc is the number of observations per 
 block (i.e., the number of sites), nbloc is the number of 
 time blocks **********************/

SEXP RestrictedPerm(SEXP nobs_bloc,SEXP nbloc,SEXP n,SEXP restricted_perm)
{
    SEXP vect_perm,vect_perm2,vect_perm3,vect_perm4,vect_perm_total,vect;
    R_len_t ind=0,j,k;
    
    nobs_bloc = PROTECT(coerceVector(nobs_bloc, INTSXP));
    nbloc = PROTECT(coerceVector(nbloc, INTSXP));
    n = PROTECT(coerceVector(n, INTSXP));
    restricted_perm = PROTECT(coerceVector(restricted_perm, INTSXP));
    
    int n_bloc=INTEGER(nbloc)[0],n_obs=INTEGER(nobs_bloc)[0];
    
    PROTECT(vect_perm = allocVector(INTSXP,n_obs));
    memset(INTEGER(vect_perm),0,n_obs*sizeof(int));  
    
    PROTECT(vect_perm4 = allocVector(INTSXP,n_obs));
    memset(INTEGER(vect_perm4),0,n_obs*sizeof(int)); 
    
    PROTECT(vect_perm2 = allocVector(INTSXP,n_bloc));
    memset(INTEGER(vect_perm2),0,n_bloc*sizeof(int));  
    
    PROTECT(vect_perm3 = allocVector(INTSXP,n_bloc));
    memset(INTEGER(vect_perm3),0,n_bloc*sizeof(int));  
    
    
    PROTECT(vect_perm_total = allocVector(INTSXP,INTEGER(n)[0]));
    memset(INTEGER(vect_perm_total),0,INTEGER(n)[0]*sizeof(int));  
    
    PROTECT(vect= allocVector(INTSXP,INTEGER(n)[0]));
    memset(INTEGER(vect),0,INTEGER(n)[0]*sizeof(int)); 
    
    for(j=0;j<INTEGER(n)[0];j++) INTEGER(vect)[j]=j;
    
    if (INTEGER(restricted_perm)[0] == 0)     vect_perm_total = sampleIntC(vect);
    
    else if (INTEGER(restricted_perm)[0] == 1) /* restricted_perm=1 space permuted and time fixed */
    {
        for (j=0;j<n_bloc;j++)
        {
            ind=n_obs*j;
            for (k=0;k<n_obs;k++)  INTEGER(vect_perm)[k]=INTEGER(vect)[ind+k];
            vect_perm4 = sampleIntC(vect_perm);
            for(k=0;k<n_obs;k++)  INTEGER(vect_perm_total)[j*n_obs+k]=INTEGER(vect_perm4)[k];
        }
    } 
    else  /* restricted _perm = 2 space fixed and time permuted */
    {  
        for (k=0;k<n_obs;k++) 
        {
            for (j=0;j<n_bloc;j++)  INTEGER(vect_perm2)[j]=INTEGER(vect)[k+n_obs*j];
            vect_perm3=sampleIntC(vect_perm2);
            for (j=0;j<n_bloc;j++) INTEGER(vect_perm_total)[k+n_obs*j]=INTEGER(vect_perm3)[j];
        }  
    }  
    UNPROTECT(10);
    return(vect_perm_total);
}

/*********** End restricted_perm**********************/

/* reorder_mat reorders rows of matrix x following order of vector vec  */

SEXP reorder_mat(SEXP x, SEXP vect)
{
    SEXP Rdim,xperm;
    R_len_t  ncol,nline,i,j,l;
    l = length(vect);   //vect=c(0,..,n-1) because indices begin with 0 in C
    x = PROTECT(coerceVector(x, REALSXP)); 
    vect = PROTECT(coerceVector(vect, INTSXP));
    Rdim = PROTECT(getAttrib(x, R_DimSymbol)); 
    ncol = INTEGER(Rdim)[1];
    nline = INTEGER(Rdim)[0];
    PROTECT(xperm = allocMatrix(REALSXP, nline, ncol));
    memset(REAL(xperm),0.0,nline*ncol*sizeof(double)); 
    
    for (i=0;i<l;i++) 
    {
        for(j=0;j<ncol;j++) 
        {
            int indice=INTEGER(vect)[i];
            REAL(xperm)[j*nline+i]=REAL(x)[j*nline+indice];  
        }
    }
    
    UNPROTECT(4);
    return(xperm);
}
/***** End reorder_mat ****/

/**** produit_dgemm computes matrix product  *****/

SEXP produit_dgemm(SEXP X, SEXP Y) {
    
    SEXP Rdim1,Rdim2,res;
    X = PROTECT(coerceVector(X, REALSXP)); 
    Y = PROTECT(coerceVector(Y, REALSXP));
    
    Rdim1 = PROTECT(getAttrib(X, R_DimSymbol)); 
    Rdim2 = PROTECT(getAttrib(Y, R_DimSymbol)); 
    
    double *xptr; 
    xptr = REAL(X);
    double *yptr; 
    yptr = REAL(Y);
    
    int *dimX; 
    dimX = INTEGER(Rdim1); 
    
    int *dimY; 
    dimY = INTEGER(Rdim2);
    
    PROTECT(res = allocMatrix(REALSXP, dimX[0], dimY[1]));
    
    double *resptr; 
    resptr = REAL(res);
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    
    F77_CALL(dgemm)(transa, transb, &dimX[0], &dimY[1], &dimX[1], &one, xptr, &dimX[0], yptr, &dimY[0], &zero, resptr, &dimX[0] FCONE FCONE); 
    //F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,x, &nrx, y, &nry, &zero, z, &nrx);
    
    UNPROTECT(5);
    return(res); 
}
/*** End produit_dgemm  *****/

/**** SS computes sums of squares  ***/
SEXP SS(SEXP x)
{
    SEXP result, Rdim;
    R_len_t I,J,i;
    x = PROTECT(coerceVector(x, REALSXP)); 
    Rdim = PROTECT(getAttrib(x, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    PROTECT(result = allocVector(REALSXP,1));
    memset(REAL(result),0.0,sizeof(double)); 
    
    for (i=0;i<I*J;i++) 
    {
        REAL(result)[0] += pow(REAL(x)[i],2); 
    }
    UNPROTECT(3);
    return(result);
}
/***  End SS ***/

/*** sti_loop computes permutation test for space-time interaction ***/
SEXP sti_loop(SEXP nperm,SEXP Y,SEXP s,SEXP tt, SEXP a, SEXP b, SEXP cc, SEXP SS_Y, SEXP Fref_AxB, SEXP proj_AxB, SEXP proj_ABAxB)
{
    R_len_t nline,i,iperm;
    
    SEXP  Rdim,vect,nPGE_AxB;
    double MS_Perm_Res=0.0,MS_Perm_AxB=0.0, Fper_AxB=0.0;
    
    nperm = PROTECT(coerceVector(nperm, INTSXP));
    s = PROTECT(coerceVector(s, INTSXP));
    tt = PROTECT(coerceVector(tt, INTSXP));
    a = PROTECT(coerceVector(a, INTSXP));
    b = PROTECT(coerceVector(b, INTSXP));
    cc = PROTECT(coerceVector(cc, INTSXP));
    SS_Y = PROTECT(coerceVector(SS_Y, REALSXP));
    Fref_AxB = PROTECT(coerceVector(Fref_AxB, REALSXP));
    proj_AxB = PROTECT(coerceVector(proj_AxB, REALSXP));
    proj_ABAxB = PROTECT(coerceVector(proj_ABAxB, REALSXP));
    Y = PROTECT(coerceVector(Y, REALSXP)); 
    
    Rdim = PROTECT(getAttrib(Y, R_DimSymbol)); 
    nline = INTEGER(Rdim)[0];
    
    
    PROTECT(nPGE_AxB = allocVector(INTSXP,1));
    memset(INTEGER(nPGE_AxB),0,sizeof(int)); 
    
    
    PROTECT(vect = allocVector(INTSXP,nline));
    memset(INTEGER(vect),0,nline*sizeof(int)); 
    
    for(i=0;i<nline;i++) INTEGER(vect)[i]=i;
    
    INTEGER(nPGE_AxB)[0]=1;
    
    SEXP SEXP_nline= PROTECT(ScalarInteger(nline));
    SEXP SEXP_zero=PROTECT(ScalarInteger(0));
    
    for(iperm=1;iperm<=INTEGER(nperm)[0];iperm++)
    { 
        /* PROTECT ..UNPROTECT added on Feb 15 2019 NM*/
        SEXP YRestricted=PROTECT(RestrictedPerm(s,tt,SEXP_nline,SEXP_zero)); 
        SEXP Yperm=PROTECT(reorder_mat(Y,YRestricted));
        
        SEXP YhatPerm_AxB=PROTECT(produit_dgemm(proj_AxB,Yperm));  
        SEXP SS_YhatPerm_AxB=PROTECT(SS(YhatPerm_AxB));
        MS_Perm_AxB=REAL(SS_YhatPerm_AxB)[0]/INTEGER(cc)[0];
        
        SEXP YhatPerm_ABAxB=PROTECT(produit_dgemm(proj_ABAxB,Yperm));
        SEXP SS_YhatPerm_ABAxB = PROTECT(SS(YhatPerm_ABAxB)); 
        MS_Perm_Res=(REAL(SS_Y)[0]-REAL(SS_YhatPerm_ABAxB)[0])/(nline-(INTEGER(a)[0]+INTEGER(b)[0]+INTEGER(cc)[0])-1);
        
        Fper_AxB=MS_Perm_AxB/MS_Perm_Res;
        
        if(Fper_AxB >= REAL(Fref_AxB)[0])  INTEGER(nPGE_AxB)[0] += 1; 
        
        UNPROTECT(6);
    }
    
    /*UNPROTECT(20);*/
    UNPROTECT(16);
    return(nPGE_AxB);
}
/******* End sti_loop  ****/


/*** s_loop computes permutation test for space ***/
SEXP s_loop(SEXP nbperm,SEXP Y,SEXP s,SEXP tt, SEXP a, SEXP b, SEXP cc, SEXP SS_Y, SEXP Fref_A, SEXP projA, SEXP projB, SEXP projAXB, SEXP projABAxB, SEXP T_fixed)
{
    R_len_t nline,iperm;
    
    SEXP Rdim, nPGE_A;
    double MS_Perm_Res=0.0,MS_Perm_AxB=0., MS_Perm_A=0.0, MS_Perm_B=0.0, Fper_A=0.0;
    
    nbperm = PROTECT(coerceVector(nbperm, INTSXP));
    s = PROTECT(coerceVector(s, INTSXP));
    tt = PROTECT(coerceVector(tt, INTSXP));
    a = PROTECT(coerceVector(a, INTSXP));
    b = PROTECT(coerceVector(b, INTSXP));
    cc = PROTECT(coerceVector(cc, INTSXP));
    T_fixed = PROTECT(coerceVector(T_fixed, INTSXP));
    SS_Y = PROTECT(coerceVector(SS_Y, REALSXP));
    Fref_A = PROTECT(coerceVector(Fref_A, REALSXP));
    projAXB = PROTECT(coerceVector(projAXB, REALSXP));
    projABAxB = PROTECT(coerceVector(projABAxB, REALSXP));
    projA = PROTECT(coerceVector(projA, REALSXP));
    projB = PROTECT(coerceVector(projB, REALSXP));
    Y = PROTECT(coerceVector(Y, REALSXP)); 
    
    Rdim = PROTECT(getAttrib(Y, R_DimSymbol)); 
    nline = INTEGER(Rdim)[0];
    
    PROTECT(nPGE_A = allocVector(INTSXP,1));
    memset(INTEGER(nPGE_A),0,sizeof(int)); 
    
    INTEGER(nPGE_A)[0]=1;
    
    SEXP SEXP_nline= PROTECT(ScalarInteger(nline));
    SEXP SEXP_un=PROTECT(ScalarInteger(1));
    
    for(iperm=1;iperm<=INTEGER(nbperm)[0];iperm++)
    {
        /* PROTECT ..UNPROTECT added on Feb 15 2019 NM */
        SEXP YRestricted = PROTECT(RestrictedPerm(s,tt,SEXP_nline,SEXP_un));	    
        SEXP Yperm=PROTECT(reorder_mat(Y,YRestricted));	    
        SEXP YhatPerm_A=PROTECT(produit_dgemm(projA,Yperm));     
        SEXP SS_YhatPerm_A = PROTECT(SS(YhatPerm_A)); 
        MS_Perm_A=REAL(SS_YhatPerm_A)[0]/INTEGER(a)[0];
        
        if(INTEGER(T_fixed)[0]==1) { // Time random factor in crossed design with interaction
            SEXP YhatPerm_AxB=PROTECT(produit_dgemm(projAXB,Yperm));
            SEXP SS_YhatPerm_AxB=PROTECT(SS(YhatPerm_AxB));
            MS_Perm_AxB=REAL(SS_YhatPerm_AxB)[0]/INTEGER(cc)[0];
            Fper_A=MS_Perm_A/MS_Perm_AxB;
            UNPROTECT(2);
            
        } else if(INTEGER(T_fixed)[0]==2) { // Time random factor in nested design
            SEXP YhatPerm_B=PROTECT(produit_dgemm(projB,Yperm));
            SEXP SS_YhatPerm_B = PROTECT(SS(YhatPerm_B)); 	
            MS_Perm_B = REAL(SS_YhatPerm_B)[0]/INTEGER(b)[0];
            Fper_A = MS_Perm_A/MS_Perm_B;
            UNPROTECT(2);
            
        } else {	
            SEXP YhatPerm_ABAxB=PROTECT(produit_dgemm(projABAxB,Yperm)); 
            SEXP SS_YhatPerm_ABAxB =PROTECT(SS(YhatPerm_ABAxB)); 
            MS_Perm_Res= (REAL(SS_Y)[0]-REAL(SS_YhatPerm_ABAxB)[0])/(nline-(INTEGER(a)[0]+INTEGER(b)[0]+INTEGER(cc)[0])-1);
            Fper_A=MS_Perm_A/MS_Perm_Res;
            UNPROTECT(2);					
        }
        
        if(Fper_A >= REAL(Fref_A)[0])    INTEGER(nPGE_A)[0]+=1;
        /* added on Feb 15 2019 NM */
        UNPROTECT(4);
    }
    /*UNPROTECT(26);*/
    UNPROTECT(18);
    return(nPGE_A);
    
}
/***** End s_loop  *****/

/*** t_loop computes permutation test for time ***/
SEXP t_loop(SEXP nperm,SEXP Y,SEXP s,SEXP tt, SEXP a, SEXP b, SEXP cc, SEXP SS_Y, SEXP Fref_B, SEXP projA, SEXP projB, SEXP projAXB, SEXP projABAxB, SEXP T_fixed)
{
    R_len_t nline,iperm; 
    
    SEXP  Rdim, nPGE_B;
    double MS_Perm_Res=0.0,MS_Perm_AxB=0., MS_Perm_A=0.0, MS_Perm_B=0.0, Fper_B=0.0;
    
    nperm = PROTECT(coerceVector(nperm, INTSXP));
    s = PROTECT(coerceVector(s, INTSXP));
    tt = PROTECT(coerceVector(tt, INTSXP));
    a = PROTECT(coerceVector(a, INTSXP));
    b = PROTECT(coerceVector(b, INTSXP));
    cc = PROTECT(coerceVector(cc, INTSXP));
    T_fixed = PROTECT(coerceVector(T_fixed, INTSXP));
    SS_Y = PROTECT(coerceVector(SS_Y, REALSXP));
    Fref_B = PROTECT(coerceVector(Fref_B , REALSXP));
    projAXB = PROTECT(coerceVector(projAXB, REALSXP));
    projABAxB = PROTECT(coerceVector(projABAxB, REALSXP));
    projA = PROTECT(coerceVector(projA, REALSXP));
    projB = PROTECT(coerceVector(projB, REALSXP));
    Y = PROTECT(coerceVector(Y, REALSXP)); 
    
    Rdim = PROTECT(getAttrib(Y, R_DimSymbol)); 
    nline = INTEGER(Rdim)[0];   
    
    PROTECT(nPGE_B = allocVector(INTSXP,1));
    memset(INTEGER(nPGE_B),0,sizeof(int)); 
    INTEGER(nPGE_B)[0]=1;
    
    SEXP SEXP_nline= PROTECT(ScalarInteger(nline));
    SEXP SEXP_deux=PROTECT(ScalarInteger(2));
    
    for(iperm=1;iperm<=INTEGER(nperm)[0];iperm++)
    {
        /* PROTECT ..UNPROTECT added on Feb 15 2019 NM */ 
        SEXP YRestricted = PROTECT(RestrictedPerm(s,tt,SEXP_nline,SEXP_deux)); 
        SEXP Yperm=PROTECT(reorder_mat(Y,YRestricted));   
        SEXP YhatPerm_B=PROTECT(produit_dgemm(projB,Yperm));
        SEXP SS_YhatPerm_B=PROTECT(SS(YhatPerm_B));
        MS_Perm_B = REAL(SS_YhatPerm_B)[0]/INTEGER(b)[0];
        
        
        if(INTEGER(T_fixed)[0]==1) { // Time random factor in crossed design with interaction
            SEXP YhatPerm_AxB=PROTECT(produit_dgemm(projAXB,Yperm));
            SEXP SS_YhatPerm_AxB=PROTECT(SS(YhatPerm_AxB));
            MS_Perm_AxB=REAL(SS_YhatPerm_AxB)[0]/INTEGER(cc)[0];  
            Fper_B=MS_Perm_B/MS_Perm_AxB;
            UNPROTECT(2);
            
        } else if(INTEGER(T_fixed)[0]==2) { // Time random factor in nested design
            SEXP YhatPerm_A=PROTECT(produit_dgemm(projA,Yperm));
            SEXP SS_YhatPerm_A=PROTECT(SS(YhatPerm_A));
            MS_Perm_A=REAL(SS_YhatPerm_A)[0]/INTEGER(a)[0];
            Fper_B=MS_Perm_B/MS_Perm_A;
            UNPROTECT(2);
            
        } else {
            SEXP YhatPerm_ABAxB=PROTECT(produit_dgemm(projABAxB,Yperm)); 	
            SEXP SS_YhatPerm_ABAxB=PROTECT(SS(YhatPerm_ABAxB));
            MS_Perm_Res= (REAL(SS_Y)[0]-REAL(SS_YhatPerm_ABAxB)[0])/(nline-(INTEGER(a)[0]+INTEGER(b)[0]+INTEGER(cc)[0])-1);
            Fper_B=MS_Perm_B/MS_Perm_Res;
            UNPROTECT(2);
        }	
        if(Fper_B >= REAL(Fref_B)[0])    INTEGER(nPGE_B)[0]+=1;
        UNPROTECT(4);		
    }
    
    /*UNPROTECT(26);*/
    UNPROTECT(18);
    return(nPGE_B);
}

/***** End t_loop  *****/
