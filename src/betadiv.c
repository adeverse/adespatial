/*****************************************************            betadiv.c                *****************************************************/
/**************************************************** Naima MADI DESS. BioInformatique UQAM ****************************************************/
/*****************************************************          27 septembre 2016          *****************************************************/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <float.h>
#include <math.h>
#include <R_ext/Applic.h> /* for dgemm */

#define FOR_RAND 1/RAND_MAX

/** debut calcul pour betadiv1   : Distance euclédienne et transformations des matrices de données **/

/// Euclidean Distance ========
SEXP euclidean (SEXP RinMatrix) {
    R_len_t i,j,I,J,k;
    SEXP Rdim, Rval;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), somme = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(Rval = allocMatrix(REALSXP, I, I));
    memset(REAL(Rval),0.0,I*I*sizeof(double)); 
    double *Rv = REAL(Rval);
    
    for (i = 0; i < I; i++) 
        for (j = i+1; j < I; j++) {
            for (k = 0; k < J; k++) somme += (rx[k*I+i]-rx[I*k+j])*(rx[k*I+i]-rx[I*k+j]);
            Rv[i*I+j] = sqrt(somme);
            somme = 0.0;
        }
        UNPROTECT(3);
    return(Rval);
}

///// chi-square ========
SEXP chisquare(SEXP RinMatrix) {
    R_len_t i,j,I,J,k;
    SEXP Rdim, Rval;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), sum_line = 0.0, sum_total = 0.0, sum_col = 0.0, tmp=0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(Rval = allocMatrix(REALSXP, I, J));
    memset(REAL(Rval),0.0,I*J*sizeof(double)); 
    double *Rv = REAL(Rval);
    
    for ( i=0;i<I*J;i++) sum_total += rx[i]; 
    for (i = 0; i < I; i++) 
    {
        for (k=0; k<J; k++) sum_line += rx[k*I+i];
        for (j = 0; j < J; j++)   
        { 
            for (k=0; k<I; k++)       sum_col += rx[j*I+k]; 
            if(sum_line<DBL_EPSILON)  sum_line = DBL_EPSILON;
            tmp=sqrt(sum_col);
            if(tmp<DBL_EPSILON)      tmp = DBL_EPSILON;
            Rv[j*I+i] = sqrt(sum_total)*rx[j*I+i]/(sum_line*tmp);
            sum_col=0.0;
        }
        sum_line = 0.0;
    }
    UNPROTECT(3);
    return(Rval);
}

///// profiles ========
SEXP profiles(SEXP RinMatrix) {
    R_len_t i,j,I,J,k;
    SEXP Rdim, Rval;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), sum_line = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(Rval = allocMatrix(REALSXP, I, J));
    memset(REAL(Rval),0.0,I*J*sizeof(double)); 
    double *Rv = REAL(Rval);
    
    for (i = 0; i < I; i++) 
    {
        for (k=0; k<J; k++)        sum_line += rx[k*I+i];
        if(sum_line<DBL_EPSILON)   sum_line = DBL_EPSILON;
        for (j = 0; j < J; j++)    Rv[j*I+i] =  rx[j*I+i]/sum_line;
        sum_line = 0.0;  
    }
    UNPROTECT(3);
    return(Rval);
}

///// chord ========
SEXP chord(SEXP RinMatrix) {
    R_len_t i,j,I,J,k;
    SEXP Rdim, Rval;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), sum_line = 0.0,tmp=0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(Rval = allocMatrix(REALSXP, I, J));
    memset(REAL(Rval),0.0,I*J*sizeof(double)); 
    double *Rv = REAL(Rval);
    
    for (i = 0; i < I; i++) 
    {
        for (k=0; k<J; k++)    sum_line += pow(rx[k*I+i],2);
        tmp=sqrt(sum_line);
        if(tmp<DBL_EPSILON)    tmp = DBL_EPSILON;
        for (j=0;j<J;j++)      Rv[j*I+i] = rx[j*I+i]/tmp;
        sum_line = 0.0;
    }
    UNPROTECT(3);
    return(Rval);
}

///// Hellinger ========
SEXP hellinger(SEXP RinMatrix) {
    R_len_t i,j,I,J,k;
    SEXP Rdim, Rval;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), sum_line = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(Rval = allocMatrix(REALSXP, I, J));
    memset(REAL(Rval),0.0,I*J*sizeof(double)); 
    
    double *Rv = REAL(Rval);
    
    for (i = 0; i < I; i++) 
    {
        for (k=0; k<J; k++)       sum_line += rx[k*I+i];
        if(sum_line<DBL_EPSILON)  sum_line = DBL_EPSILON;
        for (j = 0; j < J; j++)   Rv[j*I+i]= sqrt(rx[j*I+i]/sum_line);
        sum_line = 0.0;  
    }
    UNPROTECT(3);
    return(Rval);
}

//// transformation de la matrice de données ========
SEXP transform_mat(SEXP RinMatrix, SEXP type_transf) {
    SEXP Rval, Rdim;
    R_len_t I,J;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(Rval = allocMatrix(REALSXP, I, J));
    memset(REAL(Rval),0.0,I*J*sizeof(double)); 
    
    PROTECT(type_transf = AS_CHARACTER(type_transf));
    
    if (STRING_ELT(type_transf,0) == mkChar("profiles"))  Rval = profiles(RinMatrix);
    else if (STRING_ELT(type_transf,0) == mkChar("chisquare")) Rval = chisquare(RinMatrix); 
    else if (STRING_ELT(type_transf,0) == mkChar("chord"))     Rval = chord(RinMatrix); 
    else if (STRING_ELT(type_transf,0) == mkChar("hellinger")) Rval = hellinger(RinMatrix); 
    else if (STRING_ELT(type_transf,0) == mkChar("euclidean")) Rval = RinMatrix; 
    
    UNPROTECT(4);
    return Rval;
}
/* fin transformation */

////  Calculs pour LCBD ========
SEXP squared_diff(SEXP RinMatrix) {
    R_len_t i,j,I,J,k;
    SEXP Rdim, Rval;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), sum_col = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(Rval = allocMatrix(REALSXP, I, J));
    memset(REAL(Rval),0.0,I*J*sizeof(double)); 
    double *Rv = REAL(Rval);
    
    for (j = 0; j < J; j++)   
    { 
        for (k=0; k<I; k++)      sum_col += rx[j*I+k]; 
        for (i=0; i<I; i++)      Rv[j*I+i] = pow(rx[j*I+i]-(sum_col/I),2);
        sum_col = 0.0;
    }
    UNPROTECT(3);
    return(Rval);
}


SEXP calcul_BD(SEXP RinMatrix) {
    SEXP Rdim;
    R_len_t I,J,k, i, j;
    double sum_line = 0.0, sum_col = 0.0;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    SEXP SStotal;
    PROTECT(SStotal = allocVector(REALSXP,1));
    memset(REAL(SStotal), 0, sizeof(double));
    
    SEXP Rval_lcbd;
    PROTECT(Rval_lcbd = allocVector(REALSXP,I));
    memset(REAL(Rval_lcbd), 0.0 , sizeof(double)*I);
    
    SEXP Rval_scbd;
    PROTECT(Rval_scbd = allocVector(REALSXP,J));
    memset(REAL(Rval_scbd), 0.0 , sizeof(double)*J);
    
    SEXP BDtotal;
    PROTECT(BDtotal = allocVector(REALSXP,1));
    memset(REAL(BDtotal), 0.0, sizeof(double));
    
    SEXP Rlist = PROTECT(allocVector(VECSXP,4));
    
    // SSTOTAL et BDTOTAL
    for(k=0; k<I*J; k++) REAL(SStotal)[0] += REAL(RinMatrix)[k];
    REAL(BDtotal)[0] = REAL(SStotal)[0]/ ((double)(I) - 1);
    SET_VECTOR_ELT(Rlist,0,SStotal);
    SET_VECTOR_ELT(Rlist,1,BDtotal);
    
    // LCBD  
    if(REAL(SStotal)[0]<DBL_EPSILON) REAL(SStotal)[0]=DBL_EPSILON;
    for (i = 0; i < I; i++) 
    {
        for (k=0; k<J; k++) sum_line += REAL(RinMatrix)[k*I+i];
        REAL(Rval_lcbd)[i] = sum_line/REAL(SStotal)[0];
        sum_line = 0.0; 
    }
    SET_VECTOR_ELT(Rlist,2,Rval_lcbd);
    
    // SCBD
    for (j = 0; j < J; j++)
    {
        for (k=0; k<I; k++) sum_col += REAL(RinMatrix)[j*I+k]; 
        REAL(Rval_scbd)[j] = sum_col/REAL(SStotal)[0];
        sum_col = 0.0;
    }
    SET_VECTOR_ELT(Rlist,3,Rval_scbd);
    UNPROTECT(7);
    return Rlist;
}
/*  fin calcul LCBD  */

// Calcul de apply(x,2,sample) : permuter les valeurs de chaque colonne ========
SEXP sampleC (SEXP x)
{
    x = PROTECT(coerceVector(x, REALSXP));
    R_len_t i, j, k, I, J, irdm;
    double *rx = REAL(x), *rz, tmp, rdm;
    SEXP z, xperm, Rdim;
    
    Rdim = PROTECT(getAttrib(x, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(xperm = allocMatrix(REALSXP, I, J));
    memset(REAL(xperm),0.0,I*J*sizeof(double));  
    double *Rxperm = REAL(xperm);
    
    PROTECT(z = allocVector(REALSXP,I));
    rz = REAL(z);
    memset(rz,0.0,I*sizeof(double)); 
    
    GetRNGstate();
    
    for (j = 0 ; j <J ; j++)
    {
        for (i = 0; i< I; i++)   rz[i] = rx[j*I+i]; 
        // permuter la colonne
        for (i = 0 ; i < I ; i++)
        {
            do {
                rdm = unif_rand();                 //Sytem RNG
            } while (rdm == 1.0);                  // pour ne pas avoir un indice irdm  = I    
            irdm = (R_len_t) (rdm * I); 
            tmp = rz[irdm];
            rz[irdm]  = rz[i];
            rz[i] = tmp;
        } 
        
        // affecter les colonnes permutées à la matrice de sortie
        for (k=0 ; k<I; k++) {
            Rxperm[j*I+k] = rz[k];
            rz[k] = 0.0;
        }     	
    }
    PutRNGstate(); 
    
    UNPROTECT(4);
    return(xperm);
}

// prépare la liste de sortie pour betadiv1 : 
// 0 : SSTOTAl, 1: BDTOTAL, 2: LCBD, 3: SCBD, 4: p.LCBD, 5: distance euclédienne, 6: méthode ========

SEXP createList1(SEXP list, SEXP vect, SEXP Edist, SEXP method, SEXP perm) {
    SEXP Rlist, names;
    
    vect = PROTECT(coerceVector(vect, INTSXP));
    Edist = PROTECT(coerceVector(Edist, REALSXP));
    method = PROTECT(coerceVector(method,STRSXP));
    perm = PROTECT(coerceVector(perm, INTSXP));
    PROTECT(list);
    int I = length(vect);
    SEXP Pvect = PROTECT(allocVector(REALSXP,I));
    PROTECT(Rlist = allocVector(VECSXP,7));
    SET_VECTOR_ELT(Rlist,0,VECTOR_ELT(list,0));
    SET_VECTOR_ELT(Rlist,1,VECTOR_ELT(list,1));
    SET_VECTOR_ELT(Rlist,2,VECTOR_ELT(list,2));
    SET_VECTOR_ELT(Rlist,3,VECTOR_ELT(list,3));
    for (int k=0;k<I;k++)  REAL(Pvect)[k] = (double)INTEGER(vect)[k]/(asInteger(perm)+1);
    SET_VECTOR_ELT(Rlist,4,Pvect);
    SET_VECTOR_ELT(Rlist,5,Edist);
    SET_VECTOR_ELT(Rlist,6,method);
    PROTECT(names = allocVector(VECSXP,7));
    SET_VECTOR_ELT(names,0,mkChar("SSTOTAL"));
    SET_VECTOR_ELT(names,1,mkChar("BDTOTAL"));
    SET_VECTOR_ELT(names,2,mkChar("LCBD"));
    SET_VECTOR_ELT(names,3,mkChar("SCBD"));
    SET_VECTOR_ELT(names,4,mkChar("p.LCBD"));
    SET_VECTOR_ELT(names,5,mkChar("D"));
    SET_VECTOR_ELT(names,6,mkChar("Method"));
    setAttrib(Rlist,R_NamesSymbol,names);
    
    UNPROTECT(8);
    return(Rlist);
}

/*********************************************************************/
/*********************************************************************/

///// début betadiv1 (Permutation + calcul SSTOTAL, BDTOTAL, LCBD, SCBD et p.LCBD por groupe1) ========
SEXP betadiv1(SEXP x, SEXP method, SEXP perm) {
    R_len_t  I, J;
    PROTECT_INDEX  ixperm, ixtransform, iresultat, iLCBD, ilist, ixsquared, idist;
    SEXP Edist, Rdim, x_perm, x_transform, x_squared, resultat, LCBD, Rlist;
    
    
    x = PROTECT(coerceVector(x, REALSXP));
    
    Rdim = PROTECT(getAttrib(x, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    
    PROTECT_WITH_INDEX(Edist = allocMatrix(REALSXP, I, I),&idist);
    memset(REAL(Edist),0.0,I*I*sizeof(double)); 
    
    method = PROTECT(coerceVector(method, STRSXP)); 
    perm = PROTECT(coerceVector(perm, INTSXP));
    
    PROTECT_WITH_INDEX(x_perm = allocMatrix(REALSXP,I,J), &ixperm);
    memset(REAL(x_perm),0.0,I*J*sizeof(double)); 
    
    PROTECT_WITH_INDEX(x_transform = allocMatrix(REALSXP,I,J), &ixtransform);
    memset(REAL(x_transform),0.0,I*J*sizeof(double)); 
    
    PROTECT_WITH_INDEX(x_squared = allocMatrix(REALSXP,I,J), &ixsquared);
    memset(REAL(x_squared),0.0,I*J*sizeof(double)); 
    
    PROTECT_WITH_INDEX (LCBD= allocVector(REALSXP,I),&iLCBD);
    memset(REAL(LCBD),0.0,I*sizeof(double));
    
    PROTECT_WITH_INDEX (Rlist= allocVector(VECSXP,4),&ilist);
    SET_VECTOR_ELT(Rlist,0,ScalarReal(0.0));
    SET_VECTOR_ELT(Rlist,1,ScalarReal(0.0));
    SET_VECTOR_ELT(Rlist,2,LCBD);
    
    PROTECT_WITH_INDEX (resultat = allocVector(VECSXP,4),&iresultat);
    
    // Transform the original matrix 
    REPROTECT(x_transform = transform_mat(x,method),ixtransform);
    
    // Distance Euclidienne
    REPROTECT(Edist = euclidean(x_transform),idist);
    
    // Squared difference of the transformed matrix
    REPROTECT(x_squared = squared_diff(x_transform),ixsquared);
    
    // Calcul SSTOTAL, BDTOTAL, LCBD et SCBD
    REPROTECT(Rlist = calcul_BD(x_squared),ilist);
    
    // Save the original LCBD
    for (int k=0 ; k<I ; k++)   REAL(LCBD)[k] = REAL(VECTOR_ELT(Rlist,2))[k];
    
    SEXP CPTlcbd = PROTECT(allocVector(INTSXP,I));
    memset(INTEGER(CPTlcbd),0,I*sizeof(int));       
    
    for(int k=0;k<I;k++) INTEGER(CPTlcbd)[k] = 1;
    
    SEXP tmp = PROTECT(allocVector(INTSXP,1));
    memset(INTEGER(tmp),0,sizeof(int)); 
    
    for (R_len_t iperm=0; iperm < asInteger(perm); iperm++)
    {
        REPROTECT(x_perm = sampleC(x),ixperm);
        
        REPROTECT(x_transform = transform_mat(x_perm,method),ixtransform);
        
        // squared difference of the transformed matrix
        REPROTECT(x_squared = squared_diff(x_transform),ixsquared);
        
        REPROTECT(resultat = calcul_BD(x_squared),iresultat);
        
        for (R_len_t k=0; k<I; k++)
        {
            //if (fabs(REAL(VECTOR_ELT(resultat,2))[k]-REAL(LCBD)[k])<=sqrt(DBL_EPSILON) || REAL(VECTOR_ELT(resultat,2))[k]-REAL(LCBD)[k]>sqrt(DBL_EPSILON))
            if (REAL(VECTOR_ELT(resultat,2))[k]+sqrt(DBL_EPSILON)>=REAL(LCBD)[k]) 
            {
                INTEGER(tmp)[0] = INTEGER(CPTlcbd)[k] + 1;
                INTEGER(CPTlcbd)[k] = INTEGER(tmp)[0];
            } 
        } 
    }
    
    
    SEXP ans = PROTECT(createList1(Rlist,CPTlcbd,Edist,method,perm));
    UNPROTECT(14);
    return(ans);
}

/* fin  betadiv1 */

/* ======== ======== ======== ======== ======== ======== ======== ======== ========*/
/* ======== ======== ======== ======== ======== ======== ======== ======== ========*/

/***  debut calcul pour betadiv2  ***/
/* calcul des différentes distances */

///// manhattan distance ========
SEXP manhattan (SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim, D;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), somme = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double)); 
    
    for (i=0; i<I; i++) 
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++) somme = somme + fabs(rx[k*I+i]-rx[k*I+j]);
            REAL(D)[i*I+j] = somme;
            somme = 0.0;
        }
        UNPROTECT(3);
    return D;
}

///// canberra ========
SEXP canberra(SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim, D, p; 
    PROTECT_INDEX ip;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), somme = 0.0, denom = 0.0, nom = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I=INTEGER(Rdim)[0];
    J=INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double)); 
    
    PROTECT_WITH_INDEX(p = allocVector(INTSXP,1),&ip);
    INTEGER(p)[0] = J;
    
    for (i=0; i<I; i++)  
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++)  
            {
                if (rx[k*I+i]==0.0  &&  rx[k*I+j]==0.0)  
                { 
                    INTEGER(p)[0]=INTEGER(p)[0]-1; // Nbre des sites - les doubles zéro
                    continue;
                }
                else
                {
                    nom=fabs(rx[k*I+i]-rx[k*I+j]);
                    denom=rx[k*I+i]+rx[k*I+j];
                    if(denom<DBL_EPSILON) denom=DBL_EPSILON;
                    somme+=nom/denom; 
                }  
            }
            
            REAL(D)[i*I+j]=somme/INTEGER(p)[0];    
            INTEGER(p)[0]=J;
            somme=0.0;
        }
        
        UNPROTECT(4);
    return D;
}

////  kulczynski ========
SEXP kulczynski(SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim,D; 
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), som_min = 0.0, som_coli = 0.0, som_colj = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double)); 
    
    for (i=0; i<I; i++)  
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++) 
            {
                if (rx[k*I+i]<rx[k*I+j]) som_min += rx[k*I+i]; else som_min += rx[k*I+j]; 
                som_coli += rx[k*I+i];
                som_colj += rx[k*I+j];
            }
            if(som_coli<DBL_EPSILON) som_coli=DBL_EPSILON;
            if(som_colj<DBL_EPSILON) som_colj=DBL_EPSILON;
            REAL(D)[i*I+j] = (double)1 - 0.5 * ((som_min/som_coli) + (som_min/som_colj));        
            som_min = som_coli = som_colj = 0.0;
        }
        UNPROTECT(3);
    return D;
}

///// divergence ======== 
SEXP divergence (SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim, D, p; 
    PROTECT_INDEX ip;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), somme = 0.0, denom = 0.0, nom = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double)); 
    
    PROTECT_WITH_INDEX(p = allocVector(INTSXP,1),&ip);
    INTEGER(p)[0] = J;
    
    for (i=0; i<I; i++) 
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++)  
            {
                if (rx[k*I+i] == 0.0  &&  rx[k*I+j] == 0.0)  
                { 
                    INTEGER(p)[0] = INTEGER(p)[0]-1; // Nbre des sites - les doubles zéro
                    continue;
                }
                else
                {
                    nom = rx[k*I+i] - rx[k*I+j];
                    denom = rx[k*I+i] + rx[k*I+j];
                    if(denom<DBL_EPSILON) denom=DBL_EPSILON;
                    somme = somme + pow(nom/denom,2);  
                } 
            } 
            REAL(D)[i*I+j] = sqrt(somme/INTEGER(p)[0]);      
            INTEGER(p)[0] = J;
            somme = 0.0;   
        }
        UNPROTECT(4);
    return D;
}

///// modmean ========
SEXP modmean (SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim, D, p;
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), somme = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double)); 
    
    PROTECT(p = allocVector(INTSXP,1));
    INTEGER(p)[0] = J;
    
    for (i=0; i<I; i++) 
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++) 
            {
                if (rx[k*I+i]  == 0.0  &&  rx[k*I+j] == 0.0)  
                { 
                    INTEGER(p)[0] = INTEGER(p)[0]-1; // Nbre des sites - les doubles zéro
                    continue;
                }
                else  somme += fabs(rx[k*I+i] - rx[k*I+j]);
            }
            REAL(D)[i*I+j] = somme/INTEGER(p)[0];//INTEGER(p)[0];
            INTEGER(p)[0] = J;
            somme = 0.0;  
        }
        UNPROTECT(4);
    return D;
}

//// ruzicka ========
SEXP ruzicka(SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim,D; 
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix);
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double)); 
    
    double A=0.0, sum_all=0.0;
    
    for (i=0; i<I; i++) 
    { 
        for (j=i+1; j<I; j++) 
        {
            for(k=0; k<J; k++)
            {
                if(rx[k*I+i] < rx[k*I+j])  A += rx[k*I+i]; else  A += rx[k*I+j];
                sum_all += (rx[k*I+i]+rx[k*I+j]);
            } 
            double denom = sum_all-A;
            if (denom<DBL_EPSILON) denom=DBL_EPSILON; 
            REAL(D)[i*I+j] = (sum_all-2*A)/denom;
            A = sum_all = 0.0;
        }
    }
    UNPROTECT(3);
    return D;
}

///// whittaker ========
SEXP whittaker(SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim, D; 
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), somme = 0.0, som_coli = 0.0, som_colj = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double)); 
    
    for (i=0; i<I; i++) 
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++) 
            {
                som_coli += rx[k*I+i]; 
                som_colj += rx[k*I+j]; 
            }
            if(som_coli<DBL_EPSILON) som_coli=DBL_EPSILON;
            if(som_colj<DBL_EPSILON) som_colj=DBL_EPSILON;
            for (k = 0; k < J; k++)  somme = somme + fabs((rx[k*I+i]/som_coli - rx[k*I+j]/som_colj));   
            REAL(D)[i*I+j] = (double) 0.5 * somme;       
            
            somme = 0.0;
            som_coli = 0.0;
            som_colj = 0.0;
        }
        UNPROTECT(3);
    return D;
}

///// wishart ========
SEXP wishart(SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim,D; 
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix), prodij = 0.0, prodi = 0.0, prodj = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double)); 
    
    for (i=0; i<I; i++)  
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++) 
            {
                prodij += rx[k*I+i]*rx[k*I+j]; 
                prodi += pow(rx[k*I+i],2); 
                prodj += pow(rx[k*I+j],2);
            } 
            double denom = prodi+prodj-prodij;
            if(denom<DBL_EPSILON) denom=DBL_EPSILON;
            REAL(D)[i*I+j] = 1 - prodij/denom;        
            prodij = prodi = prodj = 0.0;
        }
        
        UNPROTECT(3);
    return D;
}

///// chao ======== 
SEXP chao_C (SEXP RMat, SEXP coeff_C, SEXP sample_C) {
    R_len_t I,J,i,j,k,ind1,ind2;
    
    SEXP Rdim,Rres, Rv1_pa, Rv2_pa;
    
    int *v1_pa, *v2_pa, *shared_sp, *sel2, *sel1;
    int  a1_j,a1_i, a2_j, a2_i, sum_shared_sp;
    
    a1_j=0;
    a1_i=0;
    a2_j=0;
    a2_i=0;
    ind1=0;
    ind2=0;
    sum_shared_sp=0;
    
    RMat = PROTECT(coerceVector(RMat, REALSXP)); 
    double *rx = REAL(RMat);
    
    coeff_C = PROTECT(coerceVector(coeff_C, STRSXP));
    sample_C = PROTECT(coerceVector(sample_C,LGLSXP));
    
    SEXP RC_i = PROTECT(allocVector(REALSXP,1));
    double *C_i = REAL(RC_i);
    C_i[0] = 0.0;
    
    SEXP RC_j = PROTECT(allocVector(REALSXP,1));
    double *C_j = REAL(RC_j);
    C_j[0] = 0.0;
    
    SEXP RS_i = PROTECT(allocVector(REALSXP,1));
    double *S_i = REAL(RS_i);
    S_i[0] = 0.0;
    
    SEXP RS_j = PROTECT(allocVector(REALSXP,1));
    double *S_j = REAL(RS_j);
    S_j[0] = 0.0;
    
    SEXP RN_i = PROTECT(allocVector(REALSXP,1));
    double *N_i = REAL(RN_i);
    N_i[0] = 0.0;
    
    SEXP RN_j = PROTECT(allocVector(REALSXP,1));
    double *N_j = REAL(RN_j);
    N_j[0] = 0.0;
    
    SEXP RU_i = PROTECT(allocVector(REALSXP,1));
    double *U_i = REAL(RU_i);
    U_i[0] = 0.0;
    
    SEXP RU_j = PROTECT(allocVector(REALSXP,1));
    double *U_j = REAL(RU_j);
    U_j[0] = 0.0;
    
    Rdim = PROTECT(getAttrib(RMat, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(Rres = allocMatrix(REALSXP, I, I));
    memset(REAL(Rres),0.0,I*I*sizeof(double)); 
    double *res = REAL(Rres);
    
    PROTECT(Rv1_pa = allocVector(INTSXP,J));
    v1_pa = INTEGER(Rv1_pa);
    
    PROTECT(Rv2_pa = allocVector(INTSXP,J));
    v2_pa = INTEGER(Rv2_pa);
    
    SEXP Rshared_sp = PROTECT(allocVector(INTSXP,J));
    shared_sp = INTEGER(Rshared_sp);
    
    SEXP Rsel1 = PROTECT(allocVector(INTSXP,J));
    sel1 = INTEGER(Rsel1);
    
    SEXP Rsel2 = PROTECT(allocVector(INTSXP,J));
    sel2 = INTEGER(Rsel2);
    
    SEXP RU = PROTECT(allocVector(REALSXP,1));
    double *U = REAL(RU);
    U[0] = 0.0;
    
    SEXP RV = PROTECT(allocVector(REALSXP,1));
    double *V = REAL(RV);
    V[0] = 0.0;
    
    
    for (i=0; i<I; i++) 
    {
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++) 
            {
                INTEGER(Rv1_pa)[k] = 0;
                INTEGER(Rv2_pa)[k] = 0;
                INTEGER(Rsel2)[k] = -1;
                INTEGER(Rsel1)[k] = -1;
                INTEGER(Rshared_sp)[k] = 0;
                
                if (rx[k*I+i] != 0.0)  v1_pa[k] = 1;  //Vector 1 in presence-absence form
                if (rx[k*I+j] != 0.0)  v2_pa[k] = 1;  //Vector 2 in presence-absence form
                
                N_i[0] += rx[k*I+i];                    //Sum of abundanceS_in vector 1
                N_j[0] += rx[k*I+j];                    //Sum of abundanceS_in vector 2
                shared_sp[k] = v1_pa[k] * v2_pa[k];     //Vector of shared_species ("pa")
                sum_shared_sp += shared_sp[k];
            } 
            
            if (sum_shared_sp < DBL_EPSILON)  res[i*I+j] = 1.0;  //sum_shared = 0
            else 
            {
                if(LOGICAL(sample_C)[0])
                {
                    for(k=0;k<J;k++) 
                    {
                        if(shared_sp[k] == 1) 
                        {
                            C_i[0] += rx[k*I+i];                  //Sum of shared_sp. abundanceS_in v1
                            C_j[0] += rx[k*I+j];                  //Sum of shared_sp. abundanceS_in v1
                            if(fabs(rx[k*I+j]-1)<DBL_EPSILON)  a1_i += 1;     // a1_i is Nb SingletonS_in v2 for cmmun species
                            if(fabs(rx[k*I+i]-1)<DBL_EPSILON)  a1_j += 1;     //a1_j is Nb SingletonS_in v1 for cmmun species
                            if(fabs(rx[k*I+j]-2)<DBL_EPSILON)  a2_i += 1;     //a2_i Nb DoubletonS_in V2 for cmmun species
                            if(fabs(rx[k*I+i]-2)<DBL_EPSILON)  a2_j += 1;     //a2_j Nb DoubletonS_in V1 for cmmun species
                        }
                        if(fabs(rx[k*I+j]-1)<DBL_EPSILON) {sel2[ind2]=k; ind2++;}  //sel2 :  singleton position in v2 for commun species
                        if(fabs(rx[k*I+i]-1)<DBL_EPSILON) {sel1[ind1]=k; ind1++;}  //sel1 :  singleton position in v1 for commun species
                    }
                    
                    if (a2_j == 0)  a2_j = 1;       
                    if (a2_i == 0)  a2_i = 1;         
                    
                    if (ind2 > 0)  
                    { 
                        for(R_len_t ind_sel=0;ind_sel<ind2;ind_sel++)   S_i[0] += rx[sel2[ind_sel]*I+i];// Sum abond in v1 for singletons in v2
                    } 
                    else   S_i[0] = 0.0;
                    
                    if (ind1 > 0) 
                    {
                        for(R_len_t ind_sel=0;ind_sel<ind1;ind_sel++)  S_j[0] += rx[sel1[ind_sel]*I+j];  // Sum abond in v2 for singletons in v1  
                    }
                    else   S_j[0] = 0.0; 
                    
                    if (N_j[0]<DBL_EPSILON) N_j[0]=DBL_EPSILON;
                    if (N_i[0]<DBL_EPSILON) N_i[0]=DBL_EPSILON;
                    
                    U_j[0] = (C_j[0]/N_j[0])+((N_i[0]-1)/N_i[0])*((double)a1_j/(2*a2_j))*(S_j[0]/N_j[0]);
                    
                    if ((U_j[0]-1)>DBL_EPSILON)  U_j[0] = 1.0;
                    
                    U_i[0]=(C_i[0]/N_i[0]) +((N_j[0]-1)/N_j[0])*( (double)a1_i/(2*a2_i))*(S_i[0]/N_i[0]);                                          
                    if ((U_i[0]-1.0)>DBL_EPSILON)  U_i[0] = 1.0;  
                    
                    if (U_j[0]<DBL_EPSILON)  U_j[0]=DBL_EPSILON;
                    if (U_i[0]<DBL_EPSILON)  U_i[0]=DBL_EPSILON;
                    
                    if (STRING_ELT(coeff_C,0) == mkChar("Jaccard"))  res[i*I+j] = 1-(U_j[0]*U_i[0])/(U_j[0]+ U_i[0]-(U_j[0]*U_i[0])); 
                    if (STRING_ELT(coeff_C,0) == mkChar("Sorensen")) res[i*I+j] = 1-(2 * U_j[0]*U_i[0]/(U_j[0]+U_i[0]));
                    if (STRING_ELT(coeff_C,0) == mkChar("Ochiai"))      res[i*I+j] = 1-(sqrt(U_j[0]* U_i[0]));  
                    if (STRING_ELT(coeff_C,0) == mkChar("Simpson"))  
                    {
                        double tmp;
                        if ((U_j[0]-U_j[0]*U_i[0]) < (U_i[0]-U_j[0]*U_i[0])) tmp=U_j[0]-(U_j[0]*U_i[0]); else  tmp=U_i[0]-(U_j[0]*U_i[0]);
                        res[i*I+j] = 1 -((U_j[0] * U_i[0])/((U_j[0]*U_i[0]) + tmp));   
                    }               
                } ///fin de sample_C = true
                else
                {
                    ////sample_C = false
                    for(k=0;k<J;k++) 
                    {
                        if(shared_sp[k] == 1) 
                        {
                            U[0] += rx[k*I+i]/N_i[0];                  //Sum of shared_sp. abundanceS_in v1
                            V[0] += rx[k*I+j]/N_j[0];                  //Sum of shared_sp. abundanceS_in v1
                        }  
                    }
                    if (STRING_ELT(coeff_C,0) == mkChar("Jaccard"))  res[i*I+j] = 1-(U[0]*V[0]/(U[0]+V[0]-U[0]*V[0])); 
                    if (STRING_ELT(coeff_C,0) == mkChar("Sorensen")) res[i*I+j] = 1-(2 * U[0]*V[0]/(U[0]+V[0]));
                    if (STRING_ELT(coeff_C,0) == mkChar("Ochiai"))      res[i*I+j] = 1-(sqrt(U[0]*V[0]));  
                    if (STRING_ELT(coeff_C,0) == mkChar("Simpson"))  
                    {
                        double tmp;
                        if ((U[0]-U[0]*V[0]) < (V[0]-U[0]*V[0])) tmp=U[0]-U[0]*V[0]; else  tmp=V[0]-U[0]*V[0];
                        res[i*I+j] = 1 -((U[0]*V[0])/(U[0]*V[0]+tmp));   
                    }               
                } //fin sample_C = false 
            }        
            
            C_i[0]=C_j[0]=S_i[0]=S_j[0]=N_i[0]=N_j[0]=U_i[0]=U_j[0]=0.0;
            U[0]=V[0]=a1_j=a1_i=a2_j=a2_i=sum_shared_sp=ind1=ind2=0;
            
        }
    }
    UNPROTECT(20);
    return (Rres);
} 

///// %difference ========
SEXP percentdiff(SEXP RinMatrix) {
    R_len_t I,J,i,j,k;
    SEXP Rdim,RD; 
    
    RinMatrix = PROTECT(coerceVector(RinMatrix, REALSXP)); 
    double *rx = REAL(RinMatrix);
    
    SEXP Rsomme = PROTECT(allocVector(REALSXP,1));
    double *somme = REAL(Rsomme);
    somme[0] = 0.0;
    
    SEXP Rsom_coli = PROTECT(allocVector(REALSXP,1));
    double *som_coli = REAL(Rsom_coli);
    som_coli[0] = 0.0;
    
    SEXP Rsom_colj = PROTECT(allocVector(REALSXP,1));
    double *som_colj = REAL(Rsom_colj);
    som_colj[0] = 0.0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    PROTECT(RD = allocMatrix(REALSXP, I, I));
    memset(REAL(RD),0.0,I*I*sizeof(double)); 
    
    
    for (i = 0; i < I; i++) 
        for (j = i+1; j < I; j++) 
        {    
            for (k = 0; k < J; k++) 
            {
                som_coli[0] = som_coli[0]+rx[k*I+i]; 
                som_colj[0] = som_colj[0]+rx[k*I+j]; 
                somme[0] += fabs(rx[k*I+i]-rx[k*I+j]);
            }
            
            if (som_coli[0]<DBL_EPSILON) som_coli[0]=DBL_EPSILON;  //pour éviter la division par 0
            if (som_colj[0]<DBL_EPSILON) som_colj[0]=DBL_EPSILON;
            
            REAL(RD)[i*I+j] = somme[0]/(som_coli[0]+som_colj[0]);
            somme[0] = 0.0;
            som_coli[0] = 0.0;
            som_colj[0] = 0.0;  
        }
        UNPROTECT(6);
    return RD;
}

///// Distances binaires jaccard, sorensen et ochiai ========
SEXP binary_D(SEXP RinMatrix, SEXP coef) {
    R_len_t I,J,i,j,k;
    SEXP Rdim, D; 
    
    coef = PROTECT(coerceVector(coef, STRSXP));
    RinMatrix = PROTECT(coerceVector(RinMatrix, INTSXP)); 
    int *rx = INTEGER(RinMatrix), sumA=0, sumB=0, sumC=0;
    
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I=INTEGER(Rdim)[0];
    J=INTEGER(Rdim)[1];
    
    PROTECT(D = allocMatrix(REALSXP, I, I));
    memset(REAL(D),0.0,I*I*sizeof(double));
    
    for (i=0; i<I; i++) 
        for (j=i+1; j<I; j++) 
        {  
            for (k=0; k<J; k++)  
            {
                if (rx[k*I+i]!=0  &&  rx[k*I+j]!=0)  sumA+=1; 
                if (rx[k*I+i]!=0  &&  rx[k*I+j]==0)  sumB+=1; 
                if (rx[k*I+i]==0  &&  rx[k*I+j]!=0)  sumC+=1; 
            }
            if (STRING_ELT(coef,0) == mkChar("jaccard"))       REAL(D)[i*I+j]=sqrt(1-(double)sumA/(sumA+sumB+sumC));  
            else if (STRING_ELT(coef,0) == mkChar("sorensen")) REAL(D)[i*I+j]=sqrt(1-(double) 2*sumA/(2*sumA+sumB+sumC));
            else if (STRING_ELT(coef,0) == mkChar("ochiai"))   REAL(D)[i*I+j]=sqrt(1-sumA/(sqrt((sumA+sumB)*(sumA+sumC)))); 
            sumA=sumB=sumC=0; 
        }
        
        UNPROTECT(4);
    return D;
}

/* fin calcul des différentes distances */

SEXP distCalcul(SEXP RinMatrix, SEXP coeff, SEXP sample_C) {
    SEXP  Rdim,Rval;
    R_len_t I;
    Rdim = PROTECT(getAttrib(RinMatrix, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    
   
    PROTECT(Rval = allocMatrix(REALSXP, I, I));
    memset(REAL(Rval),0.0,I*I*sizeof(double)); 
    
    PROTECT(coeff = AS_CHARACTER(coeff));
    sample_C = PROTECT(coerceVector(sample_C,LGLSXP));
    
    /* PROTECT ..UNPROTECT added on Feb 15 2019 NM */
    if (STRING_ELT(coeff,0) == mkChar("manhattan"))        Rval = manhattan(RinMatrix);
    else if (STRING_ELT(coeff,0) == mkChar("modmean"))     Rval = modmean(RinMatrix);
    else if (STRING_ELT(coeff,0) == mkChar("divergence"))  Rval = divergence(RinMatrix);       
    else if (STRING_ELT(coeff,0) == mkChar("canberra"))    Rval = canberra(RinMatrix);
    else if (STRING_ELT(coeff,0) == mkChar("percentdiff")) Rval = percentdiff(RinMatrix);
    else if (STRING_ELT(coeff,0) == mkChar("kulczynski"))  Rval = kulczynski(RinMatrix);
    else if (STRING_ELT(coeff,0) == mkChar("wishart"))     Rval = wishart(RinMatrix);     
    else if (STRING_ELT(coeff,0) == mkChar("whittaker"))   Rval = whittaker(RinMatrix);
    else if (STRING_ELT(coeff,0) == mkChar("ab.jaccard"))  
    {
        /* added on Feb 25 NM */
        SEXP const_coeff = PROTECT(allocVector(STRSXP, 1));
        SET_STRING_ELT(const_coeff,0,mkChar("Jaccard"));
        Rval = chao_C(RinMatrix,const_coeff,sample_C);
        UNPROTECT(1);
    } 
    else if (STRING_ELT(coeff,0) == mkChar("ab.sorensen")) 
    {
        /* added on Feb 25 NM */
        SEXP const_coeff = PROTECT(allocVector(STRSXP, 1));
        SET_STRING_ELT(const_coeff,0,mkChar("Sorensen"));
        Rval = chao_C(RinMatrix,const_coeff,sample_C);
        UNPROTECT(1);  
    } 
    else if (STRING_ELT(coeff,0) == mkChar("ab.ochiai"))   
    {
        /* added on Feb 25 NM */
        SEXP const_coeff = PROTECT(allocVector(STRSXP, 1));
        SET_STRING_ELT(const_coeff,0,mkChar("Ochiai"));
        Rval = chao_C(RinMatrix,const_coeff,sample_C);
        UNPROTECT(1);      
    } 
    else if (STRING_ELT(coeff,0) == mkChar("ab.simpson"))  
    {
        /* added on Feb 25 NM */
        SEXP const_coeff = PROTECT(allocVector(STRSXP, 1));
        SET_STRING_ELT(const_coeff,0,mkChar("Simpson"));
        Rval = chao_C(RinMatrix,const_coeff,sample_C);
        UNPROTECT(1);
    } 
    else if (STRING_ELT(coeff,0) == mkChar("ruzicka"))  Rval = ruzicka(RinMatrix);
    else if (STRING_ELT(coeff,0) == mkChar("jaccard"))     
    {
        /* added on Feb 25 NM */
        SEXP const_coeff = PROTECT(allocVector(STRSXP, 1));
        SET_STRING_ELT(const_coeff,0,mkChar("jaccard"));
        Rval = binary_D(RinMatrix,const_coeff);
        UNPROTECT(1);  
    }
    else if (STRING_ELT(coeff,0) == mkChar("sorensen"))    
    {
        /* added on Feb 25 NM */
        SEXP const_coeff = PROTECT(allocVector(STRSXP, 1));
        SET_STRING_ELT(const_coeff,0,mkChar("sorensen"));
        Rval = binary_D(RinMatrix,const_coeff);
        UNPROTECT(1);
      
    } 
    else if (STRING_ELT(coeff,0) == mkChar("ochiai"))      
    {
         /* added on Feb 25 NM */
        SEXP const_coeff = PROTECT(allocVector(STRSXP, 1));
        SET_STRING_ELT(const_coeff,0,mkChar("ochiai"));
        Rval = binary_D(RinMatrix,const_coeff);
        UNPROTECT(1);
        
    }  
UNPROTECT(4);
return Rval;   
}

///// produit matriciel pour le centrage de Gower ========
SEXP produit(SEXP X, SEXP Y) {
    
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
    
    F77_CALL(dgemm)(transa, transb, &dimX[0], &dimY[1], &dimX[1], &one, xptr, &dimX[0], yptr, &dimY[0], &zero, resptr, &dimX[0]); 
    //F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,x, &nrx, y, &nry, &zero, z, &nrx);
    
    UNPROTECT(5);
    return(res); 
}

///// centrage de Gower ========
SEXP centre_C(SEXP x) {
    PROTECT_INDEX imatc, itmp;
    R_len_t I,i,j;
    
    SEXP mat_gower, Rdim, mat_centre, mat_delta, tmp;
    
    x = PROTECT(coerceVector(x, REALSXP)); 
    
    Rdim = PROTECT(getAttrib(x, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    
    PROTECT_WITH_INDEX(mat_centre = allocMatrix(REALSXP, I, I), &imatc);
    memset(REAL(mat_centre),0.0,I*I*sizeof(double)); 
    
    PROTECT_WITH_INDEX(tmp = allocMatrix(REALSXP, I, I), &itmp);
    memset(REAL(tmp),0.0,I*I*sizeof(double)); 
    
    PROTECT(mat_delta = allocMatrix(REALSXP, I,I));
    memset(REAL(mat_delta),0.0,I*I*sizeof(double)); 
    
    PROTECT(mat_gower = allocMatrix(REALSXP, I, I));
    memset(REAL(mat_gower),0.0,I*I*sizeof(double));
    
    for(i=0;i<I;i++)
    {
        for(j=0;j<I;j++)
        {
            if(i==j)  REAL(mat_gower)[i*I+j] = (double)1-(double)1/I;  //mat_gower = 1-1/n sur diag sinon = -1/n
            else      REAL(mat_gower)[i*I+j] = (double) -1/I;
            REAL(mat_delta)[i*I+j] = -0.5 * pow(REAL(x)[i*I+j],2);  //mat_delta = -D^2/2
        }
    }
    
    REPROTECT(tmp = produit(mat_gower,mat_delta),itmp);
    REPROTECT(mat_centre = produit(tmp,mat_gower),imatc);            //mat_centre = mat_gower %*% mat_delta %*% mat_gower
    
    UNPROTECT(6);
    return  mat_centre;
}

// prépare la liste de sortie : 
// 0 : SSTOTAl, 1: BDTOTAL, 2: LCBD, 3: p.LCBD, 4: distance, 5: méthode  ========                              

SEXP createList2(SEXP s_total,SEXP btotal, SEXP LCBD_C, SEXP vect, SEXP dist, SEXP coef, SEXP perm) {
    SEXP Rlist, names;
    
    vect = PROTECT(coerceVector(vect, INTSXP));
    dist = PROTECT(coerceVector(dist, REALSXP));
    coef = PROTECT(coerceVector(coef,STRSXP));
    perm = PROTECT(coerceVector(perm, INTSXP));
    
    R_len_t  I = length(vect);
    SEXP Pvect = PROTECT(allocVector(REALSXP,I));
    
    PROTECT(Rlist = allocVector(VECSXP,6));
    SET_VECTOR_ELT(Rlist,0,s_total);
    SET_VECTOR_ELT(Rlist,1,btotal);
    SET_VECTOR_ELT(Rlist,2,LCBD_C);
    
    for (R_len_t k=0;k<I;k++)  REAL(Pvect)[k] = (double) INTEGER(vect)[k] / (asInteger(perm) + 1);
    
    SET_VECTOR_ELT(Rlist,3,Pvect);
    SET_VECTOR_ELT(Rlist,4,dist);
    SET_VECTOR_ELT(Rlist,5,coef);
    
    PROTECT(names = allocVector(VECSXP,6));
    SET_VECTOR_ELT(names,0,mkChar("SSTOTAL"));
    SET_VECTOR_ELT(names,1,mkChar("BDTOTAL"));
    SET_VECTOR_ELT(names,2,mkChar("LCBD"));
    SET_VECTOR_ELT(names,3,mkChar("p.LCBD"));
    SET_VECTOR_ELT(names,4,mkChar("D"));
    SET_VECTOR_ELT(names,5,mkChar("Method"));
    setAttrib(Rlist,R_NamesSymbol,names);
    
    UNPROTECT(7);
    return(Rlist);
}

/*********************************************************************/
/*********************************************************************/

///// début betadiv2 (Permutation + calcul SSTOTAL, BDTOTAL, LCBD, p.LCBD et la distance pour groupe2) ========
SEXP betadiv2(SEXP x, SEXP coef, SEXP perm, SEXP sqrt_D, SEXP sample_C) {
    R_len_t  I, J, k, i, j, iperm;
    PROTECT_INDEX  ixperm, iLCBD, iLCBDperm, iRD, iRDperm, idelta;
    SEXP Rdim, x_perm, LCBD_C, LCBDperm, RD, RDperm, Rdelta, SSTOTAL_C, BDTOTAL_C, SSTOTALperm, BDTOTALperm;
    
    sample_C = PROTECT(coerceVector(sample_C,LGLSXP));
    //if(LOGICAL(sample_C)[0]==NA_LOGICAL) LOGICAL(sample_C)[0]=TRUE;
    
    double stotal = 0.0;
    
    x = PROTECT(coerceVector(x, REALSXP));
    
    Rdim = PROTECT(getAttrib(x, R_DimSymbol)); 
    I = INTEGER(Rdim)[0];
    J = INTEGER(Rdim)[1];
    
    sqrt_D = PROTECT(coerceVector(sqrt_D,LGLSXP));
    
    SEXP tmp_cptlcbd = PROTECT(allocVector(INTSXP,1));
    INTEGER(tmp_cptlcbd)[0] = 0;
    
    SEXP CPTlcbd = PROTECT(allocVector(INTSXP,I));   
    for(k=0;k<I;k++) INTEGER(CPTlcbd)[k] = 1;
    
    coef = PROTECT(coerceVector(coef, STRSXP)); 
    perm = PROTECT(coerceVector(perm, INTSXP));
    
    PROTECT_WITH_INDEX(x_perm = allocMatrix(REALSXP,I,J), &ixperm);
    memset(REAL(x_perm),0.0,I*J*sizeof(double)); 
    
    PROTECT_WITH_INDEX (LCBD_C= allocVector(REALSXP,I),&iLCBD);
    memset(REAL(LCBD_C),0.0,I*sizeof(double));
    
    PROTECT_WITH_INDEX (LCBDperm= allocVector(REALSXP,I),&iLCBDperm);
    memset(REAL(LCBDperm),0.0,I*sizeof(double));
    
    PROTECT_WITH_INDEX(RD = allocMatrix(REALSXP,I,I), &iRD);
    memset(REAL(RD),0.0,I*I*sizeof(double)); 
    
    PROTECT_WITH_INDEX(RDperm = allocMatrix(REALSXP,I,I), &iRDperm);
    memset(REAL(RDperm),0.0,I*I*sizeof(double)); 
    
    PROTECT_WITH_INDEX(Rdelta = allocMatrix(REALSXP,I,I), &idelta);
    memset(REAL(Rdelta),0.0,I*I*sizeof(double)); 
    
    PROTECT (SSTOTAL_C= allocVector(REALSXP,1));
    REAL(SSTOTAL_C)[0] = 0.0;
    
    PROTECT(BDTOTAL_C= allocVector(REALSXP,1));
    REAL(BDTOTAL_C)[0] = 0.0;
    
    PROTECT (SSTOTALperm= allocVector(REALSXP,1));
    REAL(SSTOTALperm)[0] = 0.0;
    
    PROTECT(BDTOTALperm= allocVector(REALSXP,1));
    REAL(BDTOTALperm)[0] = 0.0;
    
    // Calcul de distance selon le coefficient 
    REPROTECT(RD = distCalcul(x,coef,sample_C),iRD);
    
    for (i = 0; i < I; i++)  
    {
        for (j = i+1; j < I; j++)   
        { 
            if (LOGICAL(sqrt_D)[0]==TRUE) 
            {
                double tmp = sqrt(REAL(RD)[i*I+j]);
                REAL(RD)[i*I+j] = tmp; 
            } 
            stotal += pow(REAL(RD)[i*I+j],2); 
            REAL(RD)[j*I+i] = REAL(RD)[i*I+j];  
        }
    }
    
    REAL(SSTOTAL_C)[0] =  stotal/(double)I;
    REAL(BDTOTAL_C)[0] =  REAL(SSTOTAL_C)[0] /((double)I-(double)1);
    REPROTECT(Rdelta = centre_C(RD),idelta);
    if (REAL(SSTOTAL_C)[0]<DBL_EPSILON) REAL(SSTOTAL_C)[0]=DBL_EPSILON;
    for (k=0;k<I;k++) REAL(LCBD_C)[k]=REAL(Rdelta)[k*I+k]/REAL(SSTOTAL_C)[0];
    
    for (iperm=0; iperm < INTEGER(perm)[0]; iperm++)
    {
        stotal = 0.0;
        REPROTECT(x_perm = sampleC(x),ixperm);
        REPROTECT(RDperm = distCalcul(x_perm,coef,sample_C),iRDperm);
        for (i = 0; i<I; i++)  
        {
            for (j = i+1; j<I; j++)   
            { 
                if (asInteger(sqrt_D)==1) 
                {
                    double tmp = sqrt(REAL(RDperm)[i*I+j]);
                    REAL(RDperm)[i*I+j] = tmp;
                } 
                stotal += pow(REAL(RDperm)[i*I+j],2); 
                REAL(RDperm)[j*I+i] = REAL(RDperm)[i*I+j];  
            }
        }
        REAL(SSTOTALperm)[0] = stotal/(double)I;
        REAL(BDTOTALperm)[0] = REAL(SSTOTALperm)[0]/((double)I-1);
        REPROTECT(Rdelta = centre_C(RDperm),idelta);
        
        for (k=0; k<I; k++)
        {
            if (REAL(SSTOTALperm)[0]<DBL_EPSILON) REAL(SSTOTALperm)[0]=DBL_EPSILON;
            REAL(LCBDperm)[k] = REAL(Rdelta)[k*I+k]/REAL(SSTOTALperm)[0];
            // if (fabs(REAL(LCBDperm)[k]-REAL(LCBD_C)[k])<=sqrt(DBL_EPSILON) || REAL(LCBDperm)[k]-REAL(LCBD_C)[k]>sqrt(DBL_EPSILON))
            if (REAL(LCBDperm)[k]+sqrt(DBL_EPSILON) >= REAL(LCBD_C)[k]) 
            {
                INTEGER(tmp_cptlcbd)[0]=INTEGER(CPTlcbd)[k]+1;
                INTEGER(CPTlcbd)[k]=INTEGER(tmp_cptlcbd)[0];
            } 
        }
        
    }
    
    
    SEXP ans =  PROTECT(createList2(SSTOTAL_C,BDTOTAL_C,LCBD_C,CPTlcbd,RD,coef,perm));
    UNPROTECT(19);
    return(ans);
}  /* fin betadiv2 */

