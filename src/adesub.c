#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "adesub.h"


//* =========== fonctions adesub ============= */

/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
 * Allocation de memoire dynamique pour un tableau (l1, c1)
 --------------------------------------------------*/
{
    int i, j;
    
    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
        for (i=0;i<=l1;i++) {
            if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
                return;
                for (j=0;j<i;j++) {
                    free(*(*tab+j));
                }
            }
        }
    }
    
    **(*tab) = l1;
    **(*tab+1) = c1;
}

/***********************************************************************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
 * Allocation de memoire pour un vecteur d'entiers de longueur n
 --------------------------------------------------*/
{
    if ( (*vec = (int *) calloc(n+1, sizeof(int))) != NULL) {
        **vec = n;
        return;
    } else {
        return;
    }
}

/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
 * Produit matriciel AB
 --------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (i=1;i<=lig;i++) {
        for (k=1;k<=col2;k++) {
            s = 0;
            for (j=1;j<=col;j++) {
                s = s + a[i][j] * b[j][k];
            }
            c[i][k] = s;
        }       
    }
}
/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
 * liberation de memoire pour un vecteur
 --------------------------------------------------*/
{
    
    free((char *) vec);
    
}
/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
 * Allocation de memoire dynamique pour un tableau (l1, c1)
 --------------------------------------------------*/
{
    int     i, n;
    
    n = *(*(tab));
    for (i=0;i<=n;i++) {
        free((char *) *(tab+i) );
    }
    free((char *) tab);
}
/*************************/
void aleapermutmat (double **a)
{
    /* permute au hasard les lignes du tableau a
     Manly p. 42 le tableau est modifi? */
    int lig, i,j, col, n, k;
    double z;
    
    lig = a[0][0];
    col = a[1][0];
    for (i=1; i<=lig-1; i++) {
        j=lig-i+1;
        k = (int) (j*alea ()+1);
        /*k = (int) (j*genrand()+1);*/
        if (k>j) k=j;
        for (n=1; n<=col; n++) {
            z = a[j][n];
            a[j][n]=a[k][n];
            a[k][n] = z;
        }
    }
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
 * Produit matriciel AtA
 --------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=col;j++) {
        for (k=j;k<=col;k++) {
            s = 0;
            for (i=1;i<=lig;i++) {
                s = s + a[i][k] * a[i][j];
            }
            b[j][k] = s;
            b[k][j] = s;
        }       
    }
}
/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
 * Produit matriciel B = AAt
 --------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=lig;j++) {
        for (k=j;k<=lig;k++) {
            s = 0;
            for (i=1;i<=col;i++) {
                s = s + a[j][i] * a[k][i];
            }
            b[j][k] = s;
            b[k][j] = s;
        }       
    }
}
/**************************/
double alea (void)
{
    double w;
    GetRNGstate();
    w = unif_rand();
    PutRNGstate();
    return (w);
}


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

void tabintalloc (int ***tab, int l1, int c1)
    /*--------------------------------------------------
     * Allocation de memoire dynamique pour un tableau
     * d'entiers (l1, c1)
     --------------------------------------------------*/
{
    int     i, j;
    
    *tab = (int **) calloc(l1+1, sizeof(int *));
    
    if ( *tab != NULL) {
        for (i=0;i<=l1;i++) {
            
            *(*tab+i)=(int *) calloc(c1+1, sizeof(int));        
            if ( *(*tab+i) == NULL ) {
                for (j=0;j<i;j++) {
                    free(*(*tab+j));
                }
                return;
            }
        }
    } else return;
    **(*tab) = l1;
    **(*tab+1) = c1;
    for (i=1;i<=l1;i++) {
        for (j=1;j<=c1;j++) {
            (*tab)[i][j] = 0;
        }
    }
}

/***********************************************************************/
void freeinttab (int **tab)
    /*--------------------------------------------------
     * Allocation de memoire dynamique pour un tableau
     --------------------------------------------------*/
{
    int     i, n;
    
    n = *(*(tab));
    
    for (i=0;i<=n;i++) {
        free((char *) *(tab+i) );
    }
    
    free((char *) tab);
}


