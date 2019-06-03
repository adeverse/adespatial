#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>

void taballoc (double ***tab, int l1, int c1);
void vecintalloc (int **vec, int n);
void prodmatABC (double **a, double **b, double **c);
void freeintvec (int *vec);
void freetab (double **tab);
void aleapermutmat (double **a);
void prodmatAtAB (double **a, double **b);
void prodmatAAtB (double **a, double **b);
double alea (void);
void prodatBc(double *veca, double **matB,double *vecc); 
void vecalloc (double **vec, int n);
void freevec (double *vec);
void aleapermutvec (double *a);
void tabintalloc (int ***tab, int l1, int c1);
void freeinttab (int **tab);
