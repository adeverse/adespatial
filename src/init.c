#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void buildbinary(void *, void *, void *, void *, void *, void *);
extern void forwardsel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void testglobal(void *, void *, void *, void *, void *, void *, void *);
extern void C_WRperiodogram(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cclust(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP betadiv1(SEXP, SEXP, SEXP);
extern SEXP betadiv2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP binary_D(SEXP, SEXP);
extern SEXP canberra(SEXP);
extern SEXP chao_C(SEXP, SEXP, SEXP);
extern SEXP divergence(SEXP);
extern SEXP euclidean(SEXP);
extern SEXP kulczynski(SEXP);
extern SEXP manhattan(SEXP);
extern SEXP modmean(SEXP);
extern SEXP percentdiff(SEXP);
extern SEXP ruzicka(SEXP);
extern SEXP s_loop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sti_loop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP t_loop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP transform_mat(SEXP, SEXP);
extern SEXP whittaker(SEXP);
extern SEXP wishart(SEXP);

/* FORTRAN calls */

void F77_NAME(geodistv)(double*, double*, double*, double*, double*, int*);

static const R_CMethodDef CEntries[] = {
    {"buildbinary",   (DL_FUNC) &buildbinary,    6},
    {"forwardsel",    (DL_FUNC) &forwardsel,    18},
    {"testglobal",    (DL_FUNC) &testglobal,     7},
    {"C_WRperiodogram", (DL_FUNC) &C_WRperiodogram,  9},
    {"cclust",          (DL_FUNC) &cclust,          11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"betadiv1",      (DL_FUNC) &betadiv1,       3},
    {"betadiv2",      (DL_FUNC) &betadiv2,       5},
    {"binary_D",      (DL_FUNC) &binary_D,       2},
    {"canberra",      (DL_FUNC) &canberra,       1},
    {"chao_C",        (DL_FUNC) &chao_C,         3},
    {"divergence",    (DL_FUNC) &divergence,     1},
    {"euclidean",     (DL_FUNC) &euclidean,      1},
    {"kulczynski",    (DL_FUNC) &kulczynski,     1},
    {"manhattan",     (DL_FUNC) &manhattan,      1},
    {"modmean",       (DL_FUNC) &modmean,        1},
    {"percentdiff",   (DL_FUNC) &percentdiff,    1},
    {"ruzicka",       (DL_FUNC) &ruzicka,        1},
    {"s_loop",        (DL_FUNC) &s_loop,        14},
    {"sti_loop",      (DL_FUNC) &sti_loop,      11},
    {"t_loop",        (DL_FUNC) &t_loop,        14},
    {"transform_mat", (DL_FUNC) &transform_mat,  2},
    {"whittaker",     (DL_FUNC) &whittaker,      1},
    {"wishart",       (DL_FUNC) &wishart,        1},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
    {"geodistv",      (DL_FUNC) &F77_NAME(geodistv), 6},
    {NULL, NULL, 0}
};

void R_init_adespatial(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
