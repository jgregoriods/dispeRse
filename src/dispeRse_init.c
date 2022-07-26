#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern void run_model(int* nrow, int* ncol, double* population, double* env, int* arrival, double* r, double* phi, int* start, int* x, int* y, int* iter, int* num_origins, double* t, int* terrain);

static const R_CMethodDef CEntries[] = {
    {"run_model", (DL_FUNC) &run_model, 14},
    {NULL, NULL, 0}
};

void R_init_dispeRse(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
