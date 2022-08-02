/* Registers the function run_model exported to R */

#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern void run_model(int*, int*, double*, int*, double*, int*, int*, int*, int*, int*, int*, double*, double*, double*, int*, double*, int*);

static const R_CMethodDef CEntries[] = {
    {"run_model", (DL_FUNC) &run_model, 17},
    {NULL, NULL, 0}
};

void R_init_dispeRse(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
