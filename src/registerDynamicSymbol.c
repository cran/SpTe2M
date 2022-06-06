// RegisteringDynamic Symbols

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(cvmspe)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(modcv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(spteewks)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sptellks)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sptewls)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sptewme)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"cvmspe",   (DL_FUNC) &F77_NAME(cvmspe),   14},
    {"modcv",    (DL_FUNC) &F77_NAME(modcv),    13},
    {"spteewks", (DL_FUNC) &F77_NAME(spteewks), 11},
    {"sptellks", (DL_FUNC) &F77_NAME(sptellks), 13},
    {"sptewls",  (DL_FUNC) &F77_NAME(sptewls),  17},
    {"sptewme",  (DL_FUNC) &F77_NAME(sptewme),  14},
    {NULL, NULL, 0}
};

void R_init_SpTe2M(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
