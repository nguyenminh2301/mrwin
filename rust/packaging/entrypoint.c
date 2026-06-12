// Template entrypoint for the extendr-linked mrwin package. Copy into
// `src/entrypoint.c` when wiring the Rust backend (see ../README.md).
//
// extendr's `extendr_module!` macro emits `R_init_mrwinrust_extendr`, which
// registers the `wrap__*` call routines. We forward the package init to it so
// the static library is not dropped by the linker. If the package later gains
// other native routines, register them here as well.

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_mrwinrust_extendr(DllInfo *dll);

void R_init_mrwin(DllInfo *dll) {
    R_init_mrwinrust_extendr(dll);
    R_useDynamicSymbols(dll, FALSE);
}
