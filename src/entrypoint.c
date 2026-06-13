// Package entry point for mrwin.
//
// When the package is built with the Rust backend enabled, ./configure adds
// -DMRWIN_HAVE_RUST and links the static library produced by the `mrwinrust`
// extendr crate. extendr's `extendr_module!` macro emits
// `R_init_mrwinrust_extendr`, which registers the `wrap__*` call routines; we
// forward the package init to it so the static library is retained by the
// linker and its routines become callable from R via .Call().
//
// When cargo is unavailable, ./configure omits -DMRWIN_HAVE_RUST and this file
// compiles to a no-op initializer, so the package still installs with the
// dense/sparse R backends only and backend = "rust" reports an actionable error.

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifdef MRWIN_HAVE_RUST
void R_init_mrwinrust_extendr(DllInfo *dll);
#endif

void R_init_mrwin(DllInfo *dll) {
#ifdef MRWIN_HAVE_RUST
    R_init_mrwinrust_extendr(dll);
#endif
    R_useDynamicSymbols(dll, TRUE);
    R_forceSymbols(dll, FALSE);
}
