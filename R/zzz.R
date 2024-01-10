
.onLoad <- function(libname, pkgname) {

  # CRAN OMP THREAD LIMIT
  Sys.setenv("OMP_THREAD_LIMIT" = 1)
}
.onUnload <- function(libpath) {
  library.dynam.unload("BEKKs", libpath)
  invisible(NULL)
}
