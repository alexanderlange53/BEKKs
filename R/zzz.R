.onUnload <- function(libpath) {
  library.dynam.unload("BEKKs", libpath)
  invisible(NULL)
}
