stopifnotProbability <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(any(arg < 0 | arg > 1, na.rm=TRUE)) {
      if(allowNA) {
        stop("Argument '", name, "' must be a probability between 0 and 1, or NA.", call.=FALSE)
      } else {
        stop("Argument '", name, "' must be a probability between 0 and 1.", call.=FALSE)
      }
    }
  }
}
