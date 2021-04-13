stopifNegative <- function(arg, allowNA=FALSE, allowZero=TRUE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(allowZero) {
      if(any(arg < 0, na.rm=TRUE)) {
        if(allowNA) {
          stop("Argument '", name, "' must be non-negative, or NA.", call.=FALSE)
        } else {
          stop("Argument '", name, "' must be non-negative.", call.=FALSE)
        }
      }
    } else {
      if(any(arg <= 0, na.rm=TRUE)) {
        if(allowNA) {
          stop("Argument '", name, "' must be greater than 0, or NA.", call.=FALSE)
        } else {
          stop("Argument '", name, "' must be greater than 0.", call.=FALSE)
        }
      }
    }
  }
}
