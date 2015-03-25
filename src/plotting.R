load.local <- function(file) {
 v <- load(file)
 stopifnot(length(v) == 1)
 get(v)
}

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))
