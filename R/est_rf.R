#' A recombination fraction function
#'
#' This function computes the recombination fraction matrix for backcroos, F2 and outcrossing expertimental populations.
#' @param dat an object of class \code{outcross}, \code{bc.onemap},
##' \code{f2.onemap}, \code{riself.onemap} or \code{risib.onemap}.
#' @keywords recombination
#' @export
#' @examples
#' dat<-simulate_pop(type="f2", n.ind = 250,
#'                   n.mrk = 1000, ch.len = 200,
#'                   dom43 = 20, dom51 = 20,
#'                   missing = 15, n.ch = 2,
#'                   verbose = TRUE)
#' require(onemap)
#' dat
#' r<-est_rf(dat)


est_rf<-function(dat)
{
  ## checking for correct object
  if(!any(class(dat)=="outcross"||class(dat)=="f2.onemap" || class(dat)=="bc.onemap")) stop(deparse(substitute(dat))," is not an object of class 'outcross', 'bc.onemap', 'f2.onemap', 'riself.onemap' or 'risib.onemap'")
  if (dat$n.mar<2) stop("there must be at least two markers to proceed with analysis")

  genoR<-dat$geno
  if(class(dat)=="bc.onemap")
  {
    r<-est_bc(genoR, dat$n.ind)
    return(structure(r, class="bc.rf"))
  }
  else if(class(dat)=="f2.onemap")
  {
    r<-est_f2(genoR, dat$segr.type.num, dat$n.ind)
    return(structure(r, class="f2.rf"))
  }
  else if(class(dat)=="outcross")
  r<-est_out(genoR, dat$segr.type.num, dat$n.ind)
  return(structure(r, class="out.rf"))
}
