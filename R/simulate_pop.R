#' A Simulation function
#'
#' This function simulates experimental crosses: backcroos, F2 and outcrossing populations.
#' @param type type of cross. Must be one of the following: "bc" fro backcross, "f2" for f2 and "out" fro outcross.
#' @param n.ind number of individuals
#' @param n.mrk number of markers
#' @param ch.len length of chromosome in centimorgans
#' @param missing percentage of missing data
#' @param n.ch number of chromosomes
#' @param dom43 if \code{type=="f2"}, percentage of dominant markers NOT BB, BB
#' @param dom51 if \code{type=="f2"}, percentage of dominant markers NOT AA, AA
#' @param prob if \code{type=="out"}, probability of markers types A, B1, B2, B3, C, D1 and D2
#' @param verbose logical. If \code{TRUE}, prints the progress of the simulation per chromosome.
#' @export
#' @examples
#' dat<-simulate_pop(type="f2", n.ind = 250,
#'                   n.mrk = 1000, ch.len = 200,
#'                   missing = 15, dom43 = 20,
#'                   dom51 = 20, n.ch = 2,
#'                   verbose = TRUE)
#' require(onemap)
#' dat


simulate_pop<-function(type=c("bc", "f2", "out"), n.ind, n.mrk, ch.len, missing=0, n.ch=1, dom43=0, dom51=0, prob=c(1,1,1,1,1,1,1), verbose=TRUE)
{
    type <- match.arg(type)
    if(type=="bc")
        return(sim.pop.bc(n.ind=n.ind, n.mrk=n.mrk, ch.len=ch.len, missing=missing, n.ch=n.ch, verbose=verbose))
    else if(type=="f2")
        return(sim.pop.f2(n.ind=n.ind, n.mrk=n.mrk, ch.len=ch.len, dom43=dom43, dom51=dom51, missing=missing, n.ch=n.ch, verbose=verbose))
    else if(type=="out")
        return(sim.pop.out(n.ind=n.ind, n.mrk=n.mrk, ch.len=ch.len, missing=missing, prob=prob, n.ch=n.ch, verbose=verbose))
}
