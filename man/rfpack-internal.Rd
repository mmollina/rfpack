\name{rfpack-internal}

\alias{check.type}


\alias{est_bc}
\alias{est_f2}
\alias{est_out}
\alias{imf.k}
\alias{mf.k}
\alias{sim.gam}
\alias{sim.pop.bc}
\alias{sim.pop.f2}
\alias{sim.pop.out}

\title{Internal rfpack functions}
\description{
  Functions from \pkg{rfpack} not to be directly called by the user.
}
\usage{
est_out(geno, seg_type=NULL, nind)
est_f2(geno, seg_type=NULL, nind)
est_bc(geno, nind)
mf.k (d)
imf.k (r)
sim.gam(n.mrk, r, init=0)
sim.pop.bc(n.ind, n.mrk, ch.len, missing=0, n.ch=1, verbose=TRUE)
check.type(x)
sim.pop.f2(n.ind, n.mrk, ch.len, dom43=0, dom51=0, missing=0, n.ch=1, verbose=TRUE)
sim.pop.out(n.ind, n.mrk, ch.len, missing=0, prob=c(1,1,1,1,1,1,1), n.ch=1, verbose=TRUE)
}
\author{Marcelo Mollinari}
\keyword{internal}
