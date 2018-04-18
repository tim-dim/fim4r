#-----------------------------------------------------------------------
# File    : fim4r.R
# Contents: COntinuous time ClOsed Neuron Assembly Detection for R
# Author  : Christian Borgelt
# History : 2015.08.23 file created
#-----------------------------------------------------------------------

# item appearance indicators
fim4r.apps      <- c("-", "ign", "ignore",
                     "n", "none", "neither",
                     "i", "in", "inp", "input",
                     "o",       "out", "output",
                     "a", "ante", "antecedent",
                     "c", "cons", "consequent",
                     "b", "body",
                     "h", "head",
                     "x", "io", "i&o", "o&i", "inout", "in&out",
                     "ac", "a&c", "c&a", "canda",
                     "bh", "b&h", "h&b", "both")

# target type identifiers
fim4r.targets   <- c("s", "set", "sets",
                     "a", "all", "allset", "allsets",
                     "f", "frq", "freq", "frequent",
                     "frqset", "frqsets", "freqset", "freqsets",
                     "c", "cls", "clsd", "closed",
                     "m", "max", "maxi", "maximal",
                     "g", "gen", "gens", "generators",
                     "r", "rule", "rules", "arule", "arules")

# statistics (for accretion)
fim4r.stats     <- c("x", "none",
                     "n", "c", "X2", "chi2",
                     "p", "X2pval", "chi2pval",
                     "y", "yates",
                     "t", "yatespval",
                     "i", "info",
                     "g", "gpval", "infopval",
                     "f", "fetprob",
                     "h", "fetchi2", "fetX2",
                     "m", "fetinfo",
                     "s", "fetsupp")

# evaluation measures
fim4r.evals     <- c("x", "none",
                     "b", "ldratio" )

# evaluation measures
fim4r.xevals    <- c("x", "none",
                     "o", "supp", "support",
                     "c", "conf", "confidence",
                     "d", "confdiff",
                     "l", "lift",
                     "a", "liftdiff",
                     "q", "liftquot",
                     "v", "conviction", "cvct",
                     "e", "cvctdiff",
                     "r", "cvctquot",
                     "k", "cprob", "condprob",
                     "j", "importance", "import",
                     "z", "certainty", "cert",
                     "n", "X2", "chi2",
                     "p", "X2pval", "chi2pval",
                     "y", "yates",
                     "t", "yatespval",
                     "i", "info",
                     "g", "gpval", "infopval",
                     "f", "fetprob",
                     "h", "fetchi2", "fetX2",
                     "m", "fetinfo",
                     "s", "fetsupp",
                     "b", "ldratio" )

# aggregation modes
fim4r.aggs      <- c("x", "none",
                     "m", "min", "minimum",
                     "n", "max", "maximum",
                     "a", "avg", "average")

# FIM algorithm variants
fim4r.algo.apr  <- c("a", "auto",
                     "b", "basic")
fim4r.algo.ecl  <- c("a", "auto",
                     "e", "lists",
                     "i", "tids",
                     "b", "bits",
                     "t", "table",
                     "s", "simple",
                     "r", "ranges",
                     "o", "occdlv", "lcm",
                     "d", "diffs", "diffsets")
fim4r.algo.fpg  <- c("a", "auto",
                     "s", "simple",
                     "c", "complex",
                     "d", "single",
                     "t", "topdown")
fim4r.algo.sam  <- c("a", "auto",
                     "s", "basic",
                     "b", "bsearch",
                     "d", "double",
                     "t", "tree")
fim4r.algo.rem  <- c("a", "auto",
                     "b", "basic")
fim4r.algo.carp <- c("a", "auto",
                     "t", "table",
                     "l", "lists", "tidlist", "tidlists")
fim4r.algo.ista <- c("a", "auto", "prefix", "patricia")

# surrogate data generation methods
fim4r.surrs     <- c("i", "r", "s", "h",
                     "i", "ident",  "identity",
                     "r", "random", "randomize",
                     "s", "swap",
                     "p", "perm", "permute",
                     "h", "shuffle")

# pattern set reduction function identifiers
fim4r.reds      <- c("x", "none",
                     "c", "coins",  "coins0",
                     "C",           "coins1", "coins+1",
                     "i", "items",  "items2",
                     "s", "cover",  "cover0", "covered", "covered0",
                     "S",           "cover1",            "covered1",
                     "l", "leni",   "leni0",  "lenient", "lenient0",
                     "L",           "leni1",             "lenient1",
                     "t", "strict", "strict0",
                     "T",           "strict1")

#-----------------------------------------------------------------------

fim4r.fim <- function (tracts, wgts=NULL, target="s",
                       supp=10.0, zmin=0, zmax=-1, report="a",
                       eval="x", agg="x", thresh=10.0, border=NULL)
{                               # --- wrapper for generic FIM algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp)     &&  (supp <= 100))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.xevals))
  #stopifnot(is.character(agg)    && any(agg[1]  == fim4r.aggs))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.null(border)      || is.numeric(border))
  # call the C implementation:
  r = .Call("f4r_fim", tracts, wgts, target, supp, zmin, zmax,
                       report, eval, agg, thresh, border)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.fim()

#-----------------------------------------------------------------------

fim4r.arules <- function (tracts, wgts=NULL,
                          supp=10.0, conf=80.0, zmin=1, zmax=-1,
                          report="aC", eval="x", thresh=10.0,
                          mode="", appear=NULL)
{                               # --- wrapper for generic FIM algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.numeric(supp)     &&  (supp <= 100))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.xevals))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.character(mode))
  #stopifnot(is.null(appear)
  #  || (is.list(appear) && (length(appear) >= 2)
  #  && (is.character(appear[[2]]) || is.integer(appear[[2]]))
  #  && (typeof(appear[[1]]) == typeof(tracts[[1]]))
  #  && (length(appear[[1]]) == length(appear[[2]]))))
  # call the C implementation:
  r = .Call("f4r_arules", tracts, wgts, supp, conf, zmin, zmax,
                          report, eval, thresh, mode, appear)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.arules()

#-----------------------------------------------------------------------

fim4r.apriori <- function (tracts, wgts=NULL, target="s", supp=10.0,
                           conf=80.0, zmin=0, zmax=-1, report="a",
                           eval="x", agg="x", thresh=10.0, prune=NA,
                           algo="a", mode="", border=NULL, appear=NULL)
{                               # --- wrapper for apriori algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp)     &&  (supp <= 100))
  #stopifnot(is.numeric(conf)     && (conf >= 0)  && (conf <= 100))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.xevals))
  #stopifnot(is.character(agg)    && any(agg[1]  == fim4r.aggs))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.null(prune)       || is.integer(prune))
  #stopifnot(is.character(algo)   && any(algo[1] == fim4r.algo.apr))
  #stopifnot(is.null(border)      || is.numeric(border))
  #stopifnot(is.null(appear)
  #  || (is.list(appear) && (length(appear) >= 2)
  #  && (is.character(appear[[2]]) || is.integer(appear[[2]]))
  #  && (typeof(appear[[1]]) == typeof(tracts[[1]]))
  #  && (length(appear[[1]]) == length(appear[[2]]))))
  # call the C implementation:
  r = .Call("f4r_apriori", tracts, wgts, target, supp, conf,
                           zmin, zmax, report, eval, agg, thresh,
                           prune, algo, mode, border, appear)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.apriori()

#-----------------------------------------------------------------------

fim4r.eclat <- function (tracts, wgts=NULL, target="s", supp=10.0,
                         conf=80.0, zmin=0, zmax=-1, report="a",
                         eval="x", agg="x", thresh=10.0, prune=NA,
                         algo="a", mode="", border=NULL, appear=NULL)
{                               # --- wrapper for eclat algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp))
  #stopifnot(is.numeric(conf)     && (conf >= 0)  && (conf <= 100))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.xevals))
  #stopifnot(is.character(agg)    && any(agg[1]  == fim4r.aggs))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.null(prune)       || is.integer(prune))
  #stopifnot(is.character(algo)   && any(algo[1] == fim4r.algo.ecl))
  #stopifnot(is.null(border)      || is.numeric(border))
  #stopifnot(is.null(appear)
  #  || (is.list(appear) && (length(appear) >= 2)
  #  && (is.character(appear[[2]]) || is.integer(appear[[2]]))
  #  && (typeof(appear[[1]]) == typeof(tracts[[1]]))
  #  && (length(appear[[1]]) == length(appear[[2]]))))
  # call the C implementation:
  r = .Call("f4r_eclat", tracts, wgts, target, supp, conf,
                         zmin, zmax, report, eval, agg, thresh,
                         prune, algo, mode, border, appear)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.eclat()

#-----------------------------------------------------------------------

fim4r.fpgrowth <- function (tracts, wgts=NULL, target="s", supp=10.0,
                            conf=80.0, zmin=0, zmax=-1, report="a",
                            eval="x", agg="x", thresh=10.0, prune=NA,
                            algo="a", mode="", border=NULL, appear=NULL)
{                               # --- wrapper for eclat algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp))
  #stopifnot(is.numeric(conf)     && (conf >= 0)  && (conf <= 100))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.xevals))
  #stopifnot(is.character(agg)    && any(agg[1]  == fim4r.aggs))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.null(prune)       || is.integer(prune))
  #stopifnot(is.character(algo)   && any(algo[1] == fim4r.algo.fpg))
  #stopifnot(is.null(border)      || is.numeric(border))
  #stopifnot(is.null(appear)
  #  || (is.list(appear) && (length(appear) >= 2)
  #  && (is.character(appear[[2]]) || is.integer(appear[[2]]))
  #  && (typeof(appear[[1]]) == typeof(tracts[[1]]))
  #  && (length(appear[[1]]) == length(appear[[2]]))))
  # call the C implementation:
  r = .Call("f4r_fpgrowth", tracts, wgts, target, supp, conf,
                            zmin, zmax, report, eval, agg, thresh,
                            prune, algo, mode, border, appear)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.fpgrowth()

#-----------------------------------------------------------------------

fim4r.sam <- function (tracts, wgts=NULL, target="s", supp=10.0,
                       zmin=0, zmax=-1, report="a",
                       eval="x", thresh=10.0, algo="a", mode="",
                       border=NULL)
{                               # --- wrapper for eclat algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.evals))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.character(algo)   && any(algo[1] == fim4r.algo.sam))
  #stopifnot(is.null(border)      || is.numeric(border))
  # call the C implementation:
  r = .Call("f4r_sam", tracts, wgts, target, supp, zmin, zmax,
                       report, eval, thresh, algo, mode, border)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.sam()

#-----------------------------------------------------------------------

fim4r.relim <- function (tracts, wgts=NULL, target="s", supp=10.0,
                         zmin=0, zmax=-1, report="a",
                         eval="x", thresh=10.0, algo="a", mode="",
                         border=NULL)
{                               # --- wrapper for eclat algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.evals))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.character(algo)   && any(algo[1] == fim4r.algo.rem))
  #stopifnot(is.null(border)      || is.numeric(border))
  # call the C implementation:
  r = .Call("f4r_relim", tracts, wgts, target, supp, zmin, zmax,
                         report, eval, thresh, algo, mode, border)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.relim()

#-----------------------------------------------------------------------

fim4r.carpenter <- function (tracts, wgts=NULL, target="c", supp=10.0,
                             zmin=0, zmax=-1, report="a",
                             eval="x", thresh=10.0, algo="a", mode="",
                             border=NULL)
{                               # --- wrapper for eclat algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.evals))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.character(algo) && any(algo[1] == fim4r.algo.carp))
  #stopifnot(is.null(border)      || is.numeric(border))
  # call the C implementation:
  r = .Call("f4r_carpenter", tracts, wgts, target, supp, zmin, zmax,
                             report, eval, thresh, algo, mode, border)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.carpenter()

#-----------------------------------------------------------------------

fim4r.ista <- function (tracts, wgts=NULL, target="c", supp=10.0,
                        zmin=0, zmax=-1, report="a",
                        eval="x", thresh=10.0, algo="a", mode="",
                        border=NULL)
{                               # --- wrapper for eclat algorithm
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.character(eval)   && any(eval[1] == fim4r.evals))
  #stopifnot(is.numeric(thresh))
  #stopifnot(is.character(algo)   && any(algo[1] == fim4r.algo.ista))
  #stopifnot(is.null(border)      || is.numeric(border))
  # call the C implementation:
  r = .Call("f4r_ista", tracts, wgts, target, supp, zmin, zmax,
                        report, eval, thresh, algo, mode, border)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.ista()

#-----------------------------------------------------------------------

fim4r.genpsp <- function (tracts, wgts=NULL, target="s", supp=10.0,
                          zmin=1, zmax=-1, report="|",
                          cnt=1000, surr="s", seed=0, cpus=0)
{                               # --- wrapper for pattern spectrum
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target)  && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.numeric(cnt))
  #stopifnot(is.character(surr)   && any(surr[1]   == fim4r.surrs))
  #stopifnot(is.numeric(seed)     && is.numeric(cpus))
  # call the C implementation:
  r = .Call("f4r_genpsp", tracts, wgts, target, supp, zmin, zmax,
                          report, cnt, surr, seed, cpus)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.genpsp()

#-----------------------------------------------------------------------

fim4r.estpsp <- function (tracts, wgts=NULL, target="s", supp=10.0,
                          zmin=1, zmax=-1, report="|", equiv=10000,
                          alpha=0.5, smpls=1000, seed=0)
{                               # --- wrapper for pat. spec. estimation
  # check the function arguments:
  #stopifnot(is.list(tracts)      && (length(tracts) > 0))
  #stopifnot(all(sapply(tracts,is.integer))
  #  ||      all(sapply(tracts,is.character))
  #stopifnot(is.null(wgts)
  #  || (is.integer(wgts) && (length(wgts) == length(tracts))))
  #stopifnot(is.character(target) && any(target[1] == fim4r.targets))
  #stopifnot(is.numeric(supp))
  #stopifnot(is.numeric(zmin)     &&  (zmin >= 0))
  #stopifnot(is.numeric(zmax)     && ((zmax <  0) || (zmax > zmin)))
  #stopifnot(is.character(report))
  #stopifnot(is.numeric(equiv))
  #stopifnot(is.numeric(alpha))
  #stopifnot(is.numeric(smpls))
  #stopifnot(is.numeric(seed))
  # call the C implementation:
  r = .Call("f4r_estpsp", tracts, wgts, target, supp, zmin, zmax,
                          report, equiv, alpha, smpls, seed)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.estpsp()

#-----------------------------------------------------------------------

fim4r.psp2bdr <- function (psp)
{                               # --- extract the decision border
  stopifnot(is.list(psp))
  if ((length(psp) != 3)        # if list of triplets
  || ((length(psp) == 3) && all(sapply(psp, length) == 3)
  &&  (!is.integer(psp[[1]]) || !is.integer(psp[[2]]))))
    psp <- list(as.integer(sapply(psp, "[[", 1)),
                as.integer(sapply(psp, "[[", 2)))
  stopifnot(length(psp[[1]]) == length(psp[[2]]))
  stopifnot(is.numeric(psp[[1]])  && all(psp[[1]] >= 0))
  stopifnot(is.numeric(psp[[2]])  && all(psp[[2]] >= 0))
  b <- tapply(psp[[2]], psp[[1]], max)
  z <- min(psp[[1]])            # extract the decision border
  if (z > 0) {                  # if there are some sizes missing
    b <- c(rep(Inf, z), b+1)    # complete for sizes 0, 1 etc.
    names(b)[1:z] <- 0:(z-1)    # set names of added columns
  }
  for (i in (length(b)-1):1)    # ensure a monotone border
    if (b[i+1] > b[i]) b[i] = b[i+1]
  return(b)                     # return the created border
} # fim4r.psp2bdr()

#-----------------------------------------------------------------------

fim4r.patred <- function (pats, method="S", border=NULL, addis=TRUE)
{                               # --- wrapper for pat. set reduction
  # check the function arguments:
  #stopifnot(is.list(pats))      # check the pattern types
  #stopifnot((all(sapply(pats, function(x) is.integer(x[[1]])))
  #       ||  all(sapply(pats, function(x) is.character(x[[1]]))))
  #       &&  all(sapply(pats, function(x) is.numeric(x[[2]]))))
  #stopifnot(is.character(method)  && any(method[1] == fim4r.reds))
  #stopifnot(is.null(border)       || is.numeric(border))
  #stopifnot(is.logical(addis))
  # call the C implementation:
  r = .Call("f4r_patred", pats, method, border, addis)
  if (is.null(r)) warning("an error occurred, maybe out of memory")
  return(r)                     # check for correct execution
} # fim4r.patred()
