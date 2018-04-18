\name{fim4r}
\alias{fim4r}
\alias{fim4r.fim}
\alias{fim4r.arules}
\alias{fim4r.apriori}
\alias{fim4r.eclat}
\alias{fim4r.fpgrowth}
\alias{fim4r.sam}
\alias{fim4r.relim}
\alias{fim4r.carpenter}
\alias{fim4r.ista}
\alias{fim4r.genpsp}
\alias{fim4r.estpsp}
\alias{fim4r.psp2bdr}
\alias{fim4r.patred}
\alias{fim}
\alias{arules}
\alias{apriori}
\alias{eclat}
\alias{fpgrowth}
\alias{sam}
\alias{relim}
\alias{carpenter}
\alias{ista}
\alias{genpsp}
\alias{estpsp}
\alias{psp2bdr}
\alias{patred}

\title{FIM}
\description{Frequent Item Set Mining and Association Rule Induction}

\usage{
fim4r.fim      (tracts, wgts=NULL, target="s", supp=10.0, zmin=0, zmax=-1,
                report="a", eval="x", agg="x", thresh=10.0, border=NULL)
fim4r.arules   (tracts, wgts=NULL, supp=10.0, conf=80.0, zmin=0, zmax=-1,
                report="aC", eval="x", thresh=10.0, appear=NULL)
fim4r.apriori  (tracts, wgts=NULL, target="s", supp=10.0, conf=80.0,
                zmin=0, zmax=-1, report="a", eval="x", agg="x",
                thresh=10.0, prune=NA, algo="a", mode="",
                border=NULL, appear=NULL)
fim4r.eclat    (tracts, wgts=NULL, target="s", supp=10.0, conf=80.0,
                zmin=0, zmax=-1, report="a", eval="x", agg="x",
                thresh=10.0, prune=NA, algo="a", mode="",
                border=NULL, appear=NULL)
fim4r.fpgrowth (tracts, wgts=NULL, target="s", supp=10.0, conf=80.0,
                zmin=0, zmax=-1, report="a", eval="x", agg="x",
                thresh=10.0, prune=NA, algo="a", mode="",
                border=NULL, appear=NULL)
fim4r.sam      (tracts, wgts=NULL, target="s", supp=10.0, zmin=0, zmax=-1,
                report="a", eval="x", thresh=10.0, algo="a", mode="",
                border=NULL)
fim4r.relim    (tracts, wgts=NULL, target="s", supp=10.0, zmin=0, zmax=-1,
                report="a", eval="x", thresh=10.0, algo="a", mode="",
                border=NULL)
fim4r.carpenter(tracts, wgts=NULL, target="s", supp=10.0, zmin=0, zmax=-1,
                report="a", eval="x", thresh=10.0, algo="a", mode="",
                border=NULL)
fim4r.ista     (tracts, wgts=NULL, target="s", supp=10.0, zmin=0, zmax=-1,
                report="a", eval="x", thresh=10.0, algo="a", mode="",
                border=NULL)
fim4r.genpsp   (tracts, wgts=NULL, target="s", supp=10.0,
                zmin=1, zmax=-1, report="|",
                cnt=1000, surr="s", seed=0, cpus=0)
fim4r.estpsp   (tracts, wgts=NULL, target="s", supp=10.0,
                zmin=1, zmax=-1, report="|", equiv=10000,
                alpha=0.5, smpls=1000, seed=0)
fim4r.psp2bdr  (psp)
fim4r.patred   (pats, method="S", border=NULL, addis=TRUE)
}

\arguments{
\item{tracts}{ a list of transactions, each of which must be an
               integer or string array (and the same data type for
               each list element, that is, either all list elements
               must be integer arrays or all list elements must be
               string arrays). Transactions may be empty, but there
               must be at least one transaction.}
\item{wgts}{   an optional integer array stating the weights or
               multiplicities of the transactions; may be null --
               in this case each transaction has a default weight of 1;
               if given, this array must have the same length as the
               argument \code{tracts}.}
\item{target}{ a string specifying the target type; see
               \code{fim4r.targets} for the available options.
               Note that \code{target="r"} (association rules)
               is possible only for functions \code{fim4r.apriori},
               \code{fim4r.eclat} and \code{fim4r.fpgrowth} and that
               the functions \code{fim4r.carpenter} and
               \code{fim4r.ista} only accept targets "c" (closed
               item sets) and "m" (maximal item sets).}
\item{supp}{   minimum support (number of occurrences).}
\item{zmin}{   minimum pattern size (number of items).}
\item{zmax}{   maximum pattern size (number of items);
               \code{NA} or negative means no limit.}
\item{report}{ a string indicating the values to be reported with
               an item set or an association rule (if letters) or the
               format in which a pattern spectrum is to be reported
               (if "=", "-", "|" or "#").}
\item{algo}{   a string indicating which algorithm variant to use.
               See the variables \code{fim4r.algo.apr} (apriori),
               \code{fim4r.algo.ecl} (eclat),
               \code{fim4r.algo.fpg} (fpgrowth),
               \code{fim4r.algo.sam} (sam),
               \code{fim4r.algo.rem} (relim),
               \code{fim4r.algo.carp} (carpenter),
               \code{fim4r.algo.ista} (ista)
               for the available options.}
\item{mode}{   a string the characters of which indicate operation
               mode options, generally:

               'l': do not use 16-items machine
                    (eclat, fpgrowth, sam, relim);

               'i': do not reorder conditional databases
                    (eclat, fpgrowth);

               'p': do not prune prefix/patricia tree (ista);

               't': do not use a transaction tree (apriori);

               'T': do not use a transaction tree (apriori);

               'u': do not use head-union-tail pruning
                    (eclat, fpgrowth);

               'x': do not use perfect extension pruning
                    (all except ista);

               'y': a-posteriori pruning of infrequent sets (apriori),
                    check extensions horizontally for closed/maximal
                    item sets (eclat);

               'Y': check extensions vertically for closed/maximal
                    item sets (eclat);

               'z': invalidate evaluation below expectation
                    (apriori, eclat, fpgrowth),
                    filter maximal item sets with repository
                    (carpenter, ista).}
\item{border}{ a numeric array with minimum support values per
               pattern size; index 1 refers to size 0, index 2
               to size 1 etc.}
\item{cnt}{    number of surrogate data sets to generate.}
\item{surr}{   surrogate data generation method;
               see \code{fim4r.surrs} for the available options.}
\item{seed}{   seed for random number generator; if 0, the current
               time will be used.}
\item{cpus}{   number of CPUs to use; if 0 or negative, all CPUs
               reported as available by the system will be used.}
\item{psp}{    a pattern spectrum as a list of three numeric vectors
               of equal length; the first states the pattern sizes,
               the second the support values, and the third the
               occurrence frequencies of patterns with the
               corresponding (size,support) signature; the third
               column may be missing in an argument passed to
               \code{fim4r.psp2bdr}.}
\item{equiv}{  equivalent number of surrogate data sets.}
\item{alpha}{  probability dispersion factor.}
\item{smpls}{  number of samples per item set size.}
\item{method}{ pattern set reduction method;
               see \code{fim4r.reds} for the available options.}
\item{addis}{  a logical indicating whether pairwise intersections
               should be added to the pattern set.}
}

\value{
\code{fim4r.fim} returns
a set of item sets as a list of pairs (lists with two elements).
The first element of each pair is an array of integers or strings
(depending on the item type in the argument \code{tracts}), which
are the items of a frequent item set. The second element of each
pair is a numeric array with the values selected by the argument
\code{report} (e.g. support, evaluation etc.). If the argument
\code{report} is "=", "-", "|" or "#", \code{fim4r.fim} returns
a pattern spectrum (see \code{fim4r.genpsp} and \code{fim4r.estpsp}).

\code{fim4r.arules} returns
a set of association rules as a list of triplets (lists with three
elements). The first element of each triplet is an integer or a
string (depending on the item type in the argument \code{tracts}),
which is the consequent/head item of the association rule. The
second element of each triplet is an array of integers or strings
(depending on the item type in the argument \code{tracts}), which
are the items of the antecedent/body of the association rule.
The third element of each triplet is a numeric array with the
values selected in the argument \code{report} (e.g. support,
confidence etc.).

\code{fim4r.apriori} returns
either a list of association rules (see \code{fim4r.arules},
if argument \code{target = "r"} or a list of frequent item sets
(see \code{fim4r.fim}, if argument \code{target != "r"}).
If the argument \code{report} is "=", "-", "|" or "#",
\code{fim4r.fim} returns a pattern spectrum (see \code{fim4r.genpsp}
and \code{fim4r.estpsp}).

\code{fim4r.eclat} returns
either a list of association rules (see \code{fim4r.arules},
if argument \code{target = "r"} or a list of frequent item sets
(see \code{fim4r.fim}, if argument \code{target != "r"}).
If the argument \code{report} is "=", "-", "|" or "#",
\code{fim4r.fim} returns a pattern spectrum (see \code{fim4r.genpsp}
and \code{fim4r.estpsp}).

\code{fim4r.fpgrowth} returns
either a list of association rules (see \code{fim4r.arules},
if argument \code{target = "r"} or a list of frequent item sets
(see \code{fim4r.fim}, if argument \code{target != "r"}).
If the argument \code{report} is "=", "-", "|" or "#",
\code{fim4r.fim} returns a pattern spectrum (see \code{fim4r.genpsp}
and \code{fim4r.estpsp}).

\code{fim4r.sam} returns
a list of frequent item sets (see \code{fim4r.fim}).
If the argument \code{report} is "=", "-", "|" or "#",
\code{fim4r.fim} returns a pattern spectrum (see \code{fim4r.genpsp}
and \code{fim4r.estpsp}).

\code{fim4r.relim} returns
a list of frequent item sets (see \code{fim4r.fim}).
If the argument \code{report} is "=", "-", "|" or "#",
\code{fim4r.fim} returns a pattern spectrum (see \code{fim4r.genpsp}
and \code{fim4r.estpsp}).

\code{fim4r.carpenter} returns
a list of frequent item sets (see \code{fim4r.fim}).
If the argument \code{report} is "=", "-", "|" or "#",
\code{fim4r.fim} returns a pattern spectrum (see \code{fim4r.genpsp}
and \code{fim4r.estpsp}).

\code{fim4r.ista} returns
a list of frequent item sets (see \code{fim4r.fim}).
If the argument \code{report} is "=", "-", "|" or "#",
\code{fim4r.fim} returns a pattern spectrum (see \code{fim4r.genpsp}
and \code{fim4r.estpsp}).

\code{coco4r.genpsp} and \code{coco4r.estpsp} return
a pattern spectrum, either as a list with three columns
(for \code{report = "#"} or \code{report = "|"}, default)
or as a list of numeric vectors with three elements
(for \code{report = "="} or \code{report = "-"}).
The first column/element is the pattern size (number of items),
the second column/element is the pattern support (number of occurrences)
the third column/element is the number of patterns with the signature
(size,support) given by the first two columns/elements.

\code{fim4r.psp2bdr} returns
a pattern border as a numeric array with named elements. The names
are the pattern sizes, the elements the corresponding minimum support
values a pattern of that size has to reach in order to be considered
significant. This list can be passed to the function
\code{coco4r.coconad} for the argument \code{border}.

\code{fim4r.patred} returns
a reduced pattern set in the same form as \code{coco4r.coconad}.
}

\author{Christian Borgelt \email{christian@borgelt.net}}

\examples{
library(fim4r)                  # load the FIM for R library
data(tracts)                    # load the example transactions
tracts <- tapply(as.character(tracts[,2]), tracts[,1], c)
                                # reshape into required input format
isets  <- fim4r.fim(tracts)     # find frequent item sets
rules  <- fim4r.arules(tracts)  # find association rules
}

\keyword{transaction}
\keyword{frequent pattern mining}
\keyword{frequent item set mining}
\keyword{association rule}
