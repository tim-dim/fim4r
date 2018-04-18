#!/usr/bin/Rscript

library(fim4r)                  # load the FIM for R library

#-----------------------------------------------------------------------

showpats <- function (pats)
{                               # print found frequent item sets
  for (i in 1:length(pats))
    cat(sort(pats[[i]][[1]]), sprintf("(%d)\n", pats[[i]][[2]]))
  cat(sprintf("%d pattern(s)\n", length(pats)))
}  # showpats()

#-----------------------------------------------------------------------

showrules <- function (rules)
{                               # print found association rules
  for (i in 1:length(rules)) {
    cat(rules[[i]][[1]], "<-", sort(rules[[i]][[2]]),
        sprintf("(%g,",  rules[[i]][[3]][[1]]),
        sprintf("%g)\n", rules[[i]][[3]][[2]]))
  }
  cat(sprintf("%d rule(s)\n", length(rules)))
}  # showrules()

#-----------------------------------------------------------------------

data(tracts)                    # load the example transactions
tracts <- tapply(as.character(tracts[,2]), tracts[,1], c)
apps   <- list(c("","b"),c("a","c"))

cat("------------------------------------------------------------\n")
cat("fim\n")
cat("------------------------------------------------------------\n")
showpats(fim4r.fim(tracts, supp=-2))

cat("------------------------------------------------------------\n")
cat("apriori\n")
cat("------------------------------------------------------------\n")
showpats(fim4r.apriori(tracts, supp=-2))

cat("------------------------------------------------------------\n")
cat("eclat\n")
cat("------------------------------------------------------------\n")
showpats(fim4r.eclat(tracts, supp=-2))

cat("------------------------------------------------------------\n")
cat("fpgrowth\n")
cat("------------------------------------------------------------\n")
showpats(fim4r.fpgrowth(tracts, supp=-2))

cat("------------------------------------------------------------\n")
cat("sam\n")
cat("------------------------------------------------------------\n")
showpats(fim4r.sam(tracts, supp=-2))

cat("------------------------------------------------------------\n")
cat("relim\n")
cat("------------------------------------------------------------\n")
showpats(fim4r.relim(tracts, supp=-2))

cat("------------------------------------------------------------\n")
cat("carpenter\n")
cat("------------------------------------------------------------\n")
showpats(fim4r.carpenter(tracts, supp=-2))

cat("------------------------------------------------------------\n")
cat("ista\n")
cat("------------------------------------------------------------\n")
showpats(fim4r.ista(tracts, supp=-2))

cat("------------------------------------------------------------\n")
cat("arules\n")
cat("------------------------------------------------------------\n")
showrules(fim4r.arules(tracts, supp=-2, appear=apps))

cat("------------------------------------------------------------\n")
cat("apriori\n")
cat("------------------------------------------------------------\n")
showrules(fim4r.apriori(tracts, target="r", supp=-2,
                        report="aC", appear=apps))

cat("------------------------------------------------------------\n")
cat("eclat\n")
cat("------------------------------------------------------------\n")
showrules(fim4r.eclat(tracts, target="r", supp=-2,
                      report="aC", appear=apps))

cat("------------------------------------------------------------\n")
cat("fpgrowth\n")
cat("------------------------------------------------------------\n")
showrules(fim4r.fpgrowth(tracts, target="r", supp=-2,
                         report="aC", appear=apps))
