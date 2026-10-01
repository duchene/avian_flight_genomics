# Assumes an exponential distribution in branch lengths
require(MASS)

trimLongTipBrs <- function(tr, pval = 0.01){
	exprate <- fitdistr(tr$edge.length, "exponential")$estimate
	thresh <- qexp(pval, exprate, lower.tail = F)
	if(all(tr$edge.length < thresh)){
		print("All branches have acceptable length. Returning the input tree.")
		return(tr)
	} else {
		badParentNodes <- tr$edge[which(tr$edge.length > thresh), 2]
		if(all(badParentNodes > Ntip(tr))){
			print("Only internal branches are long. Returning the input tree.")
			return(tr)
		} else {
			badTips <- badParentNodes[which(badParentNodes <= Ntip(tr))]
			print(paste0("Long tips include ", tr$tip.label[badTips], ". Returning pruned tree.", collapse = " "))
			tr <- drop.tip(tr, badTips)
			return(tr)
		}
	}
}