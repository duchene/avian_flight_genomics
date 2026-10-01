rtt.cov <- function(tre){
        if(!is.rooted(tre)) tre <- midpoint(tre)
        rtts <- pathnode(tre)[[1]]
        rttcov <- sd(rtts, na.rm = T) / mean(rtts, na.rm = T)
        return(rttcov)
}
