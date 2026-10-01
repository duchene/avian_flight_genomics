concatenate <- function(alslist, AA = F){
            allsps <- unique(unlist(lapply(alslist, function(x) rownames(x))))
            for(i in 1:length(alslist)){
                  missingsp <- allsps[which(!allsps %in% rownames(alslist[[i]]))]
                  if(length(missingsp) > 0){
                        if(AA) missingmat <- as.AAbin(matrix("-", nrow = length(missingsp), ncol = ncol(alslist[[i]]))) else missingmat <- as.DNAbin(matrix("-", nrow = length(missingsp), ncol = ncol(alslist[[i]])))
                        rownames(missingmat) <- missingsp
                        alslist[[i]] <- rbind(alslist[[i]], missingmat)[allsps, ]
                  } else {
                        alslist[[i]] <- alslist[[i]][allsps, ]
                  }
            }
            allal <- do.call("cbind", alslist)
            return(allal)
}
