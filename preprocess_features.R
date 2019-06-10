#!/usr/bin/env Rscript
source("~/repos/junkle/junkle-utils.R",chdir=TRUE)
dataDir <- "/pod/pstore/users/uzunangelov/mkl/data/tmp_awg/LIHC/"
groups <- c("exp","mirna","methylation","mutation","cnv","cn_regions","lncrna")

##these are the names of samples you will use
new_universe <- scan(paste0("all_samples_334"),what=character())
prefix <- "20180731_"
exp <- read.delim(paste0(dataDir,"exp.tab"),header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
exp <- log(exp+1)
means <- colMeans(exp)
sds <- apply(exp,2,sd)

##exclude genes that are in the bottom quartile by mean or sd
 keep <- intersect(names(means)[means>quantile(means,probs=0.25)],names(sds)[sds>quantile(sds,probs=0.25)])
 exp <- exp[,keep]


mirna <- read.delim(paste0(dataDir,"mirna.tab"),header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)


lncrna <- read.delim(paste0(dataDir,"lncrna.tab"),header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

##cnv is noisier - remove bottom 50%
cnv <- read.delim(paste0(dataDir,"cnv.tab"),header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
means <- colMeans(abs(cnv))
sds <- apply(cnv,2,sd)
 keep <- intersect(names(means)[means>quantile(means,probs=0.5)],names(sds)[sds>quantile(sds,probs=0.5)])
cnv <- cnv[,keep]


##this section is options - it add CN features that are not gene specific but rather copy number events
##you can change the epsilon to get more or less fine-grained events
library(iClusterPlus)
library(cluster)
regions <- read.delim(paste0(dataDir,"cnv_regions.tab"),header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
cn_regions <- CNregions(regions,epsilon=0.0025)
##common <- intersect(rownames(cn_regions),common)
write.df(data.frame(cn_regions,check.names=FALSE),"samples",paste0(dataDir,"cnv_regions_reduced.tab"))

methylation <- read.delim("synapse/20180613_LIHC.CHOL380_DNA.methylation.data.tsv",row.names=1,header=TRUE,check.names=FALSE)
methylation <- as.data.frame(t(methylation))
rownames(methylation) <- paste0("TCGA-",rownames(methylation))
methylation <- methylation[new_universe,]

mutation <- read.delim(paste0(dataDir,"cravat_bop/cravat_mutations.tab"),header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

counts <- colSums(mutation>0)
##use this only if you have continuous values like from VEST/CRAVAT
##means <- sapply(1:ncol(mutation),function(x) if(counts[x]>0) sum(as.numeric(mutation[,x]))/counts[x] else 0)
##keep <- intersect(names(means)[means>0.5],names(counts)[counts>5])
keep <- names(counts)[counts>5]
mutation <- mutation[,keep]


## This is all about cleaning up clinical data - left it just for reference!
##+1 for rownames
## colClasses <- rep("numeric",6)
## colClasses[1] <- "character"
## colClasses[2:3] <- "factor"
## clinical <- read.delim(paste0(dataDir,"/synapse/FreezeSamples_MasterTable_20180718_selected.tsv"),header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE,na.strings=c("","[Not Available]","[Not Evaluated]","[Discrepancy]","[Unknown]"),colClasses=colClasses)
## colnames(clinical) <- gsub('([[:punct:]])|\\s+','_',colnames(clinical))
## colnames(clinical) <- paste0(colnames(clinical),"_general")
## rownames(clinical) <- gsub("[AB]$","",rownames(clinical))

## clinical <- clinical[new_universe,]
## clinical <- as.data.frame(lapply(clinical,function(x) {
##     if(sum(is.na(x))>0) {
##         if(is.factor(x)) {
##             res <- addNA(x)
##             ##levels(res) <- c(levels(x),"unknown")
##         } else {
##             res <- x
##             res[is.na(x)] <- mean(x,na.rm=TRUE)
##         }
##     } else {
##         res <- x
##     }
##     return(res) }))
## rownames(clinical) <- new_universe



dat.comb <- collapse.data(list(general=clinical[new_universe,],exp=exp[new_universe,],cnv=cnv[new_universe,],cn_regions=cn_regions[new_universe,],mirna=mirna[new_universe,],methylation=methylation[new_universe,],mutation=mutation[new_universe,],lncrna=lncrna[new_universe,]))

save(dat.comb,file=paste0(prefix,"combined_features_filtered_mean_var.RData"))

#####################################################################

##THIS IS THE WAY TO BIN CONTINUOUS DATA
## THE WAY IT WORKS - nbr=number of bins
##idx - indices of training data
##if idx == rownames(dat) (default) then all data is considered simultaneously for binning
## if idx!=rownames(dat) the training set is idx, which is what you use to construct the bins,
##the test set is rownames(dat)-idx, which uses the bins constructed on idx to bin data



dat.comb.binned <- foreach(i=iter(groups[!grepl("methylation|mutation",groups)]),.combine=cbind)%do%{
    quantize.data(dat.comb[new_universe,grepl(paste0("_",i,"$"),colnames(dat.comb))],nbr=5,idx=new_universe)
}

colnames(methylation) <- paste0(colnames(methylation),"_methylation")
colnames(mutation) <- paste0(colnames(mutation),"_mutation")
dat.comb.binned <- cbind(dat.comb.binned,
                         clinical[new_universe,],
                         methylation[new_universe,],
                         mutation[new_universe,])
save(dat.comb.binned,file=paste0(prefix,"combined_features_binned_5_filtered_mean_var.RData"))


###########################################################################

