#!/usr/bin/env Rscript


##To be run on 10.50.100.33
## in /home/ubuntu/projects/THYM_AKLIMATE/scripts

ncpus <- 30
library(doParallel)
library(foreach)
##library(doRNG)
registerDoParallel(cores=ncpus)

library(abind)
library(ranger)
library(dplyr)
library(caret)
library(ROCR)
##library(mlr)
source("/home/ubuntu/repos/tcga_scripts/utils.R",chdir=TRUE)
source("/home/ubuntu/repos/junkle/junkle-utils.R",chdir=TRUE)
source("/home/ubuntu/repos/junkle/junkle.R",chdir=TRUE)
source("/home/ubuntu/repos/Spicer/Spicer.R",chdir=TRUE)
source("/home/ubuntu/repos/Spicer/Spicer-classify.R",chdir=TRUE)

##need to install mySQL first
##sudo apt-get install libmariadbclient-dev
##then in R
##source("https://bioconductor.org/biocLite.R")
##biocLite("FDb.InfiniumMethylation.hg19")

library(FDb.InfiniumMethylation.hg19)



##MIR is also there - about 800 MIRs
##right now probably best to not include them - I have a scheme to indluce them by putting
##a MIR in every pathway that has one of its targets, but it didn't work well when
##I tested it

suffs<-list(nomir=c("MUTA","CNVR","METH","GEXP"))

dat<-read.delim("../data/THYM_20190424.tsv",header=T,row.names=1,check.names=FALSE)
dat<-dat[,-1]
##load("../data/meth_mappers.RData")
##move this to a makefile or something - take a couple of minutes to compute and should only be done once
hm450 <- get450k()
meth.mapper <- getNearestGene(hm450)

updated.names<-foreach(i=iter(colnames(dat)),.combine=c)%dopar%{
    if(grepl("METH",i,ignore.case = TRUE)) {
        cg<-strsplit(i,":")[[1]][4]
        gsub(paste0(cg,"::"),paste0(cg,":",meth.mapper[cg,4],":"),i)
    } else {
        i
    }
}

colnames(dat)<-updated.names

homeDir <- "/home/ubuntu/remote_dirs/courtyard/pstore_migration/home_dir"
workDir<-"../models/"


p1 <- readSetList(paste0(homeDir,"/mkl/data/feature_sets/pathcomm_pathways_cleaned"))
p1 <- p1[!grepl("[[:digit:]]_SMPDB$|Z_SMPDB$",names(p1))]
p2 <- readSetList(paste0(homeDir,"/mkl/data/tmp_awg/LIHC/genomic_position_sets.listt"))
p3 <- readSetList(paste0(homeDir,"/mkl/data/feature_sets/genesigdb_human.tab"))
p4 <- readSetList(paste0(homeDir,"/mkl/data/feature_sets/msigdb_c2_c5_no_c2_cp.tab"))

pathways <- c(p1,p2,p3,p4) 
pathways <- pathways[sapply(pathways,length)<1000]

##probably a good idea to sanitize the name a bit
##names(pathways)<-gsub('[^\\w\\s]',"_",names(pathways),perl=TRUE)

labels <-as.matrix(read.delim("../labels.tsv",check.names=F,stringsAsFactors=F,header=F,row.names=1))
labels<-factor(labels[,1])

splits<- as.matrix(read.delim("../data/THYM_CVfolds_5FOLD_20190328.tsv",check.names=F,stringsAsFactors=F,header=TRUE,row.names=1))
splits<-splits[,-1]

tasks<-colnames(splits)

#you can change the seeds - this is just to match initial runs
seeds<-11*(1:length(tasks))
names(seeds)<- tasks

worker.f<-function(tasks) {

res<- foreach(i=iter(tasks))%do% {

    set.seed(seeds[i],kind = "L'Ecuyer-CMRG")
        
    idx.train<-rownames(splits)[splits[,i]==0]
    idx.test<-setdiff(rownames(splits),idx.train)

       jklm <- junkle(dat[idx.train,],suffs,labels[idx.train],pathways,NULL,list(ttype="multiclass",bin.perf=c("bacc"),importance="permutation",min.nfeat=15,ntree=1000,sample.frac=0.5,replace=FALSE,weights=NULL,oob.cv=data.frame(min.node.prop=0.01,mtry.prop=0.25,ntree=500)),list(topn=5,subsetCV=TRUE,lamb=c(-8,3),cvlen=200,type="probability"),FALSE,TRUE)
        
        save(jklm,file=paste0(workDir,"/",i,"_junkle_final_model.RData"))
        
        jklm.preds <- predict.junkle(jklm,dat[c(idx.train,idx.test),],pathways,NULL,FALSE)$preds
    jklm.preds<-apply(jklm.preds,1,function(x) colnames(jklm.preds)[which.max(x)])

    confM <- caret::confusionMatrix(
                        factor(jklm.preds,
                               levels=levels(labels)),labels[idx.test])
    bacc <- mean(unname(confM$byClass[,"Balanced Accuracy"]))        



        save(jklm.preds,confM,bacc,file=paste0(workDir,"/",i,"_junkle_final_model_stats_preds.RData"))
        
    return(bacc)


}

    return(res)
}

results<-worker.f(tasks[1:25])
names(results)<-tasks[1:25]

results <- unlist(results)

names(results) <- as.character(tasks[1:25])

write.df(data.frame(bacc=results),"split",paste0(workDir,"/splits_balanced_accuracy.tab"))
