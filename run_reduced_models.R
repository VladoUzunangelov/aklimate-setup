#!/usr/bin/env Rscript

##copyright (c) 2019 Vlado Uzunangelov
homeDir <- "/home/ubuntu/remote_dirs/courtyard/pstore_migration/home_dir"

dataDir <- "/home/ubuntu/projects/THYM_AKLIMATE_33_environment/"
target.dir <- "/home/ubuntu/projects/THYM_AKLIMATE_33_environment/models/"

    library(doParallel)
    library(foreach)
    ##library(doRNG)
    registerDoParallel(cores=15)

    library(abind)
    library(ranger)
    library(dplyr)
    library(caret)
    library(ROCR)
    library(pracma)
    source("/home/ubuntu/repos/tcga_scripts/utils.R",chdir=TRUE)
    source("/home/ubuntu/repos/junkle/junkle-utils.R",chdir=TRUE)
source("/home/ubuntu/repos/junkle/junkle.R",chdir=TRUE)
source("/home/ubuntu/repos/Spicer/Spicer.R",chdir=TRUE)
source("/home/ubuntu/repos/Spicer/Spicer-classify.R",chdir=TRUE)
    
    
    ############################################################################


    #############################

labels <-as.matrix(read.delim(paste0(dataDir,"/labels.tsv"),check.names=F,stringsAsFactors=F,header=F,row.names=1))
labels<-factor(labels[,1])

library(FDb.InfiniumMethylation.hg19)



##MIR is also there - about 800 MIRs
##right now probably best to not include them - I have a scheme to indluce them by putting
##a MIR in every pathway that has one of its targets, but it didn't work well when
##I tested it

groups<-list(nomir=c("MUTA","CNVR","METH","GEXP"))

dat<-read.delim(paste0(dataDir,"/data/THYM_20190424.tsv"),header=T,row.names=1,check.names=FALSE)
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


    labels <- labels[rownames(dat)]
p1 <- readSetList(paste0(homeDir,"/mkl/data/feature_sets/pathcomm_pathways_cleaned"))
p1 <- p1[!grepl("[[:digit:]]_SMPDB$|Z_SMPDB$",names(p1))]
p2 <- readSetList(paste0(homeDir,"/mkl/data/tmp_awg/LIHC/genomic_position_sets.listt"))
p3 <- readSetList(paste0(homeDir,"/mkl/data/feature_sets/genesigdb_human.tab"))
p4 <- readSetList(paste0(homeDir,"/mkl/data/feature_sets/msigdb_c2_c5_no_c2_cp.tab"))

pathways <- c(p1,p2,p3,p4) 
pathways <- pathways[sapply(pathways,length)<1000]
    

splits<- as.matrix(read.delim(paste0(dataDir,"/data/THYM_CVfolds_5FOLD_20190328.tsv"),check.names=F,stringsAsFactors=F,header=TRUE,row.names=1))
splits<-splits[,-1]

tasks<-colnames(splits)[1:25]


#######################################################


#
z1 <- foreach(i=iter(tasks))%dopar%{
    load(paste0(target.dir,"/",i,"_junkle_final_model.RData"))
    imps <- rank.features.jklm(jklm)
    write.df(data.frame(importance=imps),"features",paste0(target.dir,"/",i,"_aklimate_multiclass_feature_importance.tab"))
}

cutoffs <- c(5,10,20,50,100,200,500,1000,1500)
reps.list <- lapply(1:5,function(x) tasks[seq(5*(x-1)+1,5*x,1)])


acc.reduced <- foreach(i=iter(reps.list),.combine=rbind)%dopar%{
    foreach(k=iter(cutoffs),.combine=c)%do%{
        preds <- foreach(j=iter(i),.combine=c)%do%{
            idx.train<-rownames(splits)[splits[,j]==0]
            idx.test<-setdiff(rownames(splits),idx.train)

            


            imps <- read.delim(paste0(target.dir,"/",j,"_aklimate_multiclass_feature_importance.tab"),header=TRUE,row.names=1)
            dat <- mlr::createDummyFeatures(dat)
            k.adj<-min(k,dim(imps)[1])
            rf <- ranger(data=cbind(data.frame(labels=labels[idx.train]),
                             dat[idx.train,rownames(imps)[1:k.adj],drop=FALSE]),
                            dependent.variable.name="labels",
                            always.split.variables=NULL,
                            classification = TRUE,
                            sample.fraction = 0.5,
                            num.trees=3000,
                            mtry=ceiling(k.adj/5),
                            min.node.size=1,
                            case.weights=NULL,
                            num.threads=3,
                            probability=TRUE,
                            respect.unordered.factors = FALSE,
                            importance="none",
                            write.forest=TRUE,
                            keep.inbag=TRUE,
                            replace=FALSE)
            save(rf,file=paste0(target.dir,"/",j,"_cutoff_",k,"_rf_reduced_model.RData"))
            rf.preds <- predict(rf,dat[idx.test,rownames(imps)[1:k.adj]])$predictions         
            rownames(rf.preds) <- idx.test
            save(rf.preds,file=paste0(target.dir,"/",j,"_cutoff_",k,"_rf_reduced_model_predictions.RData"))

                                         confM <- caret::confusionMatrix(
                             factor(apply(rf.preds[idx.test,],1,which.max),
                                    levels=levels(labels)),labels[idx.test])
            save(confM,file=paste0(target.dir,"/",j,"_cutoff_",k,"_rf_reduced_model_stats.RData"))

        mean(unname(confM$byClass[,"Balanced Accuracy"]))        
        ##mean(unname(confM$overall["Accuracy"]))    
        }

        mean(preds)

    }
}

colnames(acc.reduced) <- cutoffs
write.df(data.frame(acc.reduced,check.names=FALSE),"repetition",paste0(target.dir,"/balanced_accuracy_reduced_feature_sets.tab"))

######################################################
acc.reduced <- foreach(i=iter(reps.list),.combine=rbind)%dopar%{
    foreach(k=iter(cutoffs),.combine=c)%do%{
        preds <- foreach(j=iter(i),.combine=c)%do%{
            load(paste0(target.dir,"/",j,"_cutoff_",k,"_rf_reduced_model_stats.RData"))
        mean(unname(confM$overall["Accuracy"]))    
        }

        mean(preds)

    }
}

colnames(acc.reduced) <- cutoffs
write.df(data.frame(acc.reduced,check.names=FALSE),"repetition",paste0(target.dir,"/accuracy_reduced_feature_sets.tab"))
