#!/usr/bin/env Rscript

##copyright (c) 2019 Vlado Uzunangelov

dataDir <- "/home/ubuntu/remote_dirs/courtyard/pstore_migration/home_dir/mkl/data/tmp_awg/LIHC"
target.dir <- getwd()
tasks <- 1:25

    library(doParallel)
    library(foreach)
    ##library(doRNG)
    registerDoParallel(cores=12)

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

    labels <-as.matrix(read.delim(paste0(dataDir,"/iclust_membership"),check.names=F,stringsAsFactors=F,header=T,row.names=1))
    load(paste0(dataDir,"/20180731_combined_features_binned_5_filtered_mean_var.RData"))
    load(paste0(dataDir,"/path_comm_pathways_methylation_mirna_extended_plus_positional_sets.RData"))

    dat <- dat.comb.binned
    labels <- labels[rownames(dat),]
    
        pathways <- path.comm.pathways.positional.sets
        pathways <- pathways[!grepl("[[:digit:]]_SMPDB$|Z_SMPDB$",names(pathways))]

        groups <- c("exp","methylation","mutation","cnv","cn_regions","lncrna")

    new_folds <- readSetList(paste0(dataDir,"/cv_folds.listt"))
        new_universe <- scan(paste0(dataDir,"/all_samples"),what=character())
        
lbls <- factor(labels)
names(lbls) <- names(labels)
        lbls <- lbls[new_universe]
        dat <- dat[new_universe,]

#######################################################


#
z1 <- foreach(i=iter(tasks))%do%{
    load(paste0(target.dir,"/split_",i,"_aklimate_multiclass_model.RData"))
    imps <- rank.features.jklm(jklm)
    write.df(data.frame(importance=imps),"features",paste0(target.dir,"/split_",i,"_aklimate_multiclass_feature_importance.tab"))
}

cutoffs <- c(5,10,20,50,100,200,500,1000,1500)
reps.list <- lapply(1:25,function(x) seq(5*(x-1)+1,5*x,1))


acc.reduced <- foreach(i=iter(reps.list),.combine=rbind)%dopar%{
    foreach(k=iter(cutoffs),.combine=c)%do%{
        preds <- foreach(j=iter(i),.combine=c)%do%{
            idx.train <- new_folds[[j]]
            idx.test <- setdiff(new_universe,idx.train)
            


            imps <- read.delim(paste0(target.dir,"/split_",j,"_aklimate_multiclass_feature_importance.tab"),header=TRUE,row.names=1)
            dat <- mlr::createDummyFeatures(dat)
            rf <- ranger(data=cbind(data.frame(labels=lbls[idx.train]),
                             dat[idx.train,rownames(imps)[1:k],drop=FALSE]),
                            dependent.variable.name="labels",
                            always.split.variables=NULL,
                            classification = TRUE,
                            sample.fraction = 0.5,
                            num.trees=3000,
                            mtry=ceiling(k/5),
                            min.node.size=1,
                            case.weights=NULL,
                            num.threads=3,
                            probability=TRUE,
                            respect.unordered.factors = FALSE,
                            importance="none",
                            write.forest=TRUE,
                            keep.inbag=TRUE,
                            replace=FALSE)
            save(rf,file=paste0(target.dir,"/split_",j,"_cutoff_",k,"_rf_reduced_model.RData"))
            rf.preds <- predict(rf,dat[idx.test,rownames(imps)[1:k]])$predictions         
            rownames(rf.preds) <- idx.test
            save(rf.preds,file=paste0(target.dir,"/split_",j,"_cutoff_",k,"_rf_reduced_model_predictions.RData"))

                                         confM <- caret::confusionMatrix(
                             factor(apply(rf.preds[idx.test,],1,which.max),
                                    levels=levels(lbls)),lbls[idx.test])
            save(confM,file=paste0(target.dir,"/split_",j,"_cutoff_",k,"_rf_reduced_model_stats.RData"))

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
            load(paste0(target.dir,"/split_",j,"_cutoff_",k,"_rf_reduced_model_stats.RData"))
        mean(unname(confM$overall["Accuracy"]))    
        }

        mean(preds)

    }
}

colnames(acc.reduced) <- cutoffs
write.df(data.frame(acc.reduced,check.names=FALSE),"repetition",paste0(target.dir,"/accuracy_reduced_feature_sets.tab"))

#######################################################################
##same as above, but with feature selection proportions derived from aklimate weights
cutoffs <- c(5,10,20,50,100,200,500,1000,1500)
reps.list <- lapply(1:25,function(x) seq(5*(x-1)+1,5*x,1))


acc.reduced.akl <- foreach(i=iter(reps.list),.combine=rbind)%dopar%{
    foreach(k=iter(cutoffs),.combine=c)%do%{
        preds <- foreach(j=iter(i),.combine=c)%do%{
            idx.train <- new_folds[[j]]
            idx.test <- setdiff(new_universe,idx.train)
            


            imps <- read.delim(paste0(target.dir,"/split_",j,"_aklimate_multiclass_feature_importance.tab"),header=TRUE,row.names=1)
            dat <- mlr::createDummyFeatures(dat)
            
            rf <- ranger(data=cbind(data.frame(labels=lbls[idx.train]),
                             dat[idx.train,rownames(imps)[1:k],drop=FALSE]),
                            dependent.variable.name="labels",
                            always.split.variables=NULL,
                            classification = TRUE,
                            sample.fraction = 0.5,
                            num.trees=3000,
                            mtry=ceiling(k/5),
                            min.node.size=1,
                         case.weights=NULL,
                         split.select.weights = imps[1:k,1]/sum(imps[1:k,1]), 
                            num.threads=3,
                            probability=TRUE,
                            respect.unordered.factors = FALSE,
                            importance="none",
                            write.forest=TRUE,
                            keep.inbag=TRUE,
                            replace=FALSE)
            save(rf,file=paste0(target.dir,"/split_",j,"_cutoff_",k,"_rf_reduced_model_aklimate_weighting.RData"))
            rf.preds <- predict(rf,dat[idx.test,rownames(imps)[1:k]])$predictions         
            rownames(rf.preds) <- idx.test
            save(rf.preds,file=paste0(target.dir,"/split_",j,"_cutoff_",k,"_rf_reduced_model_aklimate_weighting_predictions.RData"))

                                         confM <- caret::confusionMatrix(
                             factor(apply(rf.preds[idx.test,],1,which.max),
                                    levels=levels(lbls)),lbls[idx.test])
            save(confM,file=paste0(target.dir,"/split_",j,"_cutoff_",k,"_rf_reduced_model_aklimate_weighting_stats.RData"))

        mean(unname(confM$byClass[,"Balanced Accuracy"]))        
        ##mean(unname(confM$overall["Accuracy"]))    
        }

        mean(preds)

    }
}

colnames(acc.reduced.akl) <- cutoffs
write.df(data.frame(acc.reduced.akl,check.names=FALSE),"repetition",paste0(target.dir,"/balanced_accuracy_reduced_feature_sets_aklimate_weighting.tab"))





###########################################
####
overall.acc <- foreach(i=iter(reps.list),.combine=c)%dopar%{
    scores <- foreach(k=iter(i),.combine=c)%do%{

        load(paste0(target.dir,"/split_",k,"_predictions.RData"))
            idx.train <- new_folds[[k]]
            idx.test <- setdiff(new_universe,idx.train)

                                         confM <- caret::confusionMatrix(
                             factor(jklm.preds,
                                    levels=levels(lbls)),lbls[idx.test])

        mean(unname(confM$byClass[,"Balanced Accuracy"]))        
        ##mean(unname(confM$overall["Accuracy"]))    

    }
    mean(scores)
}

write.df(data.frame(overall=overall.acc),"repetition",paste0(target.dir,"/balanced_accuracy_full_aklimate.tab"))
##write.df(data.frame(overall=overall.acc),"repetition",paste0(target.dir,"/accuracy_full_aklimate.tab"))
###########################################################
##########

confM.list <- foreach(i=iter(reps.list))%dopar%{
    scores <- foreach(k=iter(i))%do%{

        load(paste0(target.dir,"/split_",k,"_predictions.RData"))
            idx.train <- new_folds[[k]]
            idx.test <- setdiff(new_universe,idx.train)

                                         confM <- caret::confusionMatrix(
                             factor(jklm.preds,
                                    levels=levels(lbls)),lbls[idx.test])

        confM$table
    }

    confM.avg <- Reduce('+',scores)
    confM.avg <- as.data.frame.matrix(confM.avg)
    rownames(confM.avg) <- colnames(confM.avg) <- paste0("iClust_",colnames(confM.avg))
    confM.avg

}

save(confM.list,file=paste0(target.dir,"/confusion_matrix_list.RData"))

confM.avg <- Reduce('+',confM.list)/length(confM.list)
confM.avg <- as.data.frame.matrix(confM.avg)

write.df(confM.avg,"",paste0(target.dir,"/confusion_matrix_average.tab"))


########################################################################3
#

props.list <- foreach(i=iter(reps.list))%dopar%{
        props <- foreach(j=iter(i))%do%{
            imps <- read.delim(paste0(target.dir,"/split_",j,"_aklimate_multiclass_feature_importance.tab"),header=TRUE,row.names=1)
            imps <- setNames(imps[,1],rownames(imps))
            coffs <- c(cutoffs,length(imps))
            types.prop <- foreach(k=iter(coffs),.combine=cbind)%do%{
                out <- rank.importance.type(groups,imps[1:k],intron="_")
                out[groups]
            }
            colnames(types.prop) <- c(as.character(coffs[-length(coffs)]),"full")
            types.prop
    }
       Reduce("+",props)/length(props)
}
avg.props <- Reduce("+",props.list)/length(props.list)
avg.props <- avg.props[order(avg.props[,"full"],decreasing=TRUE),]

write.df(data.frame(avg.props,check.names=FALSE),"data_type",paste0(target.dir,"/data_type_contributions.tab"))

#########################################################3

t1 <- foreach(i=iter(reps.list),.combine=cbind)%dopar%{
    t2 <- foreach(j=iter(i),.combine=c)%do%{
        idx.train <- new_folds[[j]]
        idx.test <- setdiff(new_universe,idx.train)
        load(paste0(target.dir,"/split_",j,"_predictions.RData"))
        t3 <- as.numeric(jklm.preds)
        names(t3) <- names(jklm.preds)
        t3
    }
    t2 <- t2[new_universe]
    t2
}

write.df(t1,"sample",paste0(target.dir,"/aklimate_cv_individual_calls_334_samples.tab"))

t1 <- foreach(i=iter(reps.list),.combine=cbind)%dopar%{
    t2 <- foreach(j=iter(i),.combine=c)%do%{
        idx.train <- new_folds[[j]]
        idx.test <- setdiff(new_universe,idx.train)
        load(paste0(target.dir,"/split_",j,"_predictions.RData"))
        jklm.preds.bin <- as.numeric(as.numeric(jklm.preds[idx.test])==as.numeric(lbls[idx.test]))
        names(jklm.preds.bin) <- idx.test
        jklm.preds.bin
    }
    t2 <- t2[new_universe]
    t2
}

stats <- cbind(ncol(t1),rowMeans(t1))
colnames(stats) <- c("count","proportion")
write.df(stats,"sample",paste0(target.dir,"/aklimate_cv_per_sample_accuracy_334_samples.tab"))

######################################################################################

cf <- 20
##cf <- 50
t1 <- foreach(i=iter(reps.list),.combine=cbind)%dopar%{
    t2 <- foreach(j=iter(i),.combine=c)%do%{
        load(paste0(target.dir,"/split_",j,"_cutoff_",cf,"_rf_reduced_model_predictions.RData"))
        t3 <- apply(rf.preds,1,function(x) as.numeric(colnames(rf.preds)[which.max(x)]))
        names(t3) <- rownames(rf.preds)
        t3
    }
    t2 <- t2[new_universe]
    t2
}

write.df(data.frame(t1),"sample",paste0(target.dir,"/rf_reduced_model_cutoff_",cf,"_individual_calls_334_samples.tab"))

t1 <- foreach(i=iter(reps.list),.combine=cbind)%dopar%{
    t2 <- foreach(j=iter(i),.combine=c)%do%{
        idx.train <- new_folds[[j]]
        idx.test <- setdiff(new_universe,idx.train)
        load(paste0(target.dir,"/split_",j,"_cutoff_",cf,"_rf_reduced_model_predictions.RData"))
        t3 <- apply(rf.preds,1,function(x) as.numeric(colnames(rf.preds)[which.max(x)]))
        names(t3) <- rownames(rf.preds)
        t3 <- as.numeric(t3[idx.test]==as.numeric(lbls[idx.test]))
        names(t3) <- idx.test
        t3
    }
    t2 <- t2[new_universe]
    t2
}
stats <- cbind(ncol(t1),rowMeans(t1))
colnames(stats) <- c("count","proportion")
write.df(stats,"sample",paste0(target.dir,"/rf_reduced_model_cutoff_",cf,"_per_sample_accuracy_334_samples.tab"))

####################################################################

