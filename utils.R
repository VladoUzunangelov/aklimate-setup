##commonly used functions
##(c) Vlado Uzunangelov 2013-2016
##uzunangelov@soe.ucsc.edu

#converts membership matrix (non-zero entries above cutoff are members) to a list of names of members
#e.g. if you  in a gene x samples 0-1 membership matrix, it will return
#a list where each entry is a vector of the names of mutated genes for each sample
matrix.to.list.membership <- function(member.mat,cutoff=0){
  bool.mat <- member.mat>cutoff
  res <- lapply(1:ncol(bool.mat),function(x,data){
    rownames(data)[data[,x]]
  },data=bool.mat)
  names(res) <- colnames(member.mat)
  #return only non-empty vectors
  res<-res[sapply(res,length)>0]
  return(res)
}

##reverse operation of matrix.to.list.membership
list.to.matrix.membership <- function(member.list){
  list.df <- lapply(1:length(member.list),function(x) {
    res <- rep(1,length(member.list[[x]]))
    names(res) <- member.list[[x]]
    res <- data.frame(res,stringsAsFactors=FALSE)
    colnames(res) <- names(member.list)[x]
    return(res)
                  })
  names(list.df) <- names(member.list)

  if(length(list.df)>1){
      out <- merge.data(list.df,fill.by=0)
  } else {
      out <- as.matrix(list.df[[1]])
  }
  return(out)
}



#the names of the columns of all data frames must be different!! - otherwise, it gets plastered together in unexpected ways!! 
#fill.by parameter - 1 for expression data (log-transformed), 0 for genomic
#data (not log-transformed)
##col.id should only be set if "id" clashes with the name of another column!!
merge.data <- function(ldfs,fill.by=1,col.id="id"){

library(plyr)
library(data.table)

##need them to share a column name so that they can get merged on it
prepped.data <- lapply(c(1:length(ldfs)),function(x,data){
    y <- data.frame(rownames(data[[x]]),data[[x]],stringsAsFactors=FALSE)
    colnames(y) <- c(col.id,colnames(data[[x]]))
    return(y)
},data=ldfs)

out.data <- join_all(prepped.data, by=col.id,type="full")
##out.data<- Reduce(function(x,y) merge(x,y,by=c("gene.id"),all=TRUE),prepped.data)
row.names(out.data) <- out.data[,1]
out.data <- out.data[,-1]
##need to use data.table because the old method of replacing NAs with zeros
##via is.na() is painfully slow
out.data <- data.table(out.data,keep.rownames=TRUE)

for (j in colnames(out.data)){
    set(out.data,which(is.na(out.data[[j]])),j,fill.by)
}

out.data <- as.data.frame(out.data)
row.names(out.data) <- out.data[,"rn"]
out.data <- subset(out.data,select=-c(rn))

return(as.matrix(out.data))
}

##############################################
#list.A - list of sample names by tissue
#list.B - list of samples names by cluster - need to format apcluster output before using it
get.intersection.table <- function(list.A,list.B){
  
func.t<-function(i,j,list.1,list.2) {
 set1<-list.1[[i]]
 set2<-list.2[[j]]
 return (length(intersect(set1,set2)))
}

vect.t<-Vectorize(func.t,vectorize.args=list("i","j"))

res<-outer(1:length(list.A),1:length(list.B),vect.t,list.1=list.A,list.2=list.B)

rownames(res) <- names(list.A)
colnames(res) <- names(list.B)
return(res)
}

###########################################
#more general version of get.intersection.table
#comparator func can now take more than two arguments, passed to main function through the ellipsis
# the output of comparator function should be a single number,though (i.e. p-value, etc.)
get.cross.table <- function(list.A,list.B,comparator.func,parallel=FALSE,...){

  
func.t<-function(i,j,list.1,list.2,...) {
 set1<-list.1[[i]]
 set2<-list.2[[j]]
 return (comparator.func(set1,set2,...))
}

if (parallel) {
    cores <- if(getDoParRegistered()) getDoParWorkers() else detectCores()-2
    vect.t<-Vectorize.parallel(func.t,vectorize.args=list("i","j"),MC.CORES = cores,...)
} else {
    vect.t<-Vectorize(func.t,vectorize.args=list("i","j"),...)
}
res<-outer(1:length(list.A),1:length(list.B),vect.t,list.1=list.A,list.2=list.B)

rownames(res) <- names(list.A)
colnames(res) <- names(list.B)
return(res)
}

##Computes a statistic for each pair of entries in two lists
##list.1 - can be a data.frame (list of columns), will represent rows in output matrix
##list.2 - can be a data.frame , will represent columns in output matrix
##comp.func - a function whose arguments are source for the (full) list and target for the individual entry (vector) being compared against each element of the list
##output is a matrix of size length(list.1)xlength(list.2) with each column representing
##the statistics computed by comp.func(list.1,list.2[[x]])

get.cross.table.fast <- function(list.1,list.2,comp.func,parallel=FALSE) {
    if (parallel) {
        cores <- if(getDoParRegistered()) getDoParWorkers() else detectCores()-2
        out <- mcmapply(comp.func,target=list.2,MoreArgs=list(source=list.1),mc.cores=cores)
    } else {
        out <- mapply(comp.func,target=list.2,MoreArgs=list(source=list.1))
    }
    return(out)
}

#parallellized version of Vectorize 
Vectorize.parallel <-
    function (FUN, vectorize.args = arg.names,simplify = TRUE, use.names = TRUE,cores=10,preschedule=TRUE) 
{
    library(parallel)
    arg.names <- as.list(formals(FUN))
    arg.names[["..."]] <- NULL
    arg.names <- names(arg.names)
    vectorize.args <- as.character(vectorize.args)
    if (!length(vectorize.args)) 
        return(FUN)
    if (!all(vectorize.args %in% arg.names)) 
        stop("must specify formal argument names to vectorize")

    collisions <- arg.names %in% c("FUN", "SIMPLIFY", "USE.NAMES", 
                                   "MC.CORES","MC.PRESCHEDULE","vectorize.args")
    if (any(collisions)) 
        stop(sQuote("FUN"), " may not have argument(s) named ", 
             paste(sQuote(arg.names[collisions]), collapse = ", "))

    FUNV <- function() {
        args <- lapply(as.list(match.call())[-1L], eval, parent.frame())
        names <- if (is.null(names(args))) 
                     character(length(args))
                 else names(args)
        dovec <- names %in% vectorize.args

        do.call("mcmapply", c(FUN = FUN, args[dovec], MoreArgs = list(args[!dovec]),SIMPLIFY = simplify, USE.NAMES = use.names,mc.cores=cores,mc.preschedule=preschedule))

    }
    formals(FUNV) <- formals(FUN)
    FUNV
}

#summarizes the data in a matrix by some cluster list
# can be done by rows or columns - default is rows
#nomatch=0 is from the description of intersection in the match() help file

# faster than converting matrix to a list of rows/columns and
#using get.cross.table with the fmatch function - how??

matrix.by.list.summary <-
function(mat,cluster.list,summary.func,dim.ind=1,...){

  library(fastmatch)
 
  res <- apply(mat,dim.ind,function(x){
    sapply(cluster.list,function(y){
        #order or arguments to (f)match matters - it returns
        #positions of first argument entries in the second one!!
               summary.func(x[fmatch(y,names(x),nomatch = 0)],...)
           })
})
    

return(t(res))
}


######################################################


##Description: Wrapper around write.table that adds the name of the first column
##Author:Vlado Uzunangelov
##Args:
##df - the data frame to write
##row.names.id - name for the first column
##out.file - name of file to write

write.df <- function(df,row.names.id='',out.file){
  output <- cbind(rownames(df),df)
  colnames(output)[1] <-row.names.id
  write.table(output ,file=out.file,quote=FALSE,append=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
}

##Description: Wrapper around read.table that returns a list of two entries:
## df - the data frame read
## colname - the name of the column that contained the  rownames
read.df <- function(in.file){
    res <- read.table(in.file,header=TRUE,row.names=NULL,check.names=FALSE,stringsAsFactors=FALSE)
    cname <- colnames(res)[1]
    rownames(res) <- res[,1]
    res <- res[,-1]
    return(list(df=res,colname=cname))
}

#########################################

##reads file in list_t format
## Preferable if the list.t is of the form object1 score1 object2 score2
##Caveat: if the name-to-score-separator appears in a name, this function will read the input wrong!!
##For list.t consisting of names only, use readSetList!
read.list.t <- function(file, gene_delim="\t",score_delim="$") {

        result <- lapply(strsplit(readLines(file), gene_delim),
                        function(x,delim) {
                                lapply(x, function(y,split) strsplit(y,split)[[1]][1],split=delim)
                        },delim=score_delim)
        result<-lapply(result,unlist)
        names(result) <- sapply(result, "[[", 1)
        result <- lapply(result, function(x) as.vector(x[-1]))

        return(result)
}



#####################################################
##reads a file in a listt format,
##e.g. Name1 member1 member2 ...
##     Name2 member1 member2 ...
##returns a list with each line representing a separate vector object named by first entry on line and having as members subsequent entries

readSetList <- function(file, delim="\t") {
  l <- readLines(file)
  l.fields <- strsplit(l, delim)
  r <- lapply(l.fields, function(x) as.vector(x[-1]))
  names(r) <- sapply(l.fields, "[[", 1)
  return(r)
}


######################################################
#reverse operation of readSetList
#setList needs to be named appropriately
writeSetList <- function(setList,out.file,delim="\t"){
    #since we are appending to a file,
    #we need to check if it exists and remove it
    if(out.file!=stdout() && file.exists(out.file)){
        file.remove(out.file)
    }
    for (i in 1:length(setList)){
        write(paste(names(setList)[i],paste(setList[[i]],collapse=delim),sep=delim),file=out.file,append=TRUE)
    }

}


#########
#methods
#difference -takes the set from the first list and subtracts elements of the corresponding set in the second list, only keeps elements of first set
#intersection - intersects corresponding sets form first and second list, only keeps sets that appear in both lists
#union - does the union of corresponding sets from first and second list, keeps union of all sets
mergeSetList <- function(setList1,setList2,method=c("union","intersection","difference")){
        method <- match.arg(method)
        common <- intersect(names(setList1),names(setList2))
        switch(method,union={
            res <- c(setList1[setdiff(names(setList1),common)],setList2[setdiff(names(setList2),common)])
            tmp <- lapply(common,function(x) union(setList1[[x]],setList2[[x]]))
            names(tmp) <- common
            res <- c(res,tmp)
            
        },intersection= {
            res <- lapply(common,function(x) intersect(setList1[[x]],setList2[[x]]))
            names(res) <- common
        }, difference= {
            res <- setList1[setdiff(names(setList1),common)]
            tmp <- lapply(common,function(x) setdiff(setList1[[x]],setList2[[x]]))
            names(tmp) <- common
            res <- c(res,tmp)
        })
        return(res[sapply(res,length)>0])       
}
####################################
#cuts a vector or a list or a data.frame into pieces of equal size (all but last one)
#returns a list of the pieces
# The slices match the type of the parent input
#n-length of slice
slice<-function(x,n) {
    #matrices are atomic, so need an extra condition
    if((is.null(dim(x)) || is.atomic(x)) && class(x)!="matrix"){
        ##lists & vectors
        N <- length(x)
        res <- lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
        
    } else{
        #data frame or matrix
        N<-dim(x)[1]
        res <- lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N),,drop=FALSE])
    }
       return(res)
}

###################################
#parallelizes the supplied function over multiple cpus, multiple nodes, or combination when cubmitted via a job scheduler
#args- list of arguments supplied to func, first argument is the one we iterate over- it has to be one of the following two options:
#1.a vector of a single variable, or a list of vectors, where each vector is a single observation of a variable (e.g. a gene set)
#2.a data.frame of variables(e.g. exploration of parameter space) where each row is a set of variable values to be used in a given run (allows for variation of multiple parameters simultaneously)
#IT WILL NOT WORK IF ITERAND IS A LIST OF DATA.FRAMES - either submit a list of scalars/vectors, or a combined data.frame for iterand

#type - nodes is trying to parallelize across multiple workstations (i.e. pk), cpus is parallelizing over multiple cores on the same workstation (i.e. cruncher), batch parallelizes the operation when computations are done via a script submitted to a scheduling program (qsub) in an openMPI framework (you still need to specify n.processes,though)
#default type is "batch"

#individual - if iterand consists of large jobs that need to be run individually, then each process should be receiving only one element of iterand at a time (individual =TRUE)
#otherwise, send many small jobs in the same run by setting individual =FALSE)
#individual should be set to FALSE only if the iterand consists of one variable only!! (i.e. a list/vector,and not a data.frame)
parallelize <- function(func,args,hosts=paste0("kkr17u",seq(11,20,1)),n.processes=length(hosts),type=c("batch","nodes","cpus"),individual=TRUE,...)
    {
        #library(parallel)
        type <- match.arg(type)

        if(type=="nodes")
           {
               #working on multiple machines
               # can have several processes on the same machine if a host
               #name is repetated multiple times
               library(parallel)
               library(doSNOW)
               cluster <- makePSOCKcluster(hosts,outfile="",...)
               registerDoSNOW(cluster)
           }
       else if(type=="cpus")
           {
               #multiple processes spawned on one machine
               library(parallel)
               library(doParallel)
               cluster <- makeCluster(n.processes,outfile="",...)
               registerDoParallel(cores=n.processes)
               
           }
       else
           {
               #running in batch with a scheduling software
               ## currently not working - getting a weird segfault in Rmpi (2016/03/11)
               library(Rmpi)

               library(doMPI)
               cluster <- startMPIcluster()
               registerDoMPI(cluster)
               n.processes <- mpi.universe.size()-1
               
           }

        length.iterand <- if (is.null(dim(args[[1]]))) length(args[[1]]) else dim(args[[1]])[1]

        if(individual){
            #each iterand is run in a separate process (large jobs)
            iterand <- slice(args[[1]],1)
        } else {
            #several iterands runs in same process(small jobs, normally not parameter space exploration - i.e. iterand is a vector, not a list)
            iterand <- slice(args[[1]],ceiling(length.iterand/n.processes))
        }
        
        res <- foreach(it=iter(iterand)) %dopar% {
           
            #if iterand is a list/vector, slice() will return elements of
            #that that list that are not named(list) or possibly named
            #incorrectly (vector); in addition, slice() strips names of
            #list entries if iter is a list

            it <- list(it)
            names(it) <- names(args)[1]

            
               #I get an error when using args[-1] if only 1 argument
              if (length(args)>1) {
                   node.args <- c(it,args[-1])
               }
               else {
                   node.args <- it
               }

           #actual work
           #need to clean memory each time just in case
          
           do.call(func,node.args)
        }

        switch(type,
               nodes={
                   invisible(stopCluster(cluster))
               },
               cpus={
                   invisible(stopCluster(cluster))
               },
               batch={
                   closeCluster(cluster)
                   invisible(capture.output(mpi.finalize()))
               })
        
        return(res)
    }
        

###################################
#quick negative Euclidean distance
#the rows should be the features you are computing pairwise similarities for (i.e. samples)
fast.neg.euclidean <- function(mat,pow=1){
    library(fields)
    res <- rdist(mat)

    res <- -res^pow
    rownames(res) <- colnames(res) <- rownames(mat)
    return(res)
}
###########################################
#input - numerical matrix
#output will have -out.range.mult and 0 as ranges
scale.to.range <- function(mat,out.range.mult=100){
    (mat-max(mat))/(max(mat)-min(mat))*out.range.mult
}

###############################
#alternative to the one above - verdict is out which one is better
#taken from Hadley Wickham's scales package

rescale <-function (x, to = c(0, 1), from = range(x, na.rm = TRUE)) 
{
    if (zero_range(from) || zero_range(to)) 
        return(rep(mean(to), length(x)))
    (x - from[1])/diff(from) * diff(to) + to[1]
}

zero_range <-
    function (x, tol = 1000 * .Machine$double.eps) 
{
    if (length(x) == 1) 
        return(TRUE)
    if (length(x) != 2) 
        stop("x must be length 1 or 2")
    if (any(is.na(x))) 
        return(NA)
    if (x[1] == x[2]) 
        return(TRUE)
    if (all(is.infinite(x))) 
        return(FALSE)
    m <- min(abs(x))
    if (m == 0) 
        return(FALSE)
    abs((x[1] - x[2])/m) < tol
}

########################################
replace.infs <- function(mat,replacement=-100000){
    mat[is.infinite(mat)] <- replacement
    return(mat)
}


#################################

#currying a function - from package roxygen
## Curry <-
##     function (FUN, ...) 
## {
##     .orig = list(...)
##     function(...) do.call(FUN, c(.orig, list(...)))
## }


Curry <- function(FUN, ...) {
      args <- match.call(expand.dots = FALSE)$...
        args$... <- as.name("...")

        env <- new.env(parent = parent.frame())

        if (is.name(FUN)) {
                fname <- FUN
            } else if (is.character(FUN)) {
                    fname <- as.name(FUN)
                } else if (is.function(FUN)){
                        fname <- as.name("FUN")
                            env$FUN <- FUN
                    } else {
                            stop("FUN not function or name of function")
                        }
        curry_call <- as.call(c(list(fname), args))

        f <- eval(call("function", as.pairlist(alist(... = )), curry_call))
        environment(f) <- env
        f
  }


###########################################
#filter a list based on list entries length
threshold.list <- function(a.list,size=15){
    
    if(length(a.list)==0) return(list())
    
    res <- a.list[sapply(a.list,function(x) length(x)>size)]
    return(res)
}

#####################################
#filter list of sets by limiting to those included in/excluded from a predefined set
#setlist - the list of sets that need filtering
#universe - the set that we are filtering the list of sets with
#include - if TRUE (default) keep all members of the sets that are also in the universe; if FALSE, keep only those that are NOT in the universe!!

restrict.list <- function(setlist,universe,include=TRUE) {
    res <- lapply(setlist,function(x){
        if(include) {
            out <- x[x%in%universe]
        }
        else {
            out <- x[!x%in%universe]
        }
        return(out)
    })
    has.memb <- sapply(res,length)>0
    return(res[has.memb])
}
####################################

#########################
##expand a list to its component objects
##if the names of the output variables are not supplied, use names of list variables
#values is a list or a vector of variables we want to expand 
expand <- function(vals,nms=NULL){
    if (is.null(nms)){
        nms=names(vals)
    }
    mapply(assign, nms, vals, MoreArgs = list(envir = parent.frame()))
    ##just so it does not return the mapply output
    invisible()
}


##############################################
#convert a clustering solution of the form list of names (cluster indices are those of the list entry) to a vector of cluster memberships

list.to.flat <- function(cluster.list){
    res <- lapply(1:length(cluster.list),function(x,clusters) {
        cl <- rep(x,length(clusters[[x]]))
        names(cl) <- clusters[[x]]
        return(cl)
    },clusters=cluster.list)
    
    return(sort(unlist(res)))
}


#############################################################
#opposite of function above
flat.to.list <- function(cluster.vector){

    ngrps <- sort(unique(cluster.vector))
    res <- lapply(ngrps,function(x) names(cluster.vector)[cluster.vector==x])
    names(res) <- ngrps
    return(res)
}

##########################################################
##filter set list
#sets need to be ordered with most important first, then next most important, etc.
## First set is always included
## cutoff - value between 0 and 1. If the overlap with previosly included sets exceeds cutoff, the evaluated set is discarded
## overlap = size(intersection)/min(size(setA),size(setB))
##if weighted=TRUE, list consists of vectors of weights normalize to sum=1
##overlap is computed on sum of common weights
filter.set.list <- function(sets,cutoff=0.65,topn=length(sets),weighted=FALSE){

    lengths <- sapply(sets,length)
    ## filter out zero-length sets
    sets <- sets[lengths>0]
    lengths <- lengths[lengths>0]
    res <- rep(TRUE,length(sets))

    for (i in 2:length(sets)) {
        counter <- 1
        keep <- TRUE
        included <- which(res[1:(i-1)])
        if (length(included)==topn) {
            res[i:length(res)] <- FALSE
            break
        }

        if(weighted) {
            while(keep && counter <= length(included) ) {
                common <- intersect(names(sets[[i]]),names(sets[[included[counter]]]))
                keep <- (sum(sets[[i]][common])/sum(sets[[i]])) <= cutoff
                ##&& (sum(sets[[included[counter]]][common])/sum(sets[[included[counter]]])) <= cutoff
                counter <- counter + 1
            }
        } else {
            while(keep && counter <= length(included) ) {
            
                keep <- sum(sets[[i]]%in%sets[[included[counter]]]) <= min(lengths[i],lengths[included[counter]])*cutoff
                counter <- counter + 1
            }
        }
        
        res[i] <- keep
    }
        
    return(sets[res])

}
