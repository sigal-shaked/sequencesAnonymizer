# require(gmp)
# require(TraMineR)
# #library(dplyr)
#
#
#
# #notice: levels are important for factorizing similar states the same way
# convert.to.shingles <- function(s,p.levels,shingle.size)
# {
#   d.seq <- s
#   d.seq[d.seq=="%"]<-NA
#   #first loop for 1-gram, second for 2-gram ect.
#   cnt<-1
#   while (cnt<=shingle.size){
#     if(cnt==1){
#       cur.levels <- as.matrix(d.seq)
#     } else{
#       next.states<- as.matrix(t(apply(cur.levels,1,function(x){c(x[cnt:length(x)],rep(NA,cnt))})))
#       cur.levels <- matrix(paste(as.matrix(cur.levels), next.states,sep="-"),nrow=nrow(next.states))
#     }
#     cnt <- cnt+1
#   }
#   matrix(as.numeric(factor(as.matrix(cur.levels),levels=p.levels)),dim(cur.levels)[1])
# }
#
# create.hash.functions <- function(max.shingle.id,num.hashes){
#   #set.seed(135)
#   find.next.prime <- function(x)
#   {
#     x <- x+1
#     while (isprime(x)<2){
#       x <- x+1
#     }
#     x
#   }
#   # Our random hash function will take the form of:
#   #   h(x) = (a*x + b) % c
#   # Where 'x' is the input value, 'a' and 'b' are random coefficients, and 'c' is
#   # a prime number just greater than maxShingleId.
#   coef.c <- find.next.prime(max.shingle.id)
#   # Generate a list of 'num.hashes' random unique coefficients for the random hash functions
#   if (max.shingle.id==1){
#     coef.a <- rep(max.shingle.id, num.hashes)
#     coef.b <- rep(max.shingle.id, num.hashes)
#   } else if(max.shingle.id>num.hashes){
#     coef.a <- sample(1:max.shingle.id, num.hashes, replace=F)
#     coef.b <- sample(1:max.shingle.id, num.hashes, replace=F)
#   } else{
#     coef.a <- sample(1:max.shingle.id, num.hashes, replace=T)
#     coef.b <- sample(1:max.shingle.id, num.hashes, replace=T)
#   }
#   list(coef.a,coef.b,coef.c)
# }
#
#
#
# calc.minhash.signatures <- function(d,hash.functions.obj,num.hashes)
# {
#   coef.a <- hash.functions.obj[[1]]
#   coef.b <- hash.functions.obj[[2]]
#   coef.c <- hash.functions.obj[[3]]
#   d.signatures <- matrix(data=NA,nrow=dim(d)[1],ncol=num.hashes)
#   for (i in 1:dim(d)[1]){
#     cur.seq <- unique(d[i,])
#     cur.seq<- cur.seq[!is.na(cur.seq)]
#     for (f in 1:num.hashes){
#       d.signatures[i,f] <- min((coef.a[f]*cur.seq+coef.b[f])%%coef.c,na.rm=T)
#     }
#   }
#   d.signatures
# }
#
#
#
# calc.lsh <- function(s1,s2=NULL,shingle.size,num.hashes)
# {
#   d1.seq<-s1
#   d2.seq<-s2
#   is.pairwise=F
#   if (is.null(d2.seq)){
#     d2.seq<-d1.seq
#     is.pairwise=T
#   }
#   d1.seq[d1.seq=="%"]<-NA
#   d2.seq[d2.seq=="%"]<-NA
#   levels <- NULL #unique(c(as.vector(t(d1.seq)),as.vector(t(d2.seq))))
#   #first loop for 1-gram, second for 2-gram ect.
#   cnt<-1
#   while (cnt<=shingle.size){
#     if(cnt==1){
#       levels1 <- as.matrix(d1.seq)
#       levels2 <- as.matrix(d2.seq)
#       cur.exclude="NA"
#       exclude=as.matrix(cur.exclude)
#     } else{
#       next1<- as.matrix(t(apply(levels1,1,function(x){c(x[cnt:length(x)],rep(NA,cnt))})))
#       next2<- as.matrix(t(apply(levels2,1,function(x){c(x[cnt:length(x)],rep(NA,cnt))})))
#       levels1 <- matrix(paste(levels1, next1,sep="-"),nrow=nrow(next1))
#       levels2 <- matrix(paste(levels2, next2,sep="-"),nrow=nrow(next2))
#       cur.exclude <- paste(cur.exclude,"NA",sep="-")
#       exclude<- rbind(exclude,cur.exclude)
#     }
#     cnt <- cnt+1
#     exclude<- rbind(exclude," <NA>")
#   }
#   col.size <- max(dim(levels1)[2],dim(levels2)[2])
#   row.size <- dim(levels1)[1] +dim(levels2)[1]
#   levels<- matrix(data=NA,row.size,col.size)
#   levels[1:dim(levels1)[1],1:dim(levels1)[2]]<- levels1[,]
#   levels[(dim(levels1)[1]+1):(dim(levels1)[1]+dim(levels2)[1]),1:dim(levels2)[2]]<- levels2[,]
#   levels <- unique(as.vector(levels))
#   levels <- levels[!is.na(levels)]
#   levels <- levels[!levels %in% as.vector(exclude)]
#   levels <- factor(levels)
#   ids <- as.numeric(levels)
#
#   max.shingle.id <- max(ids[!is.na(ids)])
#   hash.functions.obj <- create.hash.functions(max.shingle.id,num.hashes)
#   d1.shingles <- convert.to.shingles(s=d1.seq,p.levels=levels,shingle.size)
#   d2.shingles <- convert.to.shingles(s=d2.seq,p.levels=levels,shingle.size)
#   d1.signatures <- calc.minhash.signatures(d1.shingles,hash.functions.obj,num.hashes)
#   d2.signatures <- calc.minhash.signatures(d2.shingles,hash.functions.obj,num.hashes)
#
#   dist <- matrix(data=NA,dim(d1.signatures)[1],dim(d2.signatures)[1])
#   max.length <- max(dim(d1.signatures)[2],dim(d2.signatures)[2])
#
#   for (i in 1:dim(d1.signatures)[1]){
#     if (is.pairwise){
#       end.point = i-1 #for pairwise comparison
#     }else{
#       end.point = dim(d2.signatures)[1] #otherwise
#     }
#     if (end.point >0){
#       for (j in 1:end.point){
#         dist[i,j] <- 1-sum(d1.signatures[i,]==d2.signatures[j,])/max.length
#       }
#     }
#   }
#   if (is.pairwise){
#     dist <- as.vector(dist[lower.tri(dist, diag = FALSE)])  ## "melt" it into a vector
#     dist
#   } else {
#     dist
#   }
# }
