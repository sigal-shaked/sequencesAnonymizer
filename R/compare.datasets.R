#synthetic.data <- read.table("try1.csv", header=T,sep=",")
#input.file.name <- 'datasets\\last-fm Music Recommendation Datasets for Research\\lastfm-dataset-1K\\userid-timestamp-artid-artname-traid-traname.tsv'
#df<- load.transactions.dataset(input.file.name=input.file.name, header=FALSE, sep='\t',objectid.pos=objectid.pos,timestamp.pos=timestamp.pos,stateid.pos=stateid.pos,seq_id.pos=0,secs.between.sequences=secs.between.sequences,nrows=1000)

library(dplyr)
library(TraMineR)
library(MinHash)

compare.statistics.by.attribute<- function(t1,t2,relevant.col.names,by.col.names,factor.calculation=NA,seq_id.col.name,objectid.col.name){
  c.relevant.col.names <- relevant.col.names
  d1 <- t1
  d2 <- t2
  d1 <- select_(d1,.dots=relevant.col.names)
  d2 <- select_(d2,.dots=relevant.col.names)
  if(!is.na(factor.calculation)){
    d1<- mutate_(d1,factor = factor.calculation)
    d2<- mutate_(d2,factor = factor.calculation)
    c.relevant.col.names<- c(c.relevant.col.names,"factor")
  }
  summarise.seq.sentence<-paste("n_distinct(",seq_id.col.name,")")
  summarise.obj.sentence<-paste("n_distinct(",objectid.col.name,")")
  src.grouped <- group_by_(d1, .dots=by.col.names )
  src.stats<- summarise_(src.grouped, sequences=summarise.seq.sentence)
  src.stats$objects<- summarise_(src.grouped, objects=summarise.obj.sentence)$objects
  src.stats$total<- summarise(src.grouped,total=n())$total
  dest.grouped <- group_by_(d2, .dots=by.col.names )
  dest.stats<- summarise_(dest.grouped, sequences=n_distinct(seq_id.col.name),objects=n_distinct(objectid.col.name))
  dest.stats$total<- summarise(dest.grouped,total=n())$total
  compare.sets<- full_join(src.stats, dest.stats, by =by.col.names )
  compare.sets[is.na(compare.sets)] <- 0
  compare.sets<-mutate(compare.sets,sequences.diff=abs(sequences.x-sequences.y),objects.diff=abs(objects.x-objects.y),total.diff=abs(total.x-total.y))
  by.col.names_<- c(by.col.names,"sequences.diff","objects.diff","total.diff")
  compare.sets<-select_(compare.sets,.dots=by.col.names_)
  list(mean=apply(compare.sets[,setdiff(names(compare.sets),by.col.names)],2,mean),
       sd=apply(compare.sets[,setdiff(names(compare.sets),by.col.names)],2,sd))
}
# #per state and hour
# compare.statistics.by.attribute(t1=df,t2=synthetic.data,relevant.col.names=c("objectid","state_id","timestamp","seq_id"),by.col.names=c("state_id","factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%H\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
# #per state and weekday
# compare.statistics.by.attribute(t1=df,t2=synthetic.data,relevant.col.names=c("objectid","state_id","timestamp","seq_id"),by.col.names=c("state_id","factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%w\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
# #per hour factor
# compare.statistics.by.attribute(t1=df,t2=synthetic.data,relevant.col.names=c("objectid","timestamp","seq_id"),by.col.names=c("factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%H\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
# #per weekday factor
# compare.statistics.by.attribute(t1=df,t2=synthetic.data,relevant.col.names=c("objectid","timestamp","seq_id"),by.col.names=c("factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%w\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
# #per state
# compare.statistics.by.attribute(t1=df,t2=synthetic.data,relevant.col.names=c("objectid","state_id","timestamp","seq_id"),by.col.names=c("state_id"),factor.calculation=NA,seq_id.col.name="seq_id",objectid.col.name="objectid")
# #per date
# compare.statistics.by.attribute(t1=df,t2=synthetic.data,relevant.col.names=c("objectid","timestamp","seq_id"),by.col.names=c("factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%D\")",seq_id.col.name="seq_id",objectid.col.name="objectid")


mean.seq.dist.lsh <- function(d1,d2,shingle.size,num.hashes){
  op <- options(warn = (-1)) # suppress warnings
  d1.seq <- calc.seq.obj(d1)[[2]]
  d2_<- dplyr::arrange(as.data.frame(d2),seq_id,timestamp)
  d2.seq <- calc.seq.obj(d2_)[[2]]
  dist<- calc.lsh(d1.seq,d2.seq,shingle.size,num.hashes)
  op <- options(warn = (1)) # cancel suppression
  mean(apply(dist,2,min))
}
mean.obj.dist.lsh <- function(d1,d2,shingle.size,num.hashes){
  d1_ <- d1
  d1_$seq_id <- d1_$objectid
  d2_ <- d2
  d2_$seq_id <- d2_$objectid
  op <- options(warn = (-1)) # suppress warnings
  d1.seq <- calc.seq.obj(d1_)[[2]]
  d2.seq <- calc.seq.obj(d2_)[[2]]
  dist<- calc.lsh(d1.seq,d2.seq,shingle.size,num.hashes)
  op <- options(warn = (1)) # cancel suppression
  mean(apply(dist,2,min))
}
#mean.obj.dist.lsh(df,synthetic.data,shingle.size=1,num.hashes=500)
#0.1226897
#mean.seq.dist.LSH(df,synthetic.data,shingle.size=2,num.hashes=500)
#[1] 0.193931

mean.obj.dist.lcs <- function(d1,d2){
  relevant.col.names=c("objectid","seq_id","timestamp","state_id")
  seqid.col.name="seq_id"
  d.all<- rbind(d1[,relevant.col.names],d2[,relevant.col.names])
  d.all$seq_id <- d.all$objectid
  end.p1<- length(unique(d1[,seqid.col.name]))
  end.p2<- end.p1 + length(unique(d2[,seqid.col.name]))
  d.all.seq <- calc.seq.obj(d.all)[[2]]
  rm(d.all)
  dist<- calc.pairwise.lcs(d.all.seq)[1:end.p1,(end.p1+1):end.p2]
  rm(d.all.seq)
  mean(apply(dist,2,min))
}
mean.seq.dist.lcs <- function(d1,d2){
  relevant.col.names=c("objectid","seq_id","timestamp","state_id")
  seqid.col.name="seq_id"
  d.all<- rbind(d1[,relevant.col.names],d2[,relevant.col.names])
  #d.all$seq_id <- d.all$objectid
  end.p1<- length(unique(d1[,seqid.col.name]))
  end.p2<- end.p1 + length(unique(d2[,seqid.col.name]))
  d.all.seq <- calc.seq.obj(d.all)[[2]]
  rm(d.all)
  dist<- calc.pairwise.lcs(d.all.seq)[1:end.p1,(end.p1+1):end.p2]
  rm(d.all.seq)
  mean(apply(dist,2,min))
}
#
#mean.obj.dist.lcs(df,synthetic.data)
#[1] 0.1403342

# compare.frequent.itemset<- function(d1,d2,pMinSupport,seqid.col.name="seq_id",timestamp.col.name="timestamp",state.col.name="state_id",top=10){
#   d1.seqe <- seqecreate(id = as.numeric(factor(d1[,seqid.col.name])), timestamp = (as.numeric(as.POSIXct(d1[,timestamp.col.name])) - as.numeric(as.POSIXct('1970-01-01'))),event = d1[,state.col.name])
#   d2.seqe <- seqecreate(id = as.numeric(factor(d2[,seqid.col.name])), timestamp = (as.numeric(as.POSIXct(d2[,timestamp.col.name])) - as.numeric(as.POSIXct('1970-01-01'))),event = d2[,state.col.name])
#   fsubseq.1 <- seqefsub(d1.seqe, pMinSupport = pMinSupport)
#   fsubseq.2 <- seqefsub(d2.seqe, pMinSupport = pMinSupport)
#   subseq.1 <-unlist(fsubseq.1[[3]])
#   subseq.1 <-rapply(subseq.1, paste)
#   subseq.2 <-unlist(fsubseq.2[[3]])
#   subseq.2 <-rapply(subseq.2, paste)
#   top.1<- min(top,length(subseq.1))
#   top.2<- min(top,length(subseq.2))
#   rm(d1.seqe)
#   rm(d2.seqe)
#   list(intersect=length(intersect(subseq.1[1:top.1],subseq.2[1:top.2])),
#     diffA_B=length(setdiff(subseq.1[1:top.1],subseq.2[1:top.2])),
#     diffB_A=length(setdiff(subseq.2[1:top.2],subseq.1[1:top.1])))
# }
# #compare.frequent.itemset(df,synthetic.data,pMinSupport=0.25,top=20)


#compares top frequent 2-grams
compare.frequent.itemset<- function(d1,d2,top){
  tran1 <- d1
  tran1$next_seqid <- c(tran1$seq_id[2:dim(tran1)[1]],NA)
  tran1$next_state_id <- c(tran1$state_id[2:dim(tran1)[1]],NA)
  tran1$next_state_id <- ifelse(tran1$next_seqid!=tran1$seq_id, "*", tran1$next_state_id)
  tran1$next_state_id <- ifelse(is.na(tran1$next_state_id), "*", tran1$next_state_id)
  tran1 <- dplyr::group_by(tran1,state_id,next_state_id)
  tran1 <- dplyr::filter(tran1,next_state_id!="*")
  tran1 <- dplyr::summarize(tran1,cnt=n())

  top.tran1 <- dplyr::top_n(dplyr::ungroup(tran1), top, cnt)
  tran2<-as.data.frame(d2,stringsAsFactors = F)
  tran2$next_seqid <- c(tran2[2:dim(tran2)[1],"seq_id"],NA)
  tran2$next_state_id <- c((tran2$state_id[2:dim(tran2)[1]]),NA)
  tran2$next_state_id <- ifelse(tran2$next_seqid!=tran2$seq_id, "*", tran2$next_state_id)
  tran2$next_state_id <- ifelse(is.na(tran2$next_state_id), "*", tran2$next_state_id)
  tran2 <- dplyr::group_by(tran2,state_id,next_state_id)
  tran2 <- dplyr::filter(tran2,next_state_id!="*")
  tran2 <- dplyr::summarize(tran2,cnt=n())
  top.tran2 <- dplyr::top_n(dplyr::ungroup(tran2), top, cnt)
  top.tran1<-sapply(top.tran1, as.numeric)
  top.tran1<- data.frame(top.tran1, stringsAsFactors = F)
  top.tran2<-sapply(top.tran2, as.numeric)
  top.tran2<- data.frame(top.tran2, stringsAsFactors = F)
  dim(dplyr::inner_join(top.tran1,top.tran2,by=c("state_id","next_state_id")))[1]/max(dim(top.tran1)[1],dim(top.tran2)[1])
}


compare.datasets <- function(data1,data2,experiment.name){
  d1<- data1
  d2 <- data2
  max.seqid <- max(as.numeric(factor(d1$seq_id)))+ 1
  d2$seq_id <- max.seqid + as.numeric(factor(d2$seq_id))
  alphabet<- unique(c(d1$state_id,d2$state_id))
  combined <- rbind (d1[,c("objectid","seq_id","timestamp","state_id")],
                     d2[,c("objectid","seq_id","timestamp","state_id")])

  d.seq <- calc.seq.obj(combined)[[2]]
  alphabet(d.seq)<- alphabet
  cpal(d.seq)<-  alphabet
  pdf(paste0(experiment.name,'.pdf'))
  seqIplot(d.seq, group=c(rep("origin",length(unique(d1$seq_id))),rep("synthetic",length(unique(d2$seq_id)))),withlegend=F)
  dev.off()
  rm(d.seq)
  rm(d1)
  rm(d2)
  #seqIplot(d.seq, border = NA, withlegend = "right")
}
#compare.datasets(df,synthetic.data)
