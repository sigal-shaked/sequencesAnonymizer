#anonymize.type: 0 - no anonymization; 1 - model level; 2 - generation level
run.single.experiment <- function(experiment.name,d,objectid.pos,timestamp.pos,stateid.pos,secs.between.sequences,k,similarity,shingle.size=1,num.hashes=50,anonymize.type=0,min.objects=0,max.freq=1,timestamp.format='%Y-%m-%dT%H:%M:%S'){
  cur.results <- c(experiment.name,similarity,k,shingle.size,num.hashes,anonymize.type,min.objects)
  ptm <- proc.time()

  ###df<- init.transactions.dataset(d,objectid.pos=objectid.pos,timestamp.pos=timestamp.pos,stateid.pos=stateid.pos,seq_id.pos=0,secs.between.sequences=secs.between.sequences,timestamp.format)
  generator<- sequencesGenerator(d,objectid.pos=objectid.pos,timestamp.pos=timestamp.pos,stateid.pos=stateid.pos,seq_id.pos=0,secs.between.sequences=secs.between.sequences,timestamp.format)
  df <- generator$src
  init.time <- as.numeric(unlist((proc.time() - ptm)[3]))
  cur.results <- append(cur.results,init.time)

  n.objects <- length(unique(df$objectid))
  n.sequences <- length(unique(df$seq_id))
  n.records <- dim(df)[1]
  cur.results <- append(cur.results,c(n.objects,n.sequences,n.records))
  if(k==0)
  {
    cur.k=1
  } else {
    cur.k = ceiling(k/100 * n.sequences)
  }
  ptm <- proc.time()
  #####cluster data
  df.clustered <- cluster.sequences(s=df,clusters.k=cur.k,sim.measure=similarity,shingle.size=1,num.hashes=50)
  ###df.clustered <- cluster.sequences(s=generator$src,clusters.k=cur.k,sim.measure=similarity,shingle.size=1,num.hashes=50)

  cluster.time <- as.numeric(unlist((proc.time() - ptm)[3]))
  cur.results <- append(cur.results,cluster.time)

  ptm <- proc.time()
  #####model creation
  model<-build.model(clustered.data=df.clustered,min.objects=ifelse(anonymize.type==1,min.objects,0),max.freq=ifelse(anonymize.type==1,max.freq,1))
  model.size <- sum(unlist(lapply(model,nrow)))
  rm(generator)
  rm(df.clustered)
  create.model.time <- as.numeric(unlist((proc.time() - ptm)[3]))
  cur.results <- append(cur.results,create.model.time)
  cur.results <- append(cur.results,model.size)
  inverse_method=T
  ptm <- proc.time()
  #####generate data
  output.file.name <- paste0(experiment.name,'.csv')
  #registerDoParallel(detectCores())
  synthetic.data1 <- generate.synthetic.data(m=model,p.method=inverse_method,p.min.objects=ifelse(anonymize.type==2,min.objects,0),p.max.freq=ifelse(anonymize.type==2,max.freq,1))
  rm(model)
  #registerDoSEQ()
  ###write.table(synthetic.data, file = output.file.name, append = F, sep = ",",row.names = F,col.names = T)
  synthetic.obj<- sequencesGenerator(data=synthetic.data1,objectid.pos=1,timestamp.pos=3,stateid.pos=4,seq_id.pos=2,timestamp.format='%Y-%m-%d %H:%M:%S')
  synthetic.data <-synthetic.obj$src
  rm(synthetic.data1)

  data.generation.time <- as.numeric(unlist((proc.time() - ptm)[3]))
  cur.results <- append(cur.results,data.generation.time)

  ns.objects <- length(unique(synthetic.data$objectid))
  ns.sequences <- length(unique(synthetic.data$seq_id))
  ns.records <- dim(synthetic.data)[1]
  cur.results <- append(cur.results,c(ns.objects,ns.sequences,ns.records))

  #####compare datsets

  dist.lsh <- mean.seq.dist.lsh(df,synthetic.data,shingle.size=1,num.hashes=500)
  cur.results <- append(cur.results,dist.lsh)

  ###dist.lcs <- mean.seq.dist.lcs(df,synthetic.data)
  ###cur.results <- append(cur.results,dist.lcs)

  dist.obj.lsh <- mean.obj.dist.lsh(df,synthetic.data,shingle.size=1,num.hashes=500)
  cur.results <- append(cur.results,dist.obj.lsh)

  frequent.itemset <- compare.frequent.itemset(d1=df,d2=synthetic.data,top=20)
  cur.results <- append(cur.results,frequent.itemset)

  #frequent.itemset <- compare.frequent.itemset(D1=df,D2=synthetic.data,pMinSupport=0.23,top=20)
  #cur.results <- append(cur.results,as.vector(unlist(frequent.itemset)))
  #cur.results <- append(cur.results,frequent.itemset[[1]]/sum(unlist(frequent.itemset)))


  #per state and hour
  state.hour.stats <- compare.statistics.by.attribute(df,synthetic.data,relevant.col.names=c("objectid","state_id","timestamp","seq_id"),by.col.names=c("state_id","factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%H\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
  cur.results <- append(cur.results,as.vector(unlist(state.hour.stats)))

  #per state and weekday
  state.weekday.stats <- compare.statistics.by.attribute(df,synthetic.data,relevant.col.names=c("objectid","state_id","timestamp","seq_id"),by.col.names=c("state_id","factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%w\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
  cur.results <- append(cur.results,as.vector(unlist(state.weekday.stats)))

  #per hour factor
  hour.stats <- compare.statistics.by.attribute(df,synthetic.data,relevant.col.names=c("objectid","timestamp","seq_id"),by.col.names=c("factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%H\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
  cur.results <- append(cur.results,as.vector(unlist(hour.stats)))

  #per weekday factor
  weekday.stats <- compare.statistics.by.attribute(df,synthetic.data,relevant.col.names=c("objectid","timestamp","seq_id"),by.col.names=c("factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%w\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
  cur.results <- append(cur.results,as.vector(unlist(weekday.stats)))

  #per state
  state.stats <- compare.statistics.by.attribute(df,synthetic.data,relevant.col.names=c("objectid","state_id","timestamp","seq_id"),by.col.names=c("state_id"),factor.calculation=NA,seq_id.col.name="seq_id",objectid.col.name="objectid")
  cur.results <- append(cur.results,as.vector(unlist(state.stats)))

  #per date
  date.stats <- compare.statistics.by.attribute(df,synthetic.data,relevant.col.names=c("objectid","timestamp","seq_id"),by.col.names=c("factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%D\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
  cur.results <- append(cur.results,as.vector(unlist(date.stats)))


  #pdf(paste0(experiment.name,'.pdf'))
  #compare.datasets(df,synthetic.data,experiment.name)
  #dev.off()
  rm(df)
  rm(synthetic.data)
  cur.results
}
#
# ####D1 last-fm
# #####load data
# #set.seed(5)
# #registerDoParallel(detectCores())
#
# library("sequencesGenerator")
# setwd('D:\\sigal\\doctorat\\עיבוד הצעת מחקר ומצגת לבחינה')
#
# #Set parameters according to the music dataset:
# #setwd('D:\\sigal\\doctorat\\׳¢׳™׳‘׳•׳“ ׳”׳¦׳¢׳× ׳׳—׳§׳¨ ׳•׳׳¦׳’׳× ׳׳‘׳—׳™׳ ׳”')
# #input.file.name <- 'datasets\\last-fm Music Recommendation Datasets for Research\\lastfm-dataset-1K\\userid-timestamp-artid-artname-traid-traname.tsv'
# input.file.name <- 'D1.csv'
#
# #col.names<- c('userid','timestamp','musicbrainz-artist-id','artist-name','musicbrainz-track-id','track-name')
# objectid.pos<-1
# timestamp.pos<-2
# stateid.pos<-c(4,6)
# secs.between.sequences <- 140*60 #according to examined histogram
# samples <- 1
# replications<- 10
# #nrows <- 2500
# skip <- 0
# rep <- NULL
# similarity <- NULL
# k <- NULL
# exp.file.name <- 'D1_sum.csv'
# res <- c("experiment.name","similarity","k","shingle.size","num.hashes","anonymize.type","min.objects",
#          "init.time","n.objects","n.sequences","n.records","cluster.time","create.model.time","model.size","data.generation.time",
#          "dist.lsh","dist.obj.lsh","frequent.intersect",#"frequent.A_B","frequent.B_A","frequent.accuracy",
#          "state.hour.mean.sequences.diff","state.hour.mean.objects.diff","state.hour.mean.total.diff",
#          "state.hour.sd.sequences.diff","state.hour.sd.objects.diff","state.hour.sd.total.diff",
#          "state.weekday.mean.sequences.diff","state.weekday.mean.objects.diff","state.weekday.mean.total.diff",
#          "state.weekday.sd.sequences.diff","state.weekday.sd.objects.diff","state.weekday.sd.total.diff",
#          "hour.mean.sequences.diff","hour.mean.objects.diff","hour.mean.total.diff",
#          "hour.sd.sequences.diff","hour.sd.objects.diff","hour.sd.total.diff",
#          "weekday.mean.sequences.diff","weekday.mean.objects.diff","weekday.mean.total.diff",
#          "weekday.sd.sequences.diff","weekday.sd.objects.diff","weekday.sd.total.diff",
#          "state.mean.sequences.diff","state.mean.objects.diff","state.mean.total.diff",
#          "state.sd.sequences.diff","state.sd.objects.diff","state.sd.total.diff",
#          "date.mean.sequences.diff","date.mean.objects.diff","date.mean.total.diff",
#          "date.sd.sequences.diff","date.sd.objects.diff","date.sd.total.diff")
#
# write.table(t(res), file = exp.file.name, append = T, sep = ",",row.names = F,col.names = F)
# for (sample.size in c(5000)){#},15000,25000)){#c(1000,seq(5000, 25000, by = 10000))){
#   for (samp in 1: samples){
#     f <- read.delim(input.file.name, header=T, sep=',',skip=skip,nrows=sample.size)
#     skip <- skip + sample.size
#     for (similarity in c(1)){#,1)){
#       for (k in seq(0, 100, by = 20)){
#         set.seed(123)
#         all.rep.results<- matrix(data = NA, nrow = 0, ncol = length(res))
#         for (rep in 1:replications){
#           experiment.name <- paste("D1",samp,rep,sample.size,similarity,k,sep="_")
#           res <- run.single.experiment(experiment.name,d=f,objectid.pos,timestamp.pos,stateid.pos,secs.between.sequences,k,similarity,timestamp.format='%Y-%m-%dT%H:%M:%S')
#           write.table(t(res), file = exp.file.name, append = T, sep = ",",row.names = F,col.names = F)
#           all.rep.results <- rbind(all.rep.results,t(res))
#           #rm(res)
#         }
#         if (replications>1){
#           all.rep.results<- all.rep.results[,-1]
#           class(all.rep.results) <- "numeric"
#           mean.rep.results<- apply(all.rep.results,2,mean)
#           write.table(t(c(paste("mean_D1",samp,sample.size,similarity,k,sep="_"),mean.rep.results)), file = exp.file.name, append = T, sep = ",",row.names = F,col.names = F)
#           sd.rep.results<- apply(all.rep.results,2,sd)
#           write.table(t(c(paste("sd_D1",samp,sample.size,similarity,k,sep="_"),sd.rep.results)), file = exp.file.name, append = T, sep = ",",row.names = F,col.names = F)
#           rm(all.rep.results)
#         }
#       }
#     }
#     rm(f)
#   }
# }
# #registerDoSEQ()
#
#
#
#
#
#
#
#
#
#
# library("sequencesGenerator")
# setwd('D:\\sigal\\doctorat\\עיבוד הצעת מחקר ומצגת לבחינה')
#
# ####D2 d4d
# #####load data
# #set.seed(5)
# #registerDoParallel(detectCores())
# #Set parameters according to the music dataset:
# #47183543 rows
# input.file.name <- 'D2.csv'
# #colnames(f1)<- c("objectid","timestamp","state_id")
# objectid.pos<-1
# timestamp.pos<-2
# stateid.pos<-3
# secs.between.sequences <- 1750*60 #according to examined histogram based on first 100000 records
# samples <- 5
# replications <- 5
# #nrows <- 2500
# skip <- 0
# rep <- NULL
# similarity <- NULL
# k <- NULL
# exp.file.name <- 'D2_sum.csv'
# res <- c("experiment.name","similarity","k","shingle.size","num.hashes","anonymize.type","min.objects",
#          "init.time","n.objects","n.sequences","n.records","cluster.time","create.model.time","model.size","data.generation.time",
#          "dist.lsh","dist.obj.lsh","frequent.intersect",#"frequent.A_B","frequent.B_A","frequent.accuracy",
#          "state.hour.mean.sequences.diff","state.hour.mean.objects.diff","state.hour.mean.total.diff",
#          "state.hour.sd.sequences.diff","state.hour.sd.objects.diff","state.hour.sd.total.diff",
#          "state.weekday.mean.sequences.diff","state.weekday.mean.objects.diff","state.weekday.mean.total.diff",
#          "state.weekday.sd.sequences.diff","state.weekday.sd.objects.diff","state.weekday.sd.total.diff",
#          "hour.mean.sequences.diff","hour.mean.objects.diff","hour.mean.total.diff",
#          "hour.sd.sequences.diff","hour.sd.objects.diff","hour.sd.total.diff",
#          "weekday.mean.sequences.diff","weekday.mean.objects.diff","weekday.mean.total.diff",
#          "weekday.sd.sequences.diff","weekday.sd.objects.diff","weekday.sd.total.diff",
#          "state.mean.sequences.diff","state.mean.objects.diff","state.mean.total.diff",
#          "state.sd.sequences.diff","state.sd.objects.diff","state.sd.total.diff",
#          "date.mean.sequences.diff","date.mean.objects.diff","date.mean.total.diff",
#          "date.sd.sequences.diff","date.sd.objects.diff","date.sd.total.diff")
#
# write.table(t(res), file = exp.file.name, append = T, sep = ",",row.names = F,col.names = F)
# for (sample.size in 5000:5000){#seq(5000, 25000, by = 5000)){
#   for (samp in 1:1){#samples){
#     f <- read.delim(input.file.name, header=T, sep=',',skip=skip,nrows=sample.size,stringsAsFactors =F)
#     skip <- skip + sample.size
#     for (similarity in 1:1){#c(0,1)){
#       for (k in 1:1){#c(1,seq(100, 300, by = 100))){
#         set.seed(123)
#         all.rep.results<- matrix(data = NA, nrow = 0, ncol = length(res))
#         for (rep in 1:replications){#5){
#           experiment.name <- paste("D2",samp,rep,sample.size,similarity,k,sep="_")
#           res <- run.single.experiment(experiment.name,d=f,objectid.pos,timestamp.pos,stateid.pos,secs.between.sequences,k,similarity,timestamp.format='%Y-%m-%d %H:%M:%S')
#           write.table(t(res), file = exp.file.name, append = T, sep = ",",row.names = F,col.names = F)
#           all.rep.results <- rbind(all.rep.results,t(res))
#         }
#         if (replications>1){
#           all.rep.results<- all.rep.results[,-1]
#           class(all.rep.results) <- "numeric"
#           mean.rep.results<- apply(all.rep.results,2,mean)
#           write.table(t(c(paste("mean_D2",samp,sample.size,similarity,k,sep="_"),mean.rep.results)), file = exp.file.name, append = T, sep = ",",row.names = F,col.names = F)
#           sd.rep.results<- apply(all.rep.results,2,sd)
#           write.table(t(c(paste("sd_D2",samp,sample.size,similarity,k,sep="_"),sd.rep.results)), file = exp.file.name, append = T, sep = ",",row.names = F,col.names = F)
#           rm(all.rep.results)
#         }
#       }
#     }
#     rm(f)
#   }
# }
# #registerDoSEQ()
#
#
