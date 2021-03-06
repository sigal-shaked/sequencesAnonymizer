#anonymize.type: 0 - no anonymization; 1 - model level; 2 - generation level
run.single.experiment <- function(experiment.name,d,objectid.pos,timestamp.pos,stateid.pos,secs.between.sequences,k,similarity,p_eps,shingle.size=1,num.hashes=150,timestamp.format='%Y-%m-%dT%H:%M:%S'){
  cur.results <- c(experiment.name,similarity,k,shingle.size,num.hashes,p_eps)
  ptm <- proc.time()

  ###df<- init.transactions.dataset(d,objectid.pos=objectid.pos,timestamp.pos=timestamp.pos,stateid.pos=stateid.pos,seq_id.pos=0,secs.between.sequences=secs.between.sequences,timestamp.format)
  generator<- sequencesGenerator(d,objectid.pos=objectid.pos,timestamp.pos=timestamp.pos,stateid.pos=stateid.pos,seq_id.pos=0,secs.between.sequences=secs.between.sequences,timestamp.format)
  df <- generator$src
  init.time <- as.numeric(unlist((proc.time() - ptm)[3]))
  cur.results <- append(cur.results,init.time)

  n.objects <- length(unique(df[,"objectid"]))
  n.sequences <- length(unique(df[,"seq_id"]))
  n.records <- dim(df)[1]

  g<- dplyr::group_by(df,seq_id)
  g<- dplyr::summarize(g,size=n())
  mean.sequence.size <- mean(g$size)
  min.sequence.size <- min(g$size)
  max.sequence.size <- max(g$size)

  g<- dplyr::group_by(df,objectid)
  g<- dplyr::summarize(g,size=n())
  mean.object.size <- mean(g$size)
  min.object.size <- min(g$size)
  max.object.size <- max(g$size)

  g<- unique(df[,c("objectid","seq_id")])
  g<- dplyr::group_by(g,objectid)
  g<- dplyr::summarize(g,size=n())
  mean.seq.per.object <- mean(g$size)
  min.seq.per.object <- min(g$size)
  max.seq.per.object <- max(g$size)

  cur.results <- append(cur.results,c(n.objects,n.sequences,n.records,
                                      mean.sequence.size,min.sequence.size,max.sequence.size,
                                      mean.object.size,min.object.size,max.object.size,
                                      mean.seq.per.object,min.seq.per.object,max.seq.per.object))
  if(k==0)
  {
    cur.k=1
  } else {
    cur.k = ceiling(k/100 * n.sequences)
  }
  ptm <- proc.time()
  #####cluster data
  if(cur.k>1){
    df.clustered <- cluster.sequences(s=df,clusters.k=cur.k,sim.measure=similarity,shingle.size=1,num.hashes=150)
  } else{
    df.clustered <- df
    df.clustered$"cluster_id" <- 1
  }
  ###df.clustered <- cluster.sequences(s=generator$src,clusters.k=cur.k,sim.measure=similarity,shingle.size=1,num.hashes=50)

  cluster.time <- as.numeric(unlist((proc.time() - ptm)[3]))
  cur.results <- append(cur.results,cluster.time)
  cur.results <- append(cur.results,cur.k)

  min.cluster.size <- dplyr::group_by(df.clustered,cluster_id)
  min.objs.per.cluster<- min(dplyr::summarize(min.cluster.size,length(unique(objectid)))[2])
  min.seqs.per.cluster<- min(dplyr::summarize(min.cluster.size,length(unique(seq_id)))[2])
  mean.objs.per.cluster<- mean(unlist(dplyr::summarize(min.cluster.size,length(unique(objectid)))[2]))
  mean.seqs.per.cluster<- mean(unlist(dplyr::summarize(min.cluster.size,length(unique(seq_id)))[2]))
  cur.results <- append(cur.results,min.objs.per.cluster)
  cur.results <- append(cur.results,min.seqs.per.cluster)
  cur.results <- append(cur.results,mean.objs.per.cluster)
  cur.results <- append(cur.results,mean.seqs.per.cluster)

  ptm <- proc.time()
  #####model creation
  model<-build.model(clustered.data=df.clustered,c_eps=p_eps)
  model.size <- sum(unlist(lapply(model,nrow)))
  states<- length(unique(model$nocluster.nofactor.starting.state[,"state_id"]))
  supressed_pct<- sum(as.numeric(model$supression_log[,2]))/(model.size+sum(as.numeric(model$supression_log[,2]))-nrow(model$supression_log))
  min_support<-min(unlist(lapply(model,FUN=function(x){if( "total_numeric" %in% colnames(x)&& dim(x)[1]>0 ){min(x[,"total_numeric"],na.rm = T)}})))
  mean_support<-mean(unlist(lapply(model,FUN=function(x){if( "total_numeric" %in% colnames(x)&& dim(x)[1]>0 ){mean(x[,"total_numeric"],na.rm = T)}})))
  stdev_support<-sd(unlist(lapply(model,FUN=function(x){if( "total_numeric" %in% colnames(x) && dim(x)[1]>0){mean(x[,"total_numeric"],na.rm = T)}})))
  max_support<-max(unlist(lapply(model,FUN=function(x){if( "total_numeric" %in% colnames(x) && dim(x)[1]>0){max(x[,"total_numeric"],na.rm = T)}})))
  rm(generator)
  rm(df.clustered)
  create.model.time <- as.numeric(unlist((proc.time() - ptm)[3]))
  cur.results <- append(cur.results,c(create.model.time,model.size))
  cur.results <- append(cur.results,c(states,supressed_pct,min_support,mean_support,stdev_support,max_support))

  inverse_method=T
  ptm <- proc.time()
  #####generate data
  output.file.name <- paste0(experiment.name,'.csv')
  #registerDoParallel(detectCores())
  synthetic.data <- generate.synthetic.data(m=model,p.method=inverse_method)
  rm(model)
  #registerDoSEQ()
  ###write.table(synthetic.data, file = output.file.name, append = F, sep = ",",row.names = F,col.names = T)
  #synthetic.obj<- sequencesGenerator(data=synthetic.data1,objectid.pos=1,timestamp.pos=3,stateid.pos=4,seq_id.pos=2,timestamp.format='%Y-%m-%d %H:%M:%S')
  #synthetic.data <-synthetic.obj$src
  #rm(synthetic.data1)

  data.generation.time <- as.numeric(unlist((proc.time() - ptm)[3]))
  cur.results <- append(cur.results,data.generation.time)

  ns.objects <- length(unique(synthetic.data[,"objectid"]))
  ns.sequences <- length(unique(synthetic.data[,"seq_id"]))
  ns.records <- dim(synthetic.data)[1]
  cur.results <- append(cur.results,c(ns.objects,ns.sequences,ns.records))

  #####compare datsets
  if(ns.records>0){
    synthetic.data <- as.data.frame(synthetic.data,stringsAsFactors = F)
    synthetic.data$state_id  <- as.numeric(synthetic.data$state_id)
    dist.lsh <- mean.seq.dist.lsh(df,synthetic.data,shingle.size=1,num.hashes=200)
    cur.results <- append(cur.results,dist.lsh)

    ###dist.lcs <- mean.seq.dist.lcs(df,synthetic.data)
    ###cur.results <- append(cur.results,dist.lcs)

    dist.obj.lsh <- mean.obj.dist.lsh(df,synthetic.data,shingle.size=1,num.hashes=200)
    cur.results <- append(cur.results,dist.obj.lsh)

    frequent.itemset <- compare.frequent.itemset(d1=df,d2=synthetic.data,top=20)
    cur.results <- append(cur.results,frequent.itemset)

    #frequent.itemset <- compare.frequent.itemset(D1=df,D2=synthetic.data,pMinSupport=0.23,top=20)
    #cur.results <- append(cur.results,as.vector(unlist(frequent.itemset)))
    #cur.results <- append(cur.results,frequent.itemset[[1]]/sum(unlist(frequent.itemset)))

    #per state and hour
    state.hour.stats <- compare.statistics.by.attribute(t1=df,t2=synthetic.data,relevant.col.names=c("objectid","state_id","timestamp","seq_id"),by.col.names=c("state_id","factor"),factor.calculation="format(as.POSIXlt(timestamp), \"%H\")",seq_id.col.name="seq_id",objectid.col.name="objectid")
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

  } else{
    cur.results <- append(cur.results,rep(NA,39))
  }


  #pdf(paste0(experiment.name,'.pdf'))
  #compare.datasets(df,synthetic.data,experiment.name)
  #dev.off()
  rm(df)
  rm(synthetic.data)
  cur.results
}
