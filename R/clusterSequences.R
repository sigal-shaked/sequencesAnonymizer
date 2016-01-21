require(cluster)
require(MinHash)
require(TraMineR)
require(gmp)

# s - structured set of transactions
# clusters.k - number of clusters
# sim.measure - 0 for LCS, 1 for LSH
# the folowing parameters are only required for LSH (sim.measure=1)
# shingle.size - size of subsequences to be hashed (1 for a, 2 for a-b, 3 for a-b-c ect.)
# num.hashes - number of hash function used for signaturing a sequence
cluster.sequences <- function(s,clusters.k,sim.measure,shingle.size,num.hashes)
{
  res <- calc.seq.obj(s)
  sequences <- res[[2]]
  seq.idx <- res[[1]]
  #calculate distance matrix
  if (sim.measure==0){
    dist.m<- calc.pairwise.lcs(sequences)
  } else if (sim.measure==1){
    dist.m<- MinHash::calc.lsh(sequences,NULL,shingle.size,num.hashes)
  }

  #cluster sequences according to distance matrix
  clusterward <- cluster::agnes(dist.m , diss = T, method = "ward")
  rm(dist.m)
  cur.k<-min(clusters.k,length(seq.idx))
  clusters <- cutree(clusterward, k = cur.k)

  #   #ploting clusters
  #    cpal(D.seq)<-c(1:227)
  #    lab <- factor(clusters, labels = paste("Cluster", 1:clusters.k))
  #
  #   seqplot(D.seq, type="f", group=lab)
  #   seqplot(D.seq, type="i", group=lab)
  #   seqrplot(D.seq, group = lab, dist.matrix = dist,criterion="density",trep=0.35,ylab=c(0,20))#,
  #   seqdplot(D.seq, group = lab)#state distribution across time from optimal matching distance ,



  #seq.rep<- seqrep(D.seq, group = lab, dist.matrix = dist,
  #         criteria = "centrality", nrep=1)

  #seqrplot(D.seq, group = lab, dist.matrix = dist,criterion="density",trep=0.35)#,
  #         criteria = "prob", nrep=1)


  # table(clusters)
  #adding cluster id to origin data
  #seqfplot(d4d.seq, group = clusters, pbarw = T)
  clustered.seqe <- cbind(seq.idx, as.vector(clusters))
  colnames(clustered.seqe) <- c( 'id', 'cluster_id')
  clustered.transactions <- merge(s,clustered.seqe,by.y="id",by.x="seq_id",all = TRUE)
  rm(clusterward)
  rm(clusters)
  rm(clustered.seqe)
  #rename "id" column to "seq_id"
  names(clustered.transactions)[1]<- "seq_id"
  clustered.transactions
}

#df.clustered <- cluster.sequences(s=df,clusters.k=3,sim.measure=1,shingle.size=1,num.hashes=50)

calc.seq.obj <- function(d)
{
  #creating TSE object and converting it into STS object
  #d.seq_id <- as.numeric(factor(d$seq_id))
  ###d.seq_id <- d$seq_id
  #d.seq_id <-as.numeric(factor(paste(d[,"objectid"],d[,"seq_id"],sep="-")))
  d.seqe <- TraMineR::seqecreate(id = as.numeric(d[,"seq_id"]), timestamp = (as.numeric(as.POSIXct(unlist(d[,"timestamp"]),origin='1970-01-01')) - as.numeric(as.POSIXct('1970-01-01',origin='1970-01-01'))),event = d[,"state_id"])
  tmp <- as.character(d.seqe)
  open <- "("
  close <- ")"
  tmp <- gsub("[1234567890e+.]+-", replacement=",", x=tmp, fixed=FALSE)
  tmp <- gsub(",", replacement="-", x=tmp, fixed=TRUE)
  tmp <- gsub(open, replacement="", x=tmp, fixed=TRUE)
  tmp <- gsub(close, replacement="", x=tmp, fixed=TRUE)
  tmp <- gsub("--", replacement="-", x=tmp, fixed=TRUE)
  tmp <- substr(tmp ,2,nchar(tmp))
  ## Make a state sequence object from the lists of events
  d.seq <- TraMineR::seqdef(tmp, sep=",")
  rm(tmp)

  #seq.idx <- group_by(d,seq_id)
  #seq.idx <-arrange(seq.idx,timestamp)
  #seq.idx <- filter(seq.idx, timestamp == first(timestamp))$seq_id
  #list(seq.idx,d.seq)
  list(TraMineR::seqeid(d.seqe),d.seq)
}

