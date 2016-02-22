require(dplyr)
require(smoothmest)

differential.privacy<- function(data,epsilon,by.col.names,within.col.names,statistics.calculations,statistics.col.names){
  origin_size <- dim(data)[1]
  new.data<- dplyr::mutate(data,differ=(total.x/total.y)/((total.x-1)/(total.y-1)))
  if(!is.null(new.data)){
    if(length(new.data[is.na(new.data$"differ"),"differ"])[1]>0){
      new.data[is.na(new.data$"differ"),"differ"]<- 1/0
    }
  }
  removed<- dplyr::filter(new.data , differ>exp(epsilon))# | 1/differ>exp(epsilon))
  new.data<- dplyr::filter(new.data , differ<=exp(epsilon))# , 1/differ<=exp(epsilon))
  if (!is.null(epsilon)){
    while (nrow(removed)>0 & nrow(new.data)>0){
      removed.sum <- dplyr::group_by_(removed,.dots=c(by.col.names))
      removed.sum <- dplyr::summarise(removed.sum,total.y=sum(total.x))
      if(!is.null(by.col.names)){
        new.data$idx <- do.call("paste", c(new.data[,by.col.names],sep="_"))
        removed.sum$idx <- do.call("paste", c(removed.sum[,by.col.names],sep="_"))
        removed.list <- match(new.data$idx ,removed.sum$idx, nomatch =NA)
      } else{
        removed.list <- 1
      }
      if(is.null(removed.list)){
        new.data$removed<- 0
      } else{
        removed.sum<- removed.sum$"total.y"[removed.list]
        removed.sum[is.na(removed.sum)]<- 0
        new.data$removed <- removed.sum
      }
      new.data <- dplyr::mutate(new.data,total.y=total.y-removed)
      new.data<- dplyr::mutate(new.data,differ=(total.x/total.y)/((total.x-1)/(total.y-1)))
      removed<- dplyr::filter(new.data , differ>exp(epsilon))
      new.data<- dplyr::filter(new.data , differ<=exp(epsilon))
    }
  }
  names(new.data)[names(new.data) %in% paste(statistics.col.names,".x",sep="")]<- statistics.col.names
  new.data<- dplyr::mutate(new.data,total_numeric=total)
  new.data<- dplyr::mutate(new.data,total=(total/total.y))
  new.data <- dplyr::select_(new.data ,.dots=c(by.col.names,within.col.names,statistics.col.names,"total_numeric"))
  new.data[is.na(new.data)] <- 0
  list(stats=new.data,supressed = (origin_size-dim(new.data)[1]))#/origin_size)
}




calc.single.model<- function(clustered.data,relevant.col.names,factor.calculation=NA,by.col.names,within.col.names,statistics.calculations,statistics.col.names,eps = NULL){
  d <- clustered.data
  if(dim(d)[1]==0){
    if(length(within.col.names)>0){
      col.names<-c(by.col.names,within.col.names,statistics.col.names,"total_numeric")
    } else{
      col.names<-c(by.col.names,statistics.col.names,"total_numeric")
    }
    d <- matrix(data = NA, nrow = 0, ncol = length(col.names))
    colnames(d)<- col.names
    list(stats=NULL,supressed = 0)
  } else{
    d <- dplyr::select_(d,.dots=relevant.col.names)
    if(!is.na(factor.calculation)){
      d<- dplyr::mutate_(d,factor = factor.calculation)
    }

    d.level1 <- dplyr::group_by_(d,.dots=c(by.col.names,within.col.names))
    d.level1 <- dplyr::summarise_(d.level1,.dots= statistics.calculations)
    names(d.level1)<- c(by.col.names,within.col.names,statistics.col.names)

    if(length(within.col.names)>0){
      d.level0 <- dplyr::group_by_(d,.dots=by.col.names)
      d.level0 <- dplyr::summarise_(d.level0,.dots= statistics.calculations)
      names(d.level0)<- c(by.col.names,statistics.col.names)
    }else{
      d.level0 <- dplyr::group_by_(d,.dots=c())
      d.level0 <- dplyr::summarise_(d.level0,.dots= statistics.calculations)
      names(d.level0)<- c(statistics.col.names)
    }
    #d.level0 <- dplyr::summarise_(d.level0,.dots= statistics.calculations)
    #names(d.level0)<- c(by.col.names,statistics.col.names)
    if(is.null(by.col.names) | is.null(within.col.names)){
      d <- merge(d.level1, d.level0,by=NULL)
    } else{
      d <- dplyr::inner_join(d.level1, d.level0, by = by.col.names)
    }
    d[is.na(d)] <- 0
    if(!is.null(within.col.names)){
      differential.privacy(data=d,epsilon=eps,by.col.names,within.col.names,statistics.calculations,statistics.col.names)
    } else{
      names(d)[names(d) %in% paste(statistics.col.names,".x",sep="")]<- statistics.col.names
      d <- dplyr::select_(d ,.dots=c(by.col.names,within.col.names,statistics.col.names))
      d[is.na(d)] <- 0
      list(stats=d,supressed = 0)
    }
  }

}



to.valid.matrix<- function(data){
  tmp<-data
  if(!is.null(data)){
    if(is.vector(data)){
      tmp<-matrix(unlist(data),nrow=1)
      colnames(tmp)<- names(data)
    } else{
      # if(dim(data)[1]>0){
      tmp<- as.matrix(data)
      #  }
    }
  }
  tmp
}


#add sequence depandant factor (sequence duration)
calc.seq.duration<- function(data){
  d <- data
  d <- dplyr::group_by(d,seq_id)
  d <-dplyr::arrange(d,seq_id,timestamp)
  d <- dplyr::mutate(d, seq_duration_numeric = as.numeric(last(timestamp)-first(timestamp), units = "secs"))
  d<- dplyr::ungroup(d)
  duration.intervals <- dplyr::summarise(d,r0=min(seq_duration_numeric, na.rm = T) ,r1=quantile(seq_duration_numeric, p=0.25, na.rm = T),r2=quantile(seq_duration_numeric, p=0.5, na.rm = T),r3=quantile(seq_duration_numeric, p=0.75, na.rm = T),r4=max(seq_duration_numeric, na.rm = T))
  d <- dplyr::mutate(d,seq_duration =
                ifelse(seq_duration_numeric < duration.intervals$r1, 1,
                       ifelse (seq_duration_numeric>=duration.intervals$r1 & seq_duration_numeric < duration.intervals$r2,2,
                               ifelse(seq_duration_numeric >= duration.intervals$r2 & seq_duration_numeric < duration.intervals$r3,3,
                                      ifelse(seq_duration_numeric >= duration.intervals$r3,4,NA)))))
  d
}


build.model <- function(clustered.data,c_eps)
{
  s <- clustered.data
  s$next_seqid <- c(s$seq_id[2:dim(s)[1]],NA)
  s$next_state_id <- c(s$state_id[2:dim(s)[1]],NA)
  s$next_state_id <- ifelse(s$next_seqid!=s$seq_id, NA, s$next_state_id)
  s$next_state_id <- ifelse(is.na(s$next_state_id), NA, s$next_state_id)
  s$tbe <- c(s$timestamp[2:dim(s)[1]],NA)- s$timestamp
  s$tbe <- ifelse(s$next_seqid!=s$seq_id, NA, s$tbe)
  #Calculates number of states for each object
  sequences <- dplyr::group_by(s,objectid)
  sequences <- dplyr::mutate(sequences, obj_length = n())
  sequences<- dplyr::ungroup(sequences)
  #Calculates time duration of each sequence
  sequences <- calc.seq.duration(sequences)
  #add timestamp driven factors (hour & weekday)
  sequences<-dplyr::mutate(sequences,hour=as.numeric(substr(timestamp,12,13)))#,weekday =as.POSIXlt(date)$wday)


  start.of.obj <- dplyr::group_by(sequences,objectid)
  start.of.obj <- dplyr::mutate(start.of.obj,num.seq=n_distinct(seq_id))
  start.of.obj <- dplyr::filter(start.of.obj, timestamp == first(timestamp),row_number()==1)

    #calculating starting sequences per hour in each cluster
    #   new.objects.per.date<-calc.single.model(clustered.data=start.of.obj,relevant.col.names=c("objectid","num.seq","cluster_id","timestamp","seq_duration"),factor.calculation="as.Date(as.POSIXlt(timestamp))",by.col.names = c("cluster_id","factor","seq_duration"),within.col.names=c(),statistics.calculations =c("n_distinct(objectid)","mean(num.seq)","sd(num.seq)"),statistics.col.names=c("objects","mean.num.seq","sd.num.seq"))
    #   names(new.objects.per.date)[names(new.objects.per.date)=="factor"]<-"date"
 # start.of.obj
  start.list <- dplyr::group_by_(start.of.obj,.dots=c("as.Date(as.POSIXlt(timestamp))"))
  start.list <- dplyr::summarise_(start.list,.dots= c("n_distinct(objectid)","mean(num.seq)","sd(num.seq)"))
  start.list[is.na(start.list)] <- 0
  names(start.list)<- c("date","objects","mean.num.seq","sd.num.seq")
  #supression_log <- matrix(data=NA,nrow=0,ncol=2)
  #colnames(supression_log) = c("model","supressed_amount")

  #new.objects.per.date<-calc.single.model(clustered.data=start.of.obj,relevant.col.names=c("objectid","num.seq","cluster_id","timestamp","seq_duration"),factor.calculation="as.Date(as.POSIXlt(timestamp))",by.col.names = c("cluster_id","factor","seq_duration"),within.col.names=c(),statistics.calculations =c("n_distinct(objectid)","mean(num.seq)","sd(num.seq)"),statistics.col.names=c("objects","mean.num.seq","sd.num.seq"))
  #names(new.objects.per.date)[names(new.objects.per.date)=="factor"]<-"date"

  start.of.seq <- dplyr::group_by(sequences,seq_id)
  ######
  start.of.seq <- dplyr::filter(start.of.seq, timestamp == first(timestamp),row_number()==1)

  #factor.week.start.cluster
  #cur.model<-calc.single.model(clustered.data=start.of.obj,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp"),factor.calculation="as.integer(ceiling(as.POSIXlt(timestamp)$yday/4))",by.col.names = c("cluster_id","factor"),within.col.names=c(),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  cur.model<-calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%w\"))",by.col.names = c("factor"),within.col.names=c("cluster_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  factor.1.start.cluster <- to.valid.matrix(cur.model$stats)
  supression_log <- c(model="factor.1.start.cluster",supressed_amount=cur.model$supressed)
  #factor.month.start.cluster
  cur.model<-calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp"),factor.calculation="as.POSIXlt(timestamp)$mon",by.col.names = c("factor"),within.col.names=c("cluster_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  factor.2.start.cluster <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.2.start.cluster",supressed_amount=cur.model$supressed))
  #factor.year.start.cluster
  cur.model<-calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp"),factor.calculation="as.POSIXlt(timestamp)$year",by.col.names = c("factor"),within.col.names=c("cluster_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  factor.3.start.cluster <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.3.start.cluster",supressed_amount=cur.model$supressed))
  #nofactor.start.cluster
  cur.model<-calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp"),factor.calculation=NA,by.col.names = c(),within.col.names=c("cluster_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  nofactor.start.cluster <- to.valid.matrix(cur.model$stats)
   supression_log <- rbind(supression_log,c(model="nofactor.start.cluster",supressed_amount=cur.model$supressed))


   #Calculates mean sequqence length
  mean.seq_duration <- dplyr::group_by(start.of.seq,seq_duration)
  mean.seq_duration <- dplyr::summarize(mean.seq_duration, mean = mean(seq_duration_numeric),sd=sd(seq_duration_numeric),freq=n())
  mean.seq_duration<- dplyr::ungroup(mean.seq_duration)
  mean.seq_duration$freq <-mean.seq_duration$freq/sum(mean.seq_duration$freq)
  mean.seq_duration[is.na(mean.seq_duration$sd),"sd"]<- 0
  mean.start_hour <- dplyr::group_by(start.of.seq,hour)
  mean.start_hour <- dplyr::summarize(mean.start_hour, freq=n())
  mean.start_hour<- dplyr::ungroup(mean.start_hour)
  mean.start_hour$freq <-mean.start_hour$freq/sum(mean.start_hour$freq)

    #factor.week.seq_duration
  #cur.model<- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation="as.integer(ceiling(as.POSIXlt(timestamp)$yday/4))",by.col.names = c("cluster_id","factor"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  cur.model<- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%w\"))",by.col.names = c("cluster_id","factor"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  factor.1.start_info <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.1.start_info",supressed_amount=cur.model$supressed))
  #factor.month.seq_duration
  cur.model<- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation="as.POSIXlt(timestamp)$mon",by.col.names = c("cluster_id","factor"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  factor.2.start_info <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.2.start_info",supressed_amount=cur.model$supressed))
  #factor.year.seq_duration
  cur.model <- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation="as.POSIXlt(timestamp)$year",by.col.names = c("cluster_id","factor"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  factor.3.start_info <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.3.start_info",supressed_amount=cur.model$supressed))
  #nofactor.seq_duration
  cur.model<- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation=NA,by.col.names = c("cluster_id"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  nofactor.start_info <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nofactor.start_info",supressed_amount=cur.model$supressed))
  #nocluster.factor.week.seq_duration
  #cur.model<- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation="as.integer(ceiling(as.POSIXlt(timestamp)$yday/4))",by.col.names = c("factor"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  cur.model<- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%w\"))",by.col.names = c("factor"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  nocluster.factor.1.start_info <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.1.start_info",supressed_amount=cur.model$supressed))
  #nocluster.factor.month.seq_duration
  cur.model<- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation="as.POSIXlt(timestamp)$mon",by.col.names = c("factor"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  nocluster.factor.2.start_info <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.2.start_info",supressed_amount=cur.model$supressed))
  #nocluster.factor.year.seq_duration
  cur.model <- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation="as.POSIXlt(timestamp)$year",by.col.names = c("factor"),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  nocluster.factor.3.start_info <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.3.start_info",supressed_amount=cur.model$supressed))
  #nocluster.nofactor.seq_duration
  cur.model<- calc.single.model(clustered.data=start.of.seq,relevant.col.names=c("objectid","seq_id","timestamp","hour","seq_duration","seq_duration_numeric"),factor.calculation=NA,by.col.names = c(),within.col.names=c("hour","seq_duration"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(seq_duration_numeric)","sd(seq_duration_numeric)"),statistics.col.names=c("sequences","objects","total","seq_duration_numeric","dur_sd"),eps=c_eps)
  nocluster.nofactor.start_info <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.nofactor.start_info",supressed_amount=cur.model$supressed))

  #factor.dur.starting.state
  cur.model<-calc.single.model(clustered.data=sequences,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","state_id","seq_duration"),factor.calculation="seq_duration",by.col.names = c("cluster_id","factor"),within.col.names=c("state_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  factor.1.starting.state <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.1.starting.state",supressed_amount=cur.model$supressed))

  #factor.hour.starting.state
  cur.model<-calc.single.model(clustered.data=sequences,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","state_id"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%H\"))",by.col.names = c("cluster_id","factor"),within.col.names=c("state_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  factor.2.starting.state <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.2.starting.state",supressed_amount=cur.model$supressed))

  #factor.weekday.starting.state
  cur.model<-calc.single.model(clustered.data=sequences,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","state_id"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%w\"))",by.col.names = c("cluster_id","factor"),within.col.names=c("state_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  factor.3.starting.state <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.3.starting.state",supressed_amount=cur.model$supressed))

  #nofactor.starting.state
  cur.model<-calc.single.model(clustered.data=sequences,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","state_id"),factor.calculation=NA,by.col.names = c("cluster_id"),within.col.names=c("state_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  nofactor.starting.state <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nofactor.starting.state",supressed_amount=cur.model$supressed))
  #nocluster.factor.dur.starting.state
  cur.model<-calc.single.model(clustered.data=sequences,relevant.col.names=c("objectid","seq_id","timestamp","state_id","seq_duration"),factor.calculation="seq_duration",by.col.names = c("factor"),within.col.names=c("state_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  nocluster.factor.1.starting.state <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.1.starting.state",supressed_amount=cur.model$supressed))

  #nocluster.factor.hour.starting.state
  cur.model<-calc.single.model(clustered.data=sequences,relevant.col.names=c("objectid","seq_id","timestamp","state_id"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%H\"))",by.col.names = c("factor"),within.col.names=c("state_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  nocluster.factor.2.starting.state <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.2.starting.state",supressed_amount=cur.model$supressed))

  #nocluster.factor.weekday.starting.state
  cur.model<-calc.single.model(clustered.data=sequences,relevant.col.names=c("objectid","seq_id","timestamp","state_id"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%w\"))",by.col.names = c("factor"),within.col.names=c("state_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  nocluster.factor.3.starting.state <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.3.starting.state",supressed_amount=cur.model$supressed))

  #nocluster.nofactor.starting.state
  cur.model<-calc.single.model(clustered.data=sequences,relevant.col.names=c("objectid","seq_id","timestamp","state_id"),factor.calculation=NA,by.col.names=c(),within.col.names=c("state_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  nocluster.nofactor.starting.state <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.nofactor.starting.state",supressed_amount=cur.model$supressed))

  #calculate cluster transitions
  cluster.transition <- start.of.seq
  dplyr::ungroup(cluster.transition)
  cluster.transition <- dplyr::group_by(cluster.transition,objectid)
  cluster.transition <- dplyr::arrange(cluster.transition,timestamp)
  cluster.transition <- dplyr::mutate(cluster.transition,next_id=lead(cluster_id),tbe=as.numeric((lead(timestamp)-timestamp),"secs"))
  cluster.transition <- cluster.transition[!is.na(cluster.transition$next_id),]
  dplyr::ungroup(cluster.transition)
  mean.cluster.tbe <- mean(cluster.transition$tbe)
  #factor.week.cluster.transition
  #cur.model<-calc.single.model(clustered.data=cluster.transition,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","next_id","tbe"),factor.calculation="as.integer(ceiling(as.POSIXlt(timestamp)$yday/4))",by.col.names = c("cluster_id","factor"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  cur.model<-calc.single.model(clustered.data=cluster.transition,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","next_id","tbe"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%w\"))",by.col.names = c("cluster_id","factor"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  factor.1.cluster.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.1.cluster.transition",supressed_amount=cur.model$supressed))
  #factor.month.cluster.transition
  cur.model<-calc.single.model(clustered.data=cluster.transition,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","next_id","tbe"),factor.calculation="as.POSIXlt(timestamp)$mon",by.col.names = c("cluster_id","factor"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  factor.2.cluster.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.2.cluster.transition",supressed_amount=cur.model$supressed))
  #factor.year.cluster.transition
  cur.model<-calc.single.model(clustered.data=cluster.transition,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","next_id","tbe"),factor.calculation="as.POSIXlt(timestamp)$year",by.col.names = c("cluster_id","factor"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  factor.3.cluster.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.3.cluster.transition",supressed_amount=cur.model$supressed))
  #nofactor.cluster.transition
  cur.model<-calc.single.model(clustered.data=cluster.transition,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","next_id","tbe"),factor.calculation=NA,by.col.names = c("cluster_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  nofactor.cluster.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nofactor.cluster.transition",supressed_amount=cur.model$supressed))

  sequences<- dplyr::ungroup(sequences)
  #sequences1<- dplyr::filter(sequences,next_state_id!="*")
  sequences1<- dplyr::filter(sequences,!is.na(next_state_id))

  names(sequences1)[names(sequences1)=="next_state_id"]<-"next_id"
  mean.transition.tbe<- sequences1[,c("cluster_id","tbe")]
  mean.transition.tbe<- dplyr::group_by(mean.transition.tbe,cluster_id)
  mean.transition.tbe<- dplyr::summarise(mean.transition.tbe,mean=mean(tbe))
  mean.transition.tbe<- to.valid.matrix(dplyr::ungroup(mean.transition.tbe))
  #factor.dur.transition
  cur.model<-calc.single.model(clustered.data=sequences1,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","state_id","next_id","seq_duration","tbe"),factor.calculation="seq_duration",by.col.names = c("cluster_id","factor","state_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  factor.1.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.1.transition",supressed_amount=cur.model$supressed))
  #factor.hour.transition
  cur.model<-calc.single.model(clustered.data=sequences1,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","state_id","next_id","tbe"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%H\"))",by.col.names = c("cluster_id","factor","state_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  factor.2.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.2.transition",supressed_amount=cur.model$supressed))
  #factor.weekday.transition
  cur.model<-calc.single.model(clustered.data=sequences1,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","state_id","next_id","tbe"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%w\"))",by.col.names = c("cluster_id","factor","state_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  factor.3.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="factor.3.transition",supressed_amount=cur.model$supressed))
  #nofactor.transition
  cur.model<-calc.single.model(clustered.data=sequences1,relevant.col.names=c("objectid","seq_id","cluster_id","timestamp","state_id","next_id","tbe"),factor.calculation=NA,by.col.names = c("cluster_id","state_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  nofactor.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nofactor.transition",supressed_amount=cur.model$supressed))
  #nocluster.factor.dur.transition
  cur.model<-calc.single.model(clustered.data=sequences1,relevant.col.names=c("objectid","seq_id","timestamp","state_id","next_id","seq_duration","tbe"),factor.calculation="seq_duration",by.col.names = c("factor","state_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  nocluster.factor.1.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.1.transition",supressed_amount=cur.model$supressed))
  #nocluster.factor.hour.transition
  cur.model<-calc.single.model(clustered.data=sequences1,relevant.col.names=c("objectid","seq_id","timestamp","state_id","next_id","tbe"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%H\"))",by.col.names = c("factor","state_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  nocluster.factor.2.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.2.transition",supressed_amount=cur.model$supressed))
  #nocluster.factor.weekday.transition
  cur.model<-calc.single.model(clustered.data=sequences1,relevant.col.names=c("objectid","seq_id","timestamp","state_id","next_id","tbe"),factor.calculation="as.numeric(format(as.POSIXlt(timestamp), \"%w\"))",by.col.names = c("factor","state_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  nocluster.factor.3.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.factor.3.transition",supressed_amount=cur.model$supressed))
  #nocluster.nofactor.transition
  cur.model<-calc.single.model(clustered.data=sequences1,relevant.col.names=c("objectid","seq_id","timestamp","state_id","next_id","tbe"),factor.calculation=NA,by.col.names = c("state_id"),within.col.names=c("next_id"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()","mean(tbe)","sd(tbe)"),statistics.col.names=c("sequences","objects","total","tbe","sd"),eps=c_eps)
  nocluster.nofactor.transition <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="nocluster.nofactor.transition",supressed_amount=cur.model$supressed))
  sequences1<- NULL

  # list all possible transitions between states, but not in sequential order (for corrections when no transition is found due to anonymization)
  common_states <- unique(sequences[,c("objectid","seq_id","state_id")])
  common_states <- dplyr::inner_join(common_states,common_states,by=c("seq_id" = "seq_id","objectid"="objectid"))
  cur.model<-calc.single.model(clustered.data=common_states,relevant.col.names=c("objectid","seq_id","state_id.x","state_id.y"),factor.calculation=NA,by.col.names = c("state_id.x"),within.col.names=c("state_id.y"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  common_states <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="common_states",supressed_amount=cur.model$supressed))

#   common_states <- dplyr::group_by(common_states,state_id.x,state_id.y)
#   common_states <- dplyr::summarise(common_states,sequences=n())
#   total_cnt <- sum(common_states$sequences)
#   orig_size <- dim(common_states)[1]
#   common_states$differ<-(common_states$sequences/total_cnt)/((common_states$sequences-1)/(total_cnt-1))
#   common_states<- dplyr::filter(common_states , differ<=exp(c_eps))
#   total_cnt <- sum(common_states$sequences) #renormalize after possible supression
#   common_states$sequences<-common_states$sequences/total_cnt
#   supression_log <- rbind(supression_log,c(model="common_states",supressed_amount=dim(common_states)[1]))
#

  common_clusters <- unique(cluster.transition[,c("objectid","seq_id","cluster_id")])
  common_clusters <- dplyr::inner_join(common_clusters,common_clusters,by=c("objectid"="objectid"))
  names(common_clusters)[2]<- "seq_id"
  cur.model<-calc.single.model(clustered.data=common_clusters,relevant.col.names=c("objectid","seq_id","cluster_id.x","cluster_id.y"),factor.calculation=NA,by.col.names = c("cluster_id.x"),within.col.names=c("cluster_id.y"),statistics.calculations =c("n_distinct(seq_id)","n_distinct(objectid)","n()"),statistics.col.names=c("sequences","objects","total"),eps=c_eps)
  common_clusters <- to.valid.matrix(cur.model$stats)
  supression_log <- rbind(supression_log,c(model="common_clusters",supressed_amount=cur.model$supressed))

#   common_clusters <- dplyr::group_by(common_clusters,cluster_id.x,cluster_id.y)
#   common_clusters <- dplyr::summarise(common_clusters,sequences=n())
#   total_cnt <- sum(common_clusters$sequences)
#   orig_size <- dim(common_clusters)[1]
#   common_clusters$differ<-(common_clusters$sequences/total_cnt)/((common_clusters$sequences-1)/(total_cnt-1))
#   common_clusters<- dplyr::filter(common_clusters , differ<=exp(c_eps))
#   total_cnt <- sum(common_clusters$sequences) #renormalize after possible supression
#   common_clusters$sequences<-common_clusters$sequences/total_cnt
#   supression_log <- rbind(supression_log,c(model="common_clusters",supressed_amount=dim(common_clusters)[1]))

  list(
    start.list=start.list,
    factor.1.start.cluster=factor.1.start.cluster,
    factor.2.start.cluster=factor.2.start.cluster,
    factor.3.start.cluster=factor.3.start.cluster,
    nofactor.start.cluster=nofactor.start.cluster,
    factor.1.start_info=factor.1.start_info,
    factor.2.start_info=factor.2.start_info,
    factor.3.start_info=factor.3.start_info,
    nofactor.start_info=nofactor.start_info,
    nocluster.factor.1.start_info=nocluster.factor.1.start_info,
    nocluster.factor.2.start_info=nocluster.factor.2.start_info,
    nocluster.factor.3.start_info=nocluster.factor.3.start_info,
    nocluster.nofactor.start_info=nocluster.nofactor.start_info,
    factor.1.starting.state=factor.1.starting.state,
    factor.2.starting.state=factor.2.starting.state,
    factor.3.starting.state=factor.3.starting.state,
    nofactor.starting.state=nofactor.starting.state,
    nocluster.factor.1.starting.state=nocluster.factor.1.starting.state,
    nocluster.factor.2.starting.state=nocluster.factor.2.starting.state,
    nocluster.factor.3.starting.state=nocluster.factor.3.starting.state,
    nocluster.nofactor.starting.state=nocluster.nofactor.starting.state,
    factor.1.cluster.transition=factor.1.cluster.transition,
    factor.2.cluster.transition=factor.2.cluster.transition,
    factor.3.cluster.transition=factor.3.cluster.transition,
    nofactor.cluster.transition=nofactor.cluster.transition,
    factor.1.transition=factor.1.transition,
    factor.2.transition=factor.2.transition,
    factor.3.transition=factor.3.transition,
    nofactor.transition=nofactor.transition,
    nocluster.factor.1.transition=nocluster.factor.1.transition,
    nocluster.factor.2.transition=nocluster.factor.2.transition,
    nocluster.factor.3.transition=nocluster.factor.3.transition,
    nocluster.nofactor.transition=nocluster.nofactor.transition,
    common_states=common_states,
    common_clusters=common_clusters,
    mean.cluster.tbe=mean.cluster.tbe,
    mean.transition.tbe=mean.transition.tbe,
    mean.seq_duration = mean.seq_duration,
    mean.start_hour = mean.start_hour,
    supression_log=supression_log)
    #generator_obj$
#     new.objects.per.date=new.objects.per.date
#     ,factor.1.starting.state=factor.1.starting.state
#     ,factor.2.starting.state=factor.2.starting.state
#     ,factor.3.starting.state=factor.3.starting.state
#     ,factor.1.transition=factor.1.transition
#     ,factor.2.transition=factor.2.transition
#     ,factor.3.transition=factor.3.transition
#     #   ,factor.dur.counts=factor.dur.counts
#     #   ,factor.hour.counts=factor.hour.counts
#     #   ,factor.weekday.counts=factor.weekday.counts
#     ,factor.3.cluster.transition=factor.3.cluster.transition
#     ,factor.1.cluster.transition=factor.1.cluster.transition
#     ,factor.2.cluster.transition=factor.2.cluster.transition
#     ,nofactor.starting.state=nofactor.starting.state
#     ,nofactor.cluster.transition=nofactor.cluster.transition
#     ,nofactor.transition=nofactor.transition
#     ,nofactor.start_info=nofactor.start_info
#     #   ,factor.week.counts=factor.week.counts
#     #   ,factor.month.counts=factor.month.counts
#     #   ,factor.year.counts=factor.year.counts
#     ,factor.2.start_info=factor.2.start_info
#     ,factor.1.start_info=factor.1.start_info
#     ,factor.3.start_info=factor.3.start_info
#     ,supression_log)
#   #generator_obj
}

