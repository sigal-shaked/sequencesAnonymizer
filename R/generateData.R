require(dplyr)
require (psych)
require(foreach)
require(doParallel)
require(Rcpp)

build_vector_of_sizes <- function(m1_size,m2_size,m3_size,is_inverse){
  if (!is_inverse){
    total <- m1_size+m2_size+m3_size
    if (total>0){
      w<- c(m1_size/total, m2_size/total,m3_size/total)
    } else {
      w<- c(0,0,0)
    }
  } else{
    m1_inv <- ifelse(m1_size>0,1/m1_size,0)
    m2_inv <- ifelse(m2_size>0,1/m2_size,0)
    m3_inv <- ifelse(m3_size>0,1/m3_size,0)
    total <- m1_inv+m2_inv+m3_inv
    if (total>0){
      w<- c(m1_inv/total, m2_inv/total,m3_inv/total)
    } else {
      w<- c(0,0,0)
    }
  }
  w
}

to.valid.matrix<- function(data){
  if(dim(data)[1]==1){
    tmp<-matrix(unlist(data),nrow=1)
    colnames(tmp)<- names(data)
  } else{
    tmp<- as.matrix(data)
  }
  tmp
}

R_sample_start_info<- function(m,factor.1,factor.2,factor.3,method,p_cluster_id=NULL){

  cur.nofactor.start_info<- to.valid.matrix(dplyr::filter(as.data.frame(m$nofactor.start_info),cluster_id==p_cluster_id)[,c(2,3,6,7,8)])
  if(dim(cur.nofactor.start_info)[1]>0 & !is.null(p_cluster_id)){
    cur.factor.1.start_info<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.1.start_info),cluster_id==p_cluster_id,factor==factor.1)[,c(3,4,7,8,9)])
    cur.factor.2.start_info<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.2.start_info),cluster_id==p_cluster_id,factor==factor.2)[,c(3,4,7,8,9)])
    cur.factor.3.start_info<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.3.start_info),cluster_id==p_cluster_id,factor==factor.3)[,c(3,4,7,8,9)])
    cur.m1<- ifelse(dim(cur.factor.1.start_info)[1]>0,dim(cur.factor.1.start_info)[1],0)
    cur.m2<- ifelse(dim(cur.factor.2.start_info)[1]>0,dim(cur.factor.2.start_info)[1],0)
    cur.m3<- ifelse(dim(cur.factor.3.start_info)[1]>0,dim(cur.factor.3.start_info)[1],0)

    c.w<- build_vector_of_sizes(m1_size=cur.m1,
                                m2_size=cur.m2,
                                m3_size=cur.m3,
                                is_inverse=method)

    cur.start.info <- sample_start_info(
      m1=cur.factor.1.start_info,
      m2=cur.factor.2.start_info,
      m3=cur.factor.3.start_info,
      m4=cur.nofactor.start_info,
      w=c.w)
  } else{
    cur.nocluster.nofactor.start_info<- to.valid.matrix(m$nocluster.nofactor.start_info[,c(1,2,5,6,7)])
    if(dim(cur.nocluster.nofactor.start_info)[1]>0){
      cur.nocluster.factor.1.start_info<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.1.start_info),factor==factor.1)[,c(2,3,6,7,8)])
      cur.nocluster.factor.2.start_info<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.2.start_info),factor==factor.2)[,c(2,3,6,7,8)])
      cur.nocluster.factor.3.start_info<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.3.start_info),factor==factor.3)[,c(2,3,6,7,8)])
      cur.m1<- ifelse(dim(cur.nocluster.factor.1.start_info)[1]>0,nrow(cur.nocluster.factor.1.start_info),0)
      cur.m2<- ifelse(dim(cur.nocluster.factor.2.start_info)[1]>0,nrow(cur.nocluster.factor.2.start_info),0)
      cur.m3<- ifelse(dim(cur.nocluster.factor.3.start_info)[1]>0,nrow(cur.nocluster.factor.3.start_info),0)

      c.w<- build_vector_of_sizes(m1_size=cur.m1,
                                  m2_size=cur.m2,
                                  m3_size=cur.m3,
                                  is_inverse=method)

      cur.start.info <- sample_start_info(
        m1=cur.nocluster.factor.1.start_info,
        m2=cur.nocluster.factor.2.start_info,
        m3=cur.nocluster.factor.3.start_info,
        m4=cur.nocluster.nofactor.start_info,
        w=c.w)
    } else{
      if(dim(m$mean.seq_duration)[1]>0){
        cur.index<- sample(1:dim(m$mean.seq_duration)[1],1,prob=m$mean.seq_duration$freq)
        cur.numeric_duration <- rnorm(1,mean = unlist(m$mean.seq_duration[cur.index,"mean"]),sd = unlist(m$mean.seq_duration[cur.index,"sd"]))
        cur.seq_duration <-m$mean.seq_duration[cur.index,"seq_duration"]
      } else{
        cur.numeric_duration <- rnorm(1,mean = m$mean.seq_duration[2],sd = m$mean.seq_duration[3])
        cur.seq_duration <- m$mean.seq_duration[1]
      }
      if(dim(m$mean.start_hour)[1]>0){
        cur.hour<- sample(m$mean.start_hour$hour,size=1,prob=m$mean.start_hour$freq)
      } else {
        cur.hour<- m$mean.start_hour[1]
      }
      cur.start.info <- c(cur.hour,cur.seq_duration,cur.numeric_duration)
    }
  }
}

R_sample_cluster<- function(m,factor.1,factor.2,factor.3,method,prev_id=NULL){

  #first cluster
  if(is.null(prev_id)){
    cur.nofactor.start.cluster<- to.valid.matrix(as.data.frame(m$nofactor.start.cluster)[,c(1,4)])
    cur.factor.1.start.cluster<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.1.start.cluster),factor==factor.1)[,c(1,5)])
    cur.factor.2.start.cluster<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.2.start.cluster),factor==factor.2)[,c(1,5)])
    cur.factor.3.start.cluster<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.3.start.cluster),factor==factor.3)[,c(1,5)])
    cur.m1<- ifelse(dim(cur.factor.1.start.cluster)[1]>0,dim(cur.factor.1.start.cluster)[1],0)
    cur.m2<- ifelse(dim(cur.factor.2.start.cluster)[1]>0,dim(cur.factor.2.start.cluster)[1],0)
    cur.m3<- ifelse(dim(cur.factor.3.start.cluster)[1]>0,dim(cur.factor.3.start.cluster)[1],0)

    c.w<- build_vector_of_sizes(m1_size=cur.m1,
                                m2_size=cur.m2,
                                m3_size=cur.m3,
                                is_inverse=method)


    if(dim(cur.nofactor.start.cluster)[1]>0){
      cur.sampled.cluster <- sample_starting_state(
        m1=cur.factor.1.start.cluster,
        m2=cur.factor.2.start.cluster,
        m3=cur.factor.3.start.cluster,
        m4=cur.nofactor.start.cluster,
        w=c.w)
    } else{
      if (dim(m$common_clusters)[1]>0){
        cur.sampled.cluster <- sample(unlist(unique(m$common_clusters[1])),size=1)
      } else{
        cur.sampled.cluster <- NA
      }
    }
  } else{
    cur.nofactor.cluster.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$nofactor.cluster.transition),cluster_id==prev_id)[,c(2,5,6,7)])
    cur.factor.1.cluster.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.1.cluster.transition),factor==factor.1,cluster_id==prev_id)[,c(3,6,7,8)])
    cur.factor.2.cluster.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.2.cluster.transition),factor==factor.2,cluster_id==prev_id)[,c(3,6,7,8)])
    cur.factor.3.cluster.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.3.cluster.transition),factor==factor.3,cluster_id==prev_id)[,c(3,6,7,8)])
    cur.m1<- ifelse(dim(cur.factor.1.cluster.transition)[1]>0,dim(cur.factor.1.cluster.transition)[1],0)
    cur.m2<- ifelse(dim(cur.factor.2.cluster.transition)[1]>0,dim(cur.factor.2.cluster.transition)[1],0)
    cur.m3<- ifelse(dim(cur.factor.3.cluster.transition)[1]>0,dim(cur.factor.3.cluster.transition)[1],0)

    c.w<- build_vector_of_sizes(m1_size=cur.m1,
                                m2_size=cur.m2,
                                m3_size=cur.m3,
                                is_inverse=method)
    if(dim(cur.nofactor.cluster.transition)[1]>0){
      cur.sampled.cluster <- sample_transition(
        m1=cur.factor.1.cluster.transition,
        m2=cur.factor.2.cluster.transition,
        m3=cur.factor.3.cluster.transition,
        m4=cur.nofactor.cluster.transition,
        w=c.w)
    } else{
      c_cluster_id = NA
      cur.near <- dplyr::filter(as.data.frame(to.valid.matrix(m$common_clusters)),cluster_id.x ==prev_id)[,c(2,3)]

      if(dim(cur.near)[1]>1){
        c_cluster_id <- sample(as.character(cur.near[,1]),size=1,prob=as.numeric(cur.near[,2]))
      } else if(dim(cur.near)[1]==1){
        c_cluster_id <-as.character(cur.near[1])
      }
      if (!is.null(c_cluster_id)){
        if (!is.null(prev_id)){
          cur.sampled.cluster <- c(c_cluster_id,as.numeric(m$mean.cluster.tbe))
        }
        else{
          cur.sampled.cluster <- c(c_cluster_id)
        }
      }
    }
  }
  cur.sampled.cluster
}

R_sample_state<- function(m,p_cluster_id,factor.1,factor.2,factor.3,method,prev_id=NULL){
  #first state
  if(is.null(prev_id)){
    cur.nofactor.starting.state<- to.valid.matrix(dplyr::filter(as.data.frame(m$nofactor.starting.state),cluster_id==p_cluster_id)[,c(2,5)])
    if(dim(cur.nofactor.starting.state)[1]>0 & !is.null(p_cluster_id)){
      cur.factor.1.starting.state<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.1.starting.state),cluster_id==p_cluster_id,factor==factor.1)[,c(3,6)])
      cur.factor.2.starting.state<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.2.starting.state),cluster_id==p_cluster_id,factor==factor.2)[,c(3,6)])
      cur.factor.3.starting.state<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.3.starting.state),cluster_id==p_cluster_id,factor==factor.3)[,c(3,6)])
      cur.m1<- ifelse(dim(cur.factor.1.starting.state)[1]>0,dim(cur.factor.1.starting.state)[1],0)
      cur.m2<- ifelse(dim(cur.factor.2.starting.state)[1]>0,dim(cur.factor.2.starting.state)[1],0)
      cur.m3<- ifelse(dim(cur.factor.3.starting.state)[1]>0,dim(cur.factor.3.starting.state)[1],0)

      c.w<- build_vector_of_sizes(m1_size=cur.m1,
                                  m2_size=cur.m2,
                                  m3_size=cur.m3,
                                  is_inverse=method)

      cur.state <- sample_starting_state(
        m1=cur.factor.1.starting.state,
        m2=cur.factor.2.starting.state,
        m3=cur.factor.3.starting.state,
        m4=cur.nofactor.starting.state,
        w=c.w)
    } else{
      cur.nocluster.nofactor.starting.state<- to.valid.matrix(as.data.frame(m$nocluster.nofactor.starting.state[,c(1,4)]))
      if(dim(cur.nocluster.nofactor.starting.state)[1]>0){
        cur.nocluster.factor.1.starting.state<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.1.starting.state),factor==factor.1)[,c(2,5)])
        cur.nocluster.factor.2.starting.state<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.2.starting.state),factor==factor.2)[,c(2,5)])
        cur.nocluster.factor.3.starting.state<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.3.starting.state),factor==factor.3)[,c(2,5)])
        cur.m1<- ifelse(dim(cur.nocluster.factor.1.starting.state)[1]>0,nrow(cur.nocluster.factor.1.starting.state),0)
        cur.m2<- ifelse(dim(cur.nocluster.factor.2.starting.state)[1]>0,nrow(cur.nocluster.factor.2.starting.state),0)
        cur.m3<- ifelse(dim(cur.nocluster.factor.3.starting.state)[1]>0,nrow(cur.nocluster.factor.3.starting.state),0)

        c.w<- build_vector_of_sizes(m1_size=cur.m1,
                                    m2_size=cur.m2,
                                    m3_size=cur.m3,
                                    is_inverse=method)

        cur.state <- sample_starting_state(
          m1=cur.nocluster.factor.1.starting.state,
          m2=cur.nocluster.factor.2.starting.state,
          m3=cur.nocluster.factor.3.starting.state,
          m4=cur.nocluster.nofactor.starting.state,
          w=c.w)
      } else{cur.state <-sample(unlist(unique(m$common_states[1])),1)}
    }
    #not the starting state
  } else{
    cur.nofactor.transition<- to.valid.matrix((dplyr::filter(as.data.frame(m$nofactor.transition),cluster_id==p_cluster_id,state_id==prev_id))[,c(3,6,7,8)])
    if(dim(cur.nofactor.transition)[1]>0 & !is.null(p_cluster_id)){
      cur.factor.1.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.1.transition),cluster_id==p_cluster_id,factor==factor.1,state_id==prev_id)[,c(4,7,8,9)])
      cur.factor.2.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.2.transition),cluster_id==p_cluster_id,factor==factor.2,state_id==prev_id)[,c(4,7,8,9)])
      cur.factor.3.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$factor.3.transition),cluster_id==p_cluster_id,factor==factor.3,state_id==prev_id)[,c(4,7,8,9)])
      cur.m1<- ifelse(dim(cur.factor.1.transition)[1]>0,dim(cur.factor.1.transition)[1],0)
      cur.m2<- ifelse(dim(cur.factor.2.transition)[1]>0,dim(cur.factor.2.transition)[1],0)
      cur.m3<- ifelse(dim(cur.factor.3.transition)[1]>0,dim(cur.factor.3.transition)[1],0)

      c.w<- build_vector_of_sizes(m1_size=cur.m1,
                                  m2_size=cur.m2,
                                  m3_size=cur.m3,
                                  is_inverse=method)

      cur.state <- sample_transition(
        m1=cur.factor.1.transition,
        m2=cur.factor.2.transition,
        m3=cur.factor.3.transition,
        m4=cur.nofactor.transition,
        w=c.w)
    } else{
      cur.nocluster.nofactor.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.nofactor.transition),state_id==prev_id)[,c(2,5,6,7)])
      if(dim(cur.nocluster.nofactor.transition)[1]>0){
        cur.nocluster.factor.1.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.1.transition),factor==factor.1,state_id==prev_id)[,c(3,6,7,8)])
        cur.nocluster.factor.2.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.2.transition),factor==factor.2,state_id==prev_id)[,c(3,6,7,8)])
        cur.nocluster.factor.3.transition<- to.valid.matrix(dplyr::filter(as.data.frame(m$nocluster.factor.3.transition),factor==factor.3,state_id==prev_id)[,c(3,6,7,8)])
        cur.m1<- ifelse(dim(cur.nocluster.factor.1.transition)[1]>0,nrow(cur.nocluster.factor.1.transition),0)
        cur.m2<- ifelse(dim(cur.nocluster.factor.2.transition)[1]>0,nrow(cur.nocluster.factor.2.transition),0)
        cur.m3<- ifelse(dim(cur.nocluster.factor.3.transition)[1]>0,nrow(cur.nocluster.factor.3.transition),0)

        c.w<- build_vector_of_sizes(m1_size=cur.m1,
                                    m2_size=cur.m2,
                                    m3_size=cur.m3,
                                    is_inverse=method)

        cur.state <- sample_transition(
          m1=cur.nocluster.factor.1.transition,
          m2=cur.nocluster.factor.2.transition,
          m3=cur.nocluster.factor.3.transition,
          m4=cur.nocluster.nofactor.transition,
          w=c.w)
      } else{cur.state<-c(-1,0)}
    }
  }
  if(cur.state[1]==-1){
    c.state_id <- -1
    cur.near <- dplyr::filter(as.data.frame(to.valid.matrix(m$common_states)),state_id.x ==prev_id)[,c(2,3)]
    if(dim(cur.near)[1]>1){
      c.state_id <- sample(cur.near[,1],size=1,prob=cur.near[,2])
    } else if(dim(cur.near)[1]==1){
      c.state_id<-cur.near[1]
    }
    if (c.state_id!=-1){
      cur.tbe<-NULL
      if (!is.null(prev_id)){
        cur.sampled.tbe<- dplyr::filter(as.data.frame(m$mean.transition.tbe),cluster_id==p_cluster_id)[2]
        if(nrow(cur.sampled.tbe)==0){
          cur.tbe<- mean(m$mean.transition.tbe[,2])
        } else if(nrow(cur.sampled.tbe)==1){
          cur.tbe<-cur.sampled.tbe
        }
        cur.state <- c(c.state_id,as.numeric(cur.tbe))
      }
      else{
        cur.state <- c(c.state_id)
      }
    }
  }
  cur.state
}

sample.sequence.chain<- function(m,p_objID,p_seqID,p_cluster_id,p_seq_duration,p_start_time,p_end_time,p_method){
  #seq_duration
  cur.factor.1 <- p_seq_duration
  #cur.hour
  cur.factor.2 <-  as.integer(as.POSIXlt(p_start_time,origin = "1970-01-01")$hour)
  #cur.weekday
  cur.factor.3 <- as.integer(as.POSIXlt(p_start_time,origin = "1970-01-01")$wday)
  cur.timestamp <- p_start_time
  cur.state <- R_sample_state(m,p_cluster_id=p_cluster_id,factor.1=cur.factor.1,factor.2=cur.factor.2,factor.3=cur.factor.3,method=p_method)
  if(cur.state!=-1){
    cur.sequence <- as.matrix(c(p_objID,p_seqID,as.character(cur.timestamp),cur.state))
    dim(cur.sequence)<-c(1,4)
    names(cur.sequence)<-c("objectid","seq_id","timestamp","state_id")
    cur.sequence<- to.valid.matrix(cur.sequence)

    while(cur.timestamp<=p_end_time){
      cur.transition <- R_sample_state(m,p_cluster_id=p_cluster_id,factor.1=cur.factor.1,factor.2=cur.factor.2,factor.3=cur.factor.3,method=p_method,prev_id=cur.state)
      cur.state <- cur.transition[1]
      if(cur.state==-1){
        break
      }
      cur.tbe <- as.numeric(cur.transition[2])
      if(cur.tbe<=2){
        cur.tbe <-  sample(c(2,60,3600),size=1,replace = T,prob=c(0.5,0.3,0.2))
      }
      cur.timestamp <- cur.timestamp + cur.tbe
      cur.sequence<-rbind(cur.sequence,c(p_objID,p_seqID,as.character(cur.timestamp),cur.state))
    }
  }
  to.valid.matrix(cur.sequence)
}

sample.meta.chain<- function(m,p_method){
  sequences.frame<-NULL
  end.timestamp<-NULL
  for(i in 1:dim(m$start.list)[1]){
    for(cur.obj in 1:m$start.list$objects[i]){
      cur.date <- m$start.list$date[i]
      #cur.weekday
      cur.factor.1<- as.integer(as.numeric(format(as.POSIXlt(cur.date,origin = "1970-01-01"), "%w")))
      #cur.month
      cur.factor.2<- as.integer(as.POSIXlt(cur.date,origin = "1970-01-01")$mon)
      #cur.year
      cur.factor.3<- as.integer(as.POSIXlt(cur.date,origin = "1970-01-01")$year)
      mean.num.seq <- as.numeric(unlist(m$start.list$mean.num.seq[i]))
      sd.num.seq <- as.numeric(unlist(m$start.list$sd.num.seq[i]    ))
      cur.num.seq <- max(1,round(rnorm(n=1, mean=mean.num.seq, sd=sd.num.seq),0))
      start.timestamp<- as.POSIXlt(paste(cur.date,"00:00:01",sep=""),origin = "1970-01-01")
      for (cur.seq in 1:cur.num.seq){
        if (cur.seq==1){
          c.cluster_id <- R_sample_cluster(m,factor.1=cur.factor.1,factor.2=cur.factor.2,factor.3=cur.factor.3,method=p_method,prev_id=NULL)
        } else{
          cur.cluster_id <- R_sample_cluster(m,factor.1=cur.factor.1,factor.2=cur.factor.2,factor.3=cur.factor.3,method=p_method,prev_id=c.cluster_id )
          cur.cluster_tbe <- as.numeric(cur.cluster_id[2])
          c.cluster_id <- cur.cluster_id[1]
        }
        cur.start_info <- R_sample_start_info(m,factor.1=cur.factor.1,factor.2=cur.factor.2,factor.3=cur.factor.3,method=p_method,p_cluster_id=c.cluster_id)
        cur.hour<-cur.start_info[1]
        if (cur.seq==1){
          start.timestamp<- as.POSIXlt(paste(substr(start.timestamp,1,11),cur.hour,":00:00",sep=""),origin = "1970-01-01")+runif(1,1,3599)
        } else{
          start.timestamp<- start.timestamp + max(cur.cluster_tbe,cur.seq.duration.numeric+1)
          start.timestamp<- as.POSIXlt(paste(substr(start.timestamp,1,11),cur.hour,substr(start.timestamp,14,19),sep=""),origin = "1970-01-01")
        }
        if (!is.null(end.timestamp)){
          while( start.timestamp<= end.timestamp){
            start.timestamp<- start.timestamp + 60*60*24 #add one day
          }
        }
        #cur.weekday
        cur.factor.1<- as.integer(as.numeric(format(as.POSIXlt(start.timestamp,origin = "1970-01-01"), "%w")))
        #cur.month
        cur.factor.2<- as.integer(as.POSIXlt(start.timestamp,origin = "1970-01-01")$mon)
        #cur.year
        cur.factor.3<- as.integer(as.POSIXlt(start.timestamp,origin = "1970-01-01")$year)

        cur.seq.duration<-cur.start_info[2]
        cur.seq.duration.numeric<-as.numeric(cur.start_info[3])
        end.timestamp <-  start.timestamp +cur.seq.duration.numeric
        if(is.null(sequences.frame)){
          sequences.frame <- as.matrix(c(i,cur.obj,cur.seq,c.cluster_id,cur.seq.duration,as.character(start.timestamp),as.character(end.timestamp)))
          dim(sequences.frame)<- c(1,7)
          names(sequences.frame)<- c("dateNo","objectID","seqID","clusterID","seq_duration","startTime","endTime")
          sequences.frame<-to.valid.matrix(sequences.frame)
        } else{
          sequences.frame<- rbind(sequences.frame,c(i,cur.obj,cur.seq,c.cluster_id,cur.seq.duration,as.character(start.timestamp),as.character(end.timestamp)))
        }
      }
    }
  }
  sequences.frame[,c("objectID")]<- as.numeric(factor(paste0(sequences.frame[,c("dateNo")],"-",sequences.frame[,c("objectID")])))
  sequences.frame[,c("seqID")]<- as.numeric(factor(paste0(sequences.frame[,c("dateNo")],"-",sequences.frame[,c("objectID")],"-",sequences.frame[,c("seqID")])))
  apply(sequences.frame,2,unlist)
}


generate.synthetic.data <- function(m,p.method){
  doParallel::registerDoParallel(parallel::detectCores())
  sequences.frame<-sample.meta.chain(m,p.method)
  require(foreach)

  synthetic.data<-foreach(i = 1:dim(sequences.frame)[1],.packages=c("dplyr","psych","Rcpp","sequencesAnonymizer"),.combine="rbind") %dopar% {
    #for(i in 1:dim(sequences.frame)[1])  {
    cur.cluster_id <- as.integer(sequences.frame[i,"clusterID"])
    cur.seq_duration <- as.integer(sequences.frame[i,"seq_duration"])
    cur.objectID <- sequences.frame[i,"objectID"]
    cur.seqID <- as.integer(sequences.frame[i,"seqID"])
    cur.startTime <- as.POSIXlt(unlist(sequences.frame[i,"startTime"]),origin = "1970-01-01")
    cur.endTime <- as.POSIXlt(unlist(sequences.frame[i,"endTime"]),origin = "1970-01-01")
    cur.bulk<-  sample.sequence.chain(m,
                                      p_objID=cur.objectID,
                                      p_seqID=cur.seqID,
                                      p_cluster_id=cur.cluster_id,
                                      p_seq_duration=cur.seq_duration,
                                      p_start_time=cur.startTime,
                                      p_end_time=cur.endTime,
                                      p_method=p.method)
    if(dim(cur.bulk)[1]==0){
      cur.bulk<-c("*","*","*","*")
      dim(cur.bulk)<-c(1,4)
      names(cur.bulk)<-c("objectid","seq_id","timestamp","state_id")
      cur.sequence<- to.valid.matrix(cur.bulk)
    }
    cur.bulk
  }

  foreach::registerDoSEQ()
  synthetic.data<-as.data.frame(synthetic.data,stringsAsFactors = F)
  synthetic.data<-apply(synthetic.data,2,unlist)
  synthetic.data<-synthetic.data[synthetic.data[,"objectid"]!='*',]
  synthetic.data[order(synthetic.data[,"objectid"],synthetic.data[,"timestamp"]),]
}


