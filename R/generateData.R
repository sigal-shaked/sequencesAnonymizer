
require(dplyr)
require (psych)
require(foreach)
require(doParallel)
require(Rcpp)
#setwd('D:\\sigal\\doctorat\\thesis\\sequqencesGenerator\\')
#D:\sigal\doctorat\thesis\sequencesGenerator
#D:/sigal/doctorat/thesis/sequencesGenerator/sequencesGenerator/src/sample_transition.cpp
#sourceCpp("src\\sample_transition.cpp")
#sourceCpp("src\\sample_start_info.cpp")
#sourceCpp("src\\sample_starting_state.cpp")
#sourceCpp("src\\extract_models_sizes.cpp")

#setwd('D:\\sigal\\doctorat\\thesis')




#       .First.lib <- function(lib, pkg)
#       {
#         library.dynam("LIBNAME", pkg, lib)
#       }
#dyn.load("extract_models_sizes.cpp")
#sourceCpp("code\\sequence-generator\\c\\sample_methods.cpp")
#importFrom(Rcpp, sourceCpp)
#compileAttributes()

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


# # generates a sequence and write it into the output file
# generate.transitions <- function(m,p.obj,p.seqid,p.cluster_id,p.seq_duration,p.timestamp,p.state,p.objects,p.freq,end.timestamp,method){#,min.objects,max.freq){
#   cur.cluster_id <- p.cluster_id
#   cur.seq_duration <- p.seq_duration
#   cur.timestamp <- p.timestamp
#   cur.state <- p.state
#
#   cur.sequence <- as.matrix(t(c(p.obj,p.seqid, as.character(as.POSIXct(cur.timestamp,origin = "1970-01-01")), cur.state,0)))
#   names(cur.sequence)<- c("objectid","seqid","timestamp","stateid","tbe")
#
#   #if(p.objects < min.objects || p.freq> max.freq){
#   cur.sequence<- as.matrix(cur.sequence[0,])
#   #}
#   prev.states <- as.matrix(t(c(cur.state,0)))[0,] #for sampling previous states when no transition is available
#   cntr = 1 # for counting states in the current sequence
#   c.prev_id <- cur.state
#   while (!is.na(cur.timestamp) & cur.timestamp<=end.timestamp){
#
#     cur.hour <-  as.integer(as.POSIXlt(cur.timestamp,origin = "1970-01-01")$hour)
#     cur.weekday <- as.integer(as.POSIXlt(cur.timestamp,origin = "1970-01-01")$wday)
#
#     #if(nrow(m$factor.dur.counts)>1){
#     if(nrow(m$factor.1.transition)>0){
#       #cur.m1<- m$factor.dur.counts[,1:3]
#       #cur.m1<- max(m$factor.1.transition[m$factor.1.transition[,"cluster_id"]==cur.cluster_id & m$factor.1.transition[,"factor"]==cur.seq_duration & m$factor.1.transition[,"state_id"]==c.prev_id,"sequences"])
#       cur.m1<-select_(dplyr::filter(as.data.frame(m$factor.1.transition),cluster_id==cur.cluster_id,factor==cur.seq_duration,state_id==c.prev_id),"sequences")
#       cur.m1<- ifelse(nrow(cur.m1)>0,max(cur.m1),0)
#     } else{cur.m1=0}
#     #     } else{
#     #       #cur.m1<- matrix(unlist(m$factor.dur.counts[,1:3]),nrow=1)
#     #     }
#     if(nrow(m$factor.2.transition)>0){
#       #cur.m2<- max(m$factor.2.transition[m$factor.2.transition[,"cluster_id"]==cur.cluster_id & m$factor.2.transition[,"factor"]==cur.hour & m$factor.2.transition[,"state_id"]==c.prev_id,"sequences"])
#       cur.m2<-select_(dplyr::filter(as.data.frame(m$factor.2.transition),cluster_id==cur.cluster_id,factor==cur.hour,state_id==c.prev_id),"sequences")
#       cur.m2<- ifelse(nrow(cur.m2)>0,max(cur.m2),0)
#     }else{cur.m2=0}
#     #     if(nrow(m$factor.hour.counts)>1){
#     #       cur.m2<- matrix(unlist(m$factor.hour.counts[,1:3]),nrow=1)
#     #     } else{
#     #       cur.m2<- matrix(m$factor.hour.counts,nrow=1,ncol=3)
#     #     }
#     if(nrow(m$factor.3.transition)>0){
#       #cur.m2<- max(m$factor.3.transition[m$factor.3.transition[,"cluster_id"]==cur.cluster_id & m$factor.3.transition[,"factor"]==cur.weekday & m$factor.3.transition[,"state_id"]==c.prev_id,"sequences"])
#       cur.m3<-select_(dplyr::filter(as.data.frame(m$factor.3.transition),cluster_id==cur.cluster_id,factor==cur.weekday,state_id==c.prev_id),"sequences")
#       cur.m3<- ifelse(nrow(cur.m3)>0,max(cur.m3),0)
#     }else{cur.m3=0}
#     #     if(nrow(m$factor.weekday.counts)>1){
#     #       cur.m3<- m$factor.weekday.counts[,1:3]
#     #     } else{
#     #       cur.m3<- matrix(unlist(m$factor.weekday.counts[,1:3]),nrow=1)
#     #     }
#
#     #     c.w<- extract_models_sizes( m1=cur.m1,
#     #                               m2=cur.m2,
#     #                               m3=cur.m3,
#     #                               cur_cluster_id=cur.cluster_id ,
#     #                               cur_factor1=cur.seq_duration,
#     #                               cur_factor2=cur.hour,
#     #                               cur_factor3=cur.weekday,
#     #                               inverse_method=method)
#     c.w<- build_vector_of_sizes(m1_size=cur.m1,
#                                 m2_size=cur.m2,
#                                 m3_size=cur.m3,
#                                 is_inverse=method)
#
#     cur.state.transition <- sample_transition(  m1=m$factor.1.transition,
#                                                 m2=m$factor.2.transition,
#                                                 m3= m$factor.3.transition,
#                                                 m4=m$nofactor.transition,
#                                                 cur_cluster_id=cur.cluster_id,
#                                                 cur_factor1=cur.seq_duration,
#                                                 cur_factor2=cur.hour,
#                                                 cur_factor3=cur.weekday,
#                                                 prev_id=c.prev_id,
#                                                 w=c.w)
#     if(is.na(as.integer(cur.state.transition[1])) | is.na(cur.state.transition[2]) | is.na(cur.state.transition[3])){
#       #print(paste0(c.prev_id,"->cur.state:",cur.state))
#       break
#     }
#     ##//cur.state.transition <- sample.transition(m,type=0,p.cluster_id=cur.cluster_id,p.factor1=cur.seq_duration,p.factor2=cur.hour,p.factor3=cur.weekday,p.c.prev_id=c.prev_id,method)
#     #cur.tbe <- ifelse(is.numeric(cur.state.transition[2]),as.numeric(cur.state.transition[2]),60)
#     cur.tbe <- cur.state.transition[2]
#     if(cur.tbe<=2){
#       cur.tbe <-  sample(c(0,60,3600),size=1,replace = T,prob=c(0.5,0.3,0.2))
#     }
#     cur.timestamp <- cur.timestamp  + cur.tbe
#     cur.state <- as.integer(cur.state.transition[1])
#     if ( cur.state== -1 | is.na(as.integer(cur.state.transition[1])) | is.na(cur.state.transition[2]) | is.na(cur.state.transition[3])){
#       cur.state.transition<- sample_prev_states(prev.states)#,min.objects,max.freq)
#       if ( cur.state.transition[1]== -1){break}
#     }
#
#     #     if ( cur.state== -1 & !is.null(cur.sequence)){
#     #       if(dim(cur.sequence)[1]>1){
#     #         prev.states<- cur.sequence[1:(dim(cur.sequence)[1]-1),]
#     #         if(is.vector(prev.states)){
#     #           cur.state <- as.integer(unlist(prev.states[4]))
#     # #           if (!is.integer(cur.state)){
#     # #             print("non integer inner1.transition")
#     # #           }
#     #           cur.tbe <- as.numeric(unlist(prev.states[6]))
#     #         } else{
#     #           ids<- c(1:dim(prev.states)[1])
#     #           freqs <- ids/sum(ids)
#     #           id <- sample(ids,1,prob=freqs)
#     #           cur.state <- as.numeric(unlist(prev.states[id,4]))
#     # #           if (!is.integer(cur.state)){
#     # #             print("non integer inner2.transition")
#     # #           }
#     #           cur.tbe <- as.numeric(unlist(prev.states[id,6]))
#     #         }
#     #         #cur.tbe <- ifelse(is.numeric(cur.tbe),as.numeric(cur.tbe),0)
#     #         #cur.tbe <- as.numeric(cur.tbe)
#     #         if(cur.tbe==0){
#     #           cur.tbe <-  sample(c(0,60,3600),size=1,replace = T,prob=c(0.5,0.3,0.2))
#     #         }
#     #         cur.timestamp <- cur.timestamp  + cur.tbe
#     #       } else{break}
#     #     }
#     if (!is.na(as.character(as.POSIXct(cur.timestamp,origin = "1970-01-01"))) &!is.na(cur.state) & cur.timestamp<=end.timestamp){
#       #if ( cur.state.transition[3]>min.objects & cur.state.transition[4]<=max.freq){
#       cur.sequence <- rbind(cur.sequence, unlist(c(p.obj,p.seqid,as.character(as.POSIXct(cur.timestamp,origin = "1970-01-01")),  cur.state ,cur.tbe)))
#       prev.states <- rbind(prev.states,as.numeric(cur.sequence[cntr,4:5]))
#       cntr=cntr+1
#       # }
#       c.prev_id <- cur.state
#     }
#   }
#   cur.sequence
# }

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
    } else{NULL}
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

    cur.sampled.cluster <- sample_starting_state(
      m1=cur.factor.1.start.cluster,
      m2=cur.factor.2.start.cluster,
      m3=cur.factor.3.start.cluster,
      m4=cur.nofactor.start.cluster,
      w=c.w)
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

    cur.sampled.cluster <- sample_transition(
      m1=cur.factor.1.cluster.transition,
      m2=cur.factor.2.cluster.transition,
      m3=cur.factor.3.cluster.transition,
      m4=cur.nofactor.cluster.transition,
      w=c.w)
  }
  if(cur.sampled.cluster[1]==-1){
    c_cluster_id = -1
    cur.near <- dplyr::filter(as.data.frame(to.valid.matrix(m$common_clusters)),cluster_id.x ==prev_id)[,c(2,3)]

    if(dim(cur.near)[1]>1){
      c_cluster_id <- sample(as.character(cur.near[,1]),size=1,prob=as.numeric(cur.near[,2]))
    } else if(dim(cur.near)[1]==1){
      c_cluster_id <-as.character(cur.near[1])
    }
    if (c_cluster_id!=-1){
      if (!is.null(prev_id)){
        cur.sampled.cluster <- c(c_cluster_id,as.numeric(m$mean.cluster.tbe))
      }
      else{
        cur.sampled.cluster <- c(c_cluster_id)
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
      cur.nocluster.nofactor.starting.state<- to.valid.matrix(m$nocluster.nofactor.starting.state[,c(1,4)])
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
      } else{NULL}
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
#         if(is.null(cur.tbe)){
#           cur.tbe<- 0
#         }
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
#       if(is.null(cur.tbe)){
#         cur.tbe<-0
#       }
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
      #cur.factor.1<- as.integer(ceiling(as.POSIXlt(cur.date,origin = "1970-01-01")$yday/4))        cur.factor.1<- as.integer(as.numeric(format(as.POSIXlt(start.timestamp), "%w")))
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
        #c.cluster_id  = cur.cluster_id[,1]
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
        #cur.factor.1<- as.integer(ceiling(as.POSIXlt(start.timestamp,origin = "1970-01-01")$yday/4))
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
  sequences.frame
}
#
# calc.sequences.per.object <- function (m,c.timestamp,c.num.seq,c.cluster_id,c.factor1,c.factor2,c.factor3,method) {
#
#   cur.timestamp <- c.timestamp
#   cur.cluster_id <- c.cluster_id
#   cur.factor1 <- c.factor1
#   cur.factor2 <- c.factor2
#   cur.factor3 <- c.factor3
#   sequences<- NULL
#   end.timestamp <- NULL
#   #cur.weekday<- as.POSIXlt(cur.timestamp,origin = "1970-01-01")$wday
#   for (i in 1:c.num.seq){
#     if (cur.cluster_id==-1 ){ #if no cluster transition was found, sample cluster among previous clusters
#       if(!is.null(sequences) & dim(sequences)[1]>1){
#         prev.clusters<- sequences[1:(dim(sequences)[1]-1),]
#         if(is.vector(prev.clusters)){
#           cur.cluster_id <- unlist(prev.clusters[1])
#           cur.cluster_tbe <- as.numeric(unlist(prev.clusters[5]))
#         } else {
#           ids<- c(1:dim(prev.clusters)[1])
#           freqs <- ids/sum(ids)
#           id <- ifelse(dim(prev.clusters)[1]==1,ids,sample(ids,1,prob=freqs))
#           cur.cluster_id <- unlist(prev.clusters[id,1])
#           cur.cluster_tbe <- as.numeric(unlist(prev.clusters[id,5]))
#         }
#       } else{break}
#     } else{ #sample cluster and tbe
#       prev.cluster_id <- cur.cluster_id
#
#       if(nrow(m$factor.1.cluster.transition)>0){
#         cur.m1<- select_(dplyr::filter(as.data.frame(m$factor.1.cluster.transition),cluster_id==prev.cluster_id,factor==cur.factor1),"sequences")
#         cur.m1<- ifelse(nrow(cur.m1)>0,max(cur.m1),0)
#       }else{cur.m1=0}
#
#       if(nrow(m$factor.2.cluster.transition)>0){
#         cur.m2<- select_(dplyr::filter(as.data.frame(m$factor.2.cluster.transition),cluster_id==prev.cluster_id,factor==cur.factor2),"sequences")
#         cur.m2<- ifelse(nrow(cur.m2)>0,max(cur.m2),0)
#       }else{cur.m2=0}
#
#       if(nrow(m$factor.3.cluster.transition)>0){
#         cur.m3<- select_(dplyr::filter(as.data.frame(m$factor.3.cluster.transition),cluster_id==prev.cluster_id,factor==cur.factor3),"sequences")
#         cur.m3<- ifelse(nrow(cur.m3)>0,max(cur.m3),0)
#       }else{cur.m3=0}
#
#       c.w<- build_vector_of_sizes(m1_size=cur.m1,
#                                   m2_size=cur.m2,
#                                   m3_size=cur.m3,
#                                   is_inverse=method)
#
#       #             if(nrow(m$factor.week.counts)>1){
#       #               cur.m1<- as.matrix(m$factor.week.counts[,1:3])
#       #             } else{
#       #               #cur.m1<- rbind(c(-1,-1,-1),matrix(m$factor.week.counts[1:3],nrow=1))
#       #               cur.m1<- matrix(unlist(m$factor.week.counts[,1:3]),nrow=1)
#       #             }
#       #             if(nrow(m$factor.month.counts)>1){
#       #               cur.m2<- as.matrix(m$factor.month.counts[,1:3])
#       #             } else{
#       #               cur.m2<- matrix(unlist(m$factor.month.counts[,1:3]),nrow=1)
#       #
#       #              # cur.m2<- as.matrix(t(m$factor.month.counts[1:3]))
#       #             }
#       #             if(nrow(m$factor.year.counts)>1){
#       #               cur.m3<- as.matrix(m$factor.year.counts[,1:3])
#       #             } else{
#       #               cur.m3<- matrix(unlist(m$factor.year.counts[,1:3]),nrow=1)
#       #             }
#       #
#       #       c.w<- extract_models_sizes( m1=cur.m1,
#       #                                 m2=cur.m2,
#       #                                 m3=cur.m3,
#       #                                 cur_cluster_id=as.integer(prev.cluster_id),
#       #                                 cur_factor1=cur.factor1,
#       #                                 cur_factor2=cur.factor2,
#       #                                 cur_factor3=cur.factor3,
#       #                                 inverse_method=method)
#       ##//set.seed(10)
#       continue=F
#       if(nrow(m$nofactor.cluster.transition)>0){
#         cluster.transition <-sample_transition(m1= m$factor.1.cluster.transition,
#                                                m2= m$factor.2.cluster.transition,
#                                                m3= m$factor.3.cluster.transition,
#                                                m4= m$nofactor.cluster.transition,
#                                                cur_factor1=cur.factor1,
#                                                cur_factor2=cur.factor2,
#                                                cur_factor3=cur.factor3,
#                                                prev_id=as.integer(prev.cluster_id),
#                                                w=c.w,
#                                                calc_cluster_transitions=T)
#         if(cluster.transition[1]!=-1){continue=T}
#       }
#       if(continue){
#         cur.cluster_id <-as.integer(cluster.transition[1])
#         cur.cluster_tbe <- cluster.transition[2]
#       }
#       ##//cluster.transition <- sample.transition(m,1,p.cluster_id=prev.cluster_id ,p.factor1=cur.factor1,p.factor2=cur.factor2,p.factor3=cur.factor3,p.prev_id=prev.cluster_id ,method=method)
#       ##//cur.cluster_tbe <- as.numeric(unlist(cluster.transition[2]))
#       ##//cur.cluster_id <- as.numeric(unlist(cluster.transition[1]))
#     }
#     #sample starting hour and sequence duration
#     #set.seed(6)
#     if(nrow(m$nofactor.start_info)>0 & !prev.cluster_id!=-1 & !is.na(prev.cluster_id)){
#       cur.start.info <- sample_start_info(
#         m1=m$factor.1.start_info,
#         m2=m$factor.2.start_info,
#         m3=m$factor.3.start_info,
#         m4=m$nofactor.start_info,
#         cur_cluster_id=as.integer(prev.cluster_id),
#         cur_factor1=cur.factor1,
#         cur_factor2=cur.factor2,
#         cur_factor3=cur.factor3,
#         w=c.w)
#     } else{break}
#     #cur.start.info <- sample.transition(m,type=2,p.cluster_id=prev.cluster_id,p.factor1=cur.factor1,p.factor2=cur.factor2,p.factor3=cur.factor3,method=method)
#     cur.start.hour <- as.integer(cur.start.info[1])
#     if (is.na(cur.start.hour) || cur.start.hour==-1){
#       cur.start.hour <-  sample(0:23,size=1,replace = T)
#     }
#     start.timestamp<- as.POSIXlt(paste(substr(cur.timestamp,1,11),cur.start.hour,substr(cur.timestamp,14,19),sep=""),origin = "1970-01-01")
#     if (!is.null(end.timestamp)){
#       while( start.timestamp<= end.timestamp){
#         start.timestamp<- start.timestamp + 60*60*24#+1 #add one day
#       }
#     }
#     cur.timestamp<- start.timestamp
#     cur.seq_duration <- as.integer(cur.start.info[2])
#     cur.seq_duration_numeric <- as.numeric(cur.start.info[3])
#     #end.date <- as.POSIXct(paste(substr(cur.timestamp,1,11),"23:59:59",sep=""),origin = "1970-01-01")
#     end.timestamp <-  start.timestamp +cur.seq_duration_numeric
#     sequences <- rbind(sequences,c(prev.cluster_id,cur.seq_duration,as.character(start.timestamp,format="%Y-%m-%d %H:%M:%S"),as.character(end.timestamp,format="%Y-%m-%d %H:%M:%S"),cur.cluster_tbe))
#     #prev.seq_duration_numeric <- cur.seq_duration_numeric
#     #next.date <- as.POSIXct(substr(cur.timestamp,1,11),origin = "1970-01-01")+3600*24+1
#     cur.timestamp <- cur.timestamp + max(cur.seq_duration_numeric,cur.cluster_tbe)
#     #day/4
#     cur.factor1 <- as.integer(ceiling(as.POSIXlt(cur.timestamp,origin = "1970-01-01")$yday/4))
#     #month
#     cur.factor2 <- as.integer(as.POSIXlt(cur.timestamp,origin = "1970-01-01")$mon )
#     #year
#     cur.factor3<- as.integer(as.POSIXlt(cur.timestamp,origin = "1970-01-01")$year )
#     #weekday
#     #cur.weekday<- as.POSIXlt(cur.timestamp,origin = "1970-01-01")$wday
#   }
#   #sequences
#   if(!is.null(sequences)){
#     as.matrix(sequences)
#   } else {NULL}
# }

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
    cur.startTime <- as.POSIXlt(sequences.frame[i,"startTime"],origin = "1970-01-01")
    cur.endTime <- as.POSIXlt(sequences.frame[i,"endTime"],origin = "1970-01-01")
    #tmp<-
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
#     if (i==1){
#       synthetic.data<-tmp
#     } else{
#           synthetic.data<-rbind(synthetic.data,tmp)
#     }
  }


#   new.objects.per.date<- m$new.objects.per.date
#   synthetic.data<-NULL
#   #synthetic.data<- foreach(i = 1:dim(new.objects.per.date)[1],.packages=c("dplyr","psych","foreach"),.combine="rbind") %do% {
#   for(i in 1:dim(new.objects.per.date)[1]){
#     #for (i in 1:dim(new.objects.per.hour)[1]){
#     objects <- new.objects.per.date$objects[i]
#     #obj.bulk <- foreach(cur.obj = 1:objects,.combine="rbind",.export=c("calc.sequences.per.object","sample.transition","generate.transitions","extract.models.sizes","combine")) %do% {
#     obj.bulk <- NULL
#     for(cur.obj in 1:objects){
#       cur.cluster_id <- as.integer(unlist(new.objects.per.date$cluster_id[i]))
#       cur.hour <- new.objects.per.date$hour[i]
#       cur.date <- new.objects.per.date$date[i]
#       cur.weekday <- as.integer(as.POSIXlt(cur.date,origin = "1970-01-01")$wday)
#       cur.week <- as.integer(ceiling(as.POSIXlt(cur.date,origin = "1970-01-01")$yday/4))
#       cur.month <- as.integer(as.POSIXlt(cur.date,origin = "1970-01-01")$mon)
#       cur.year <- as.integer(as.POSIXlt(cur.date,origin = "1970-01-01")$year)
#       cur.seq_duration <- as.integer(unlist(new.objects.per.date$seq_duration[i]))
#       mean.num.seq <- as.numeric(unlist(new.objects.per.date$mean.num.seq[i]))
#       sd.num.seq <- as.numeric(unlist(new.objects.per.date$sd.num.seq[i]    ))
#       cur.num.seq <-
#         max(1,round(rnorm(n=1, mean=mean.num.seq, sd=sd.num.seq),0))
#       cur.timestamp <- as.POSIXlt(cur.date,origin = "1970-01-01")+runif(1,1,3599) #+ as.numeric(cur.hour)*3600 + round(runif(1,0,3600),0)
#       sequences.frame <- calc.sequences.per.object(m,c.timestamp=cur.timestamp,c.num.seq=cur.num.seq,c.cluster_id=cur.cluster_id,c.factor1=cur.week,c.factor2=cur.month,c.factor3=cur.year,method=p.method)
#       if (is.null(sequences.frame)){
#         NULL
#       } else{
#         #names(sequences.frame)<- c("cluster_id","seq_duration","start_timestamp","end_timestamp")
#         #cur.num.seq=5
#         # seq.bulk<-foreach(cur.seq = 1:cur.num.seq,.packages=c("dplyr","psych","Rcpp","sequencesGenerator","psych"),.combine="rbind",.export=c("generate.transitions","extract_models_sizes","sample_transition","sample_starting_state","sample_start_info")) %dopar% {
#         seq.bulk<-NULL
#         seq.bulk<-foreach(cur.seq = 1:cur.num.seq,.packages=c("dplyr","psych","Rcpp","sequencesGenerator","psych"),.export=c("build_vector_of_sizes"),.combine="rbind") %dopar% {
#           #seq.bulk<-foreach(cur.seq = 1:cur.num.seq,.packages=c("dplyr","psych","Rcpp","sequencesGenerator","psych"),.combine="rbind") %dopar% {
#           #seq.bulk<-NULL
#           #for(cur.seq in 1:cur.num.seq){
#           cur.seqid = paste(cur.obj,i,cur.seq,sep="_")
#           cur.obj.id = paste(cur.obj,i,sep="_")
#           #, .export=c("sample.starting.state", "generate.transitions")) %dopar% {
#           #for (cur.seq in 1:cur.num.seq){
#           #cur.seq<- cur.seq+1
#
#           if(nrow(m$factor.1.starting.state)>0){
#             cur.m1<-select_(dplyr::filter(as.data.frame(m$factor.1.starting.state),cluster_id==cur.cluster_id,factor==cur.seq_duration),"sequences")
#             cur.m1<- ifelse(nrow(cur.m1)>0,max(cur.m1),0)
#           }else{cur.m1=0}
#           if(nrow(m$factor.2.starting.state)>0){
#             #cur.m2<- max(m$factor.2.transition[m$factor.2.transition[,"cluster_id"]==cur.cluster_id & m$factor.2.transition[,"factor"]==cur.hour & m$factor.2.transition[,"state_id"]==c.prev_id,"sequences"])
#             cur.m2<-select_(dplyr::filter(as.data.frame(m$factor.2.starting.state),cluster_id==cur.cluster_id,factor==cur.hour),"sequences")
#             cur.m2<- ifelse(nrow(cur.m2)>0,max(cur.m2),0)
#           }else{cur.m2=0}
#
#           if(nrow(m$factor.3.starting.state)>0){
#             #cur.m2<- max(m$factor.3.transition[m$factor.3.transition[,"cluster_id"]==cur.cluster_id & m$factor.3.transition[,"factor"]==cur.weekday & m$factor.3.transition[,"state_id"]==c.prev_id,"sequences"])
#             cur.m3<-select_(dplyr::filter(as.data.frame(m$factor.3.starting.state),cluster_id==cur.cluster_id,factor==cur.weekday),"sequences")
#             cur.m3<- ifelse(nrow(cur.m3)>0,max(cur.m3),0)
#           }else{cur.m3=0}
#
#           c.w<- build_vector_of_sizes(m1_size=cur.m1,
#                                       m2_size=cur.m2,
#                                       m3_size=cur.m3,
#                                       is_inverse=p.method)
#
#           #         if(nrow(m$factor.dur.counts)>1){
#           #           c_m1<- m$factor.dur.counts[,1:3]
#           #         } else{
#           #           c_m1<- matrix(unlist(m$factor.dur.counts[,1:3]),nrow=1)
#           #         }
#           #         if(nrow(m$factor.hour.counts)>1){
#           #           c_m2<- m$factor.hour.counts[,1:3]
#           #         } else{
#           #           c_m2<- matrix(unlist(m$factor.hour.counts[,1:3]),nrow=1)
#           #         }
#           #         if(nrow(m$factor.weekday.counts)>1){
#           #           c_m3<- m$factor.weekday.counts[,1:3]
#           #         } else{
#           #           c_m3<- matrix(unlist(m$factor.weekday.counts[,1:3]),nrow=1)
#           #         }
#           #
#           #         cur.w<- extract_models_sizes( m1=c_m1,
#           #                                   m2=c_m2,
#           #                                   m3=c_m3,
#           #                                   cur_cluster_id=as.integer(sequences.frame[cur.seq,1]),
#           #                                   cur_factor1=as.integer(sequences.frame[cur.seq,2]),
#           #                                   cur_factor2=as.integer(as.POSIXlt(strptime(sequences.frame[cur.seq,3],"%Y-%m-%d %H:%M"),origin = "1970-01-01")$hour),
#           #                                   cur_factor3=as.integer(as.POSIXlt(strptime(sequences.frame[cur.seq,3],"%Y-%m-%d %H:%M"),origin = "1970-01-01")$wday),
#           #                                   inverse_method=p.method)
#           start.state=-1
#           if(nrow(m$nofactor.starting.state)>0){
#             start.state <- sample_starting_state(  m1=as.matrix(m$factor.1.starting.state),
#                                                    m2=as.matrix(m$factor.2.starting.state),
#                                                    m3=as.matrix(m$factor.3.starting.state),
#                                                    m4=as.matrix(m$nofactor.starting.state),
#                                                    cur_cluster_id=as.integer(sequences.frame[cur.seq,1]),
#                                                    cur_factor1=as.integer(sequences.frame[cur.seq,2]),
#                                                    cur_factor2=as.integer(as.POSIXlt(strptime(sequences.frame[cur.seq,3],"%Y-%m-%d %H:%M"),origin = "1970-01-01")$hour),
#                                                    cur_factor3=as.integer(as.POSIXlt(strptime(sequences.frame[cur.seq,3],"%Y-%m-%d %H:%M"),origin = "1970-01-01")$wday),
#                                                    w=c.w)
#           }
#           if (start.state!=-1){
#             cur.state <- as.integer(start.state[1])
#             cur.objects <- as.integer(start.state[2])
#             #cur.objects <- ifelse(is.na(cur.objects),0,cur.objects)
#             cur.freq <-  start.state[3]
#
#             cur.sequence <- generate.transitions(m,
#                                                  p.obj=cur.obj.id,
#                                                  p.seqid = cur.seqid,
#                                                  p.cluster_id=as.numeric(unlist(sequences.frame[cur.seq,1])),
#                                                  p.seq_duration=as.numeric(unlist(sequences.frame[cur.seq,2])),
#                                                  p.timestamp=as.POSIXlt(strptime(sequences.frame[cur.seq,3],"%Y-%m-%d %H:%M:%S"),origin = "1970-01-01"),
#                                                  p.state=cur.state,
#                                                  p.objects=cur.objects,
#                                                  p.freq=cur.freq,
#                                                  end.timestamp=as.POSIXlt(strptime(sequences.frame[cur.seq,4],"%Y-%m-%d %H:%M:%S"),origin = "1970-01-01"),
#                                                  method=p.method)#,
#             #min.objects = p.min.objects,
#             #max.freq=p.max.freq)
#             #cur.sequence
#             as.matrix(cur.sequence)
#           }
#           #####else {NULL}
#
#           #          if(is.null(seq.bulk)){
#           #            seq.bulk<- cur.sequence}
#           #          else {seq.bulk<- rbind(seq.bulk,cur.sequence)}
#           #cur.seqid <- cur.seqid+1
#         }
#         #seq.bulk
#         if(is.null(obj.bulk)){
#           obj.bulk <- seq.bulk
#         }else {
#           obj.bulk<- rbind(obj.bulk,seq.bulk)
#         }
#       }
#       #obj.bulk
#       if(is.null(synthetic.data)){
#         synthetic.data <- obj.bulk
#       }else if(!is.null(obj.bulk)){
#         if (nrow(obj.bulk)>0){
#           synthetic.data<- rbind(synthetic.data,obj.bulk)
#         }
#       }
#     }
#     rm(obj.bulk)
#   }
  foreach::registerDoSEQ()
#   if(!is.null(synthetic.data)){
#     synthetic.data1 <- as.data.frame(synthetic.data,stringsAsFactors = F)
#     colnames(synthetic.data1) <- c("objectid","seq_id","timestamp","state_id")
#     rm(synthetic.data)
#     synthetic.data1
#   } else {synthetic.data}
  #synthetic.data
  synthetic.data<-as.data.frame(synthetic.data,stringsAsFactors = F)
  synthetic.data<-apply(synthetic.data,2,unlist)
  synthetic.data<-synthetic.data[synthetic.data[,"objectid"]!='*',]
  synthetic.data[order(synthetic.data[,"objectid"],synthetic.data[,"timestamp"]),]

  ###dplyr::arrange(as.data.frame(to.valid.matrix(synthetic.data)),objectid,seq_id,timestamp)
}
#}

# output.file.name <- 'try1.csv'
# #registerDoParallel(detectCores())
# synthetic.data <- generate.synthetic.data(m=M,method="inverse.weighted.mean",min.objects=0,max.freq=1)
# write.table(synthetic.data, file = output.file.name, append = F, sep = ",",row.names = F,col.names = T)
# #registerDoSEQ()
#
#


