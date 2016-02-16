sequencesGenerator <- function(f,objectid.pos=1,seq_id.pos=0,timestamp.pos=2,stateid.pos=3,secs.between.sequences=NA,timestamp.format='%Y-%m-%dT%H:%M:%S',...) UseMethod("sequencesGenerator")

# sequencesGenerator.default <- function(f,objectid.pos=1,seq_id.pos=0,timestamp.pos=2,stateid.pos=3,secs.between.sequences=NA,timestamp.format='%Y-%m-%dT%H:%M:%S',...) {
#   seqGenerator<-NULL
#   df<- init.transactions.dataset(f,objectid.pos=objectid.pos,timestamp.pos=timestamp.pos,stateid.pos=stateid.pos,seq_id.pos=seq_id.pos,secs.between.sequences=secs.between.sequences,timestamp.format=timestamp.format)
#   seqGenerator$src <- df
#   seqGenerator$model <- NA
#   class(seqGenerator) <- "sequencesGenerator"
#   seqGenerator
# }

print.sequencesGenerator <- function (x, ...) {
  x$src
  x$states_dictionary
}
add_states <- function (x, ...) {
  x$src$state <- with(x$states_dictionary, value[match(x$src$state_id,key)])
  x
}
#library(compiler)
#enableJIT(1)

#Load sequential dataset and standartize its structure.
#input.file.name - path to the source dataset
#header - true if the source dataset contains headers, otherwise false
#sep - seperator in source file(dataset)
#objectid.pos - index of object-id columns in the source dataset (if more than one column exists they are concatenated)
#timestamp.pos - index of timestamp column in the source dataset (a single position only)
#stateid.pos - index of state-id columns in the source dataset (if more than one column exists they are concatenated)
#seq_id.pos - index of a unique sequence-id column in the source dataset (a single position only). 0 if no sequence id column exists in the source dataset
#nrows = rows limit (for sampling the source dataset)
#
#Returns a processed structure of the data with the following fields:
#objectid  state_id  timestamp	hour	date	seq_id	next_objectid	next_seqid	next_state_id	tbe


#loading data from a text file into a <objectid,timestamp,stateid> representation
#init.transactions.dataset <- function(data,objectid.pos,timestamp.pos,stateid.pos,seq_id.pos,timestamp.format,secs.between.sequences=NA)
sequencesGenerator.default <- function(data,objectid.pos=1,seq_id.pos=0,timestamp.pos=2,stateid.pos=3,secs.between.sequences=NA,timestamp.format='%Y-%m-%d %H:%M:%S',...) {
  seqGenerator<-NULL

  d <- as.data.frame(data[(!is.na(data[,objectid.pos]) && !is.na(data[,timestamp.pos]) &&  !is.na(data[,stateid.pos])),])
  m<- d[,0]
  if (length(objectid.pos)>1) {
    m$objectid <- unlist(lapply(as.data.frame(do.call(paste0, d[,objectid.pos])), as.numeric))
  }  else {
    m$objectid <- unlist(lapply(as.data.frame(d[,c(objectid.pos)]), as.numeric))
  }
  #concatenate state fields if more then one exist
  if (length(stateid.pos)>1){
    values<- as.data.frame(do.call(paste0, d[,stateid.pos]))
    m$state_id <- unlist(lapply(values,as.numeric))
    #m$state_id <- unlist(lapply(as.data.frame(do.call(paste0, d[,stateid.pos])),as.numeric))
  }else {
    values<- d[,stateid.pos]
    m$state_id <- unlist(lapply(values,as.numeric))
    #m$state_id <- unlist(lapply(d[,stateid.pos],as.numeric))
    ####to preserve a list of state ids (factor ids)
  }
  seqGenerator$states_dictionary<- unique(cbind(m$state_id,values))
  names(seqGenerator$states_dictionary)<-c("key","value")
  m$timestamp<- as.POSIXct(unlist(d[,timestamp.pos]),origin = "1970-01-01",format=timestamp.format)

  m <- m[order(m$objectid,m$timestamp),]

  if (seq_id.pos>0) {
    m$seq_id<- d[,seq_id.pos]
  } else { #seperate sequences by secs.between.sequences
    m$next_objectid <- c(m$objectid[2:dim(m)[1]],NA)
    m$tbe <- c(m$timestamp[2:dim(m)[1]],NA)- m$timestamp
    m$tbe <- ifelse(m$next_objectid!=m$objectid, NA, m$tbe)

    seq.ids<- m$tbe>secs.between.sequences
    seq.ids<- as.numeric(seq.ids + m$next_objectid!=m$objectid)
    seq.ids[is.na(seq.ids)]<-1
    seq.ids<- c(0,seq.ids[1:length(seq.ids)-1])
    m$seq_id<-as.numeric(factor(cumsum(seq.ids)))
  }
  seqGenerator$src <- m[!is.na(m$seq_id),c("objectid","seq_id","timestamp","state_id")]
  class(seqGenerator) <- "sequencesGenerator"
  seqGenerator
}
