require(TraMineR)
#library(dplyr)

calc.seq.obj <- function(d)
{
  #creating TSE object and converting it into STS object
  #d.seq_id <- as.numeric(factor(d$seq_id))
  d.seq_id <- d$seq_id

  d.seqe <- seqecreate(id = as.numeric(factor(d.seq_id)), timestamp = (as.numeric(as.POSIXct(d$timestamp)) - as.numeric(as.POSIXct('1970-01-01'))),event = d$state_id)
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
  d.seq <- seqdef(tmp, sep=",")
  rm(tmp)

  #seq.idx <- group_by(d,seq_id)
  #seq.idx <-arrange(seq.idx,timestamp)
  #seq.idx <- filter(seq.idx, timestamp == first(timestamp))$seq_id
  #list(seq.idx,d.seq)
  list(seqeid(d.seqe),d.seq)
}



calc.pairwise.lcs <- function(d.seq)
{
  #calculate matrix of distances according to Least Common Subsequence measure
  d.lcs <- seqdist(d.seq, method = "LCS",with.missing=TRUE, norm=TRUE)
  d.lcs
}
