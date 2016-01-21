#include <Rcpp.h>
using namespace Rcpp;



////#define COMPILING_SEQUENCESGENERATOR
#include <iostream>
#include <map>
#include <sstream>
#include <string>


// Returns the next state, the minimal number of objects that transit from prev_id to this state and the estimated tbe (time between prev_id and next state)
// according to models calculated by factor1,factor2 and factor3, or with no factor at all when relevant factors contain no matchings.
// Input parameters: three models m2-m3 by possible factors, model m4 with no factor, current clusterid, previous state_id (prev_id) and factor values, and weights w of the different factors


// [[Rcpp::export]]
NumericVector sample_transition(NumericMatrix  m1, NumericMatrix m2, NumericMatrix m3,NumericMatrix m4,NumericVector w=NumericVector::create(0.33,0.33,0.33)){

  std::multimap<int,NumericVector > states;
  int  pnext_id, ptotal,  ptbe_mean, ptbe_sd;
  pnext_id=0; ptotal=1; ptbe_mean=2; ptbe_sd=3;
/*  if (calc_cluster_transitions){
    pcluster_id=-1;pprev_id=0; pfactor=1; pnext_id=2; ptotal=5; pobjects=4;  ptbe_mean=6; ptbe_sd=7;
  }
  else {
    pcluster_id=0;pprev_id=2; pfactor=1; pnext_id=3; ptotal=6; pobjects=5;  ptbe_mean=7; ptbe_sd=8;
  }*/
  int nrow = m1.nrow(),key=0,r=0;
  double /*cnt=0,*/val=0,total=0,tbe_mean=0,tbe_sd=0,ctbe_mean=0,ctbe_sd=0,ctbe=0,cnt_tbe_mean=0,cnt_tbe_sd=0;

  if((nrow>0) & (w[0]>0)){
    for (int i = 0; i < nrow; i++) {
   //   if ((calc_cluster_transitions | (m1(i,pcluster_id)==cur_cluster_id)) && m1(i,pfactor)==cur_factor1&& m1(i,pprev_id)==prev_id && !NumericVector::is_na(m1(i,pnext_id))&& !NumericVector::is_na(m1(i,ptotal))) {
        states.insert(std::pair<int,NumericVector > (m1(i,pnext_id),(NumericVector::create(w[0],m1(i,ptotal),m1(i,ptbe_mean),m1(i,ptbe_sd)))));
        //Rcpp::Rcout <<" i:"<<i<<" pnext_id:"<<pnext_id<< std::endl ;
   //   }
    }
  }

  nrow = m2.nrow();
  if((nrow>0) & (w[1]>0)){
    for (int i = 0; i < nrow; i++) {
     // if ((calc_cluster_transitions | (m2(i,pcluster_id)==cur_cluster_id)) && m2(i,pfactor)==cur_factor2&& m2(i,pprev_id)==prev_id && !NumericVector::is_na(m2(i,pnext_id))&& !NumericVector::is_na(m2(i,ptotal))) {
        states.insert(std::pair<int,NumericVector > (m2(i,pnext_id),(NumericVector::create(w[1],m2(i,ptotal),m2(i,ptbe_mean),m2(i,ptbe_sd)))));
      //}
    }
  }
  nrow = m3.nrow();
  if((nrow>0) & (w[2]>0)){
    for (int i = 0; i < nrow; i++) {
      //if ((calc_cluster_transitions | (m3(i,pcluster_id)==cur_cluster_id)) && m3(i,pfactor)==cur_factor3&& m3(i,pprev_id)==prev_id && !NumericVector::is_na(m3(i,pnext_id))&& !NumericVector::is_na(m3(i,ptotal))) {
        states.insert(std::pair<int,NumericVector > (m3(i,pnext_id),(NumericVector::create(w[2],m3(i,ptotal),m3(i,ptbe_mean),m3(i,ptbe_sd)))));
      //}
    }
  }

  nrow = m4.nrow();
  if (states.empty()){
    if(nrow>0){
     /* if (calc_cluster_transitions){
        pcluster_id=-1;pprev_id=0; pfactor=-1; pnext_id=1; ptotal=4; pobjects=3;  ptbe_mean=5; ptbe_sd=6;
      }
      else {
        pcluster_id=0;pprev_id=1; pfactor=-1; pnext_id=2; ptotal=5; pobjects=4;  ptbe_mean=6; ptbe_sd=7;
      }*/
      for (int i = 0; i < nrow; i++) {
        //if ((calc_cluster_transitions | (m4(i,pcluster_id)==cur_cluster_id)) && m4(i,pprev_id)==prev_id && !NumericVector::is_na(m4(i,pnext_id))&& !NumericVector::is_na(m4(i,ptotal))) {
          states.insert(std::pair<int,NumericVector > (m4(i,pnext_id),(NumericVector::create(1,m4(i,ptotal),m4(i,ptbe_mean),m4(i,ptbe_sd)))));
       // }
      }
    }
  }


  NumericMatrix out(states.size(),5);

  if (!states.empty()){
    std::multimap<int,NumericVector >::iterator it = states.begin();
    //for (std::multimap<int,NumericVector >::const_iterator it = states.begin(); it != states.end(); ++it) {
    while(it != states.end()){
      key = (*it).first;
      val=((*it).second)[0]*((*it).second)[1];
     // objects= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]*((*it).second)[2]:0);
      //cnt_objects= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
      tbe_mean=(!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]*((*it).second)[2]:0);
      cnt_tbe_mean= (!NumericVector::is_na(((*it).second)[3])? ((*it).second)[0]:0);
      tbe_sd=  (!NumericVector::is_na(((*it).second)[3])? ((*it).second)[0]*((*it).second)[3]:0);
      cnt_tbe_sd= (!NumericVector::is_na(((*it).second)[3])? ((*it).second)[0]:0);

      //cnt=((*it).second)[0];
      ++it;
      while (it != states.end() && key == it->first) {
        val += ((*it).second)[0]*((*it).second)[1];
        //objects +=(!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]*((*it).second)[2]:0);
        //cnt_objects+= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
        tbe_mean+= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]*((*it).second)[2]:0);
        cnt_tbe_mean+= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
        tbe_sd+= (!NumericVector::is_na(((*it).second)[3])? ((*it).second)[0]*((*it).second)[3]:0);
        cnt_tbe_sd+= (!NumericVector::is_na(((*it).second)[3])? ((*it).second)[0]:0);
        //cnt += ((*it).second)[0];
        ++it;
      }
      out(r,0)=key;
      out(r,1)=val;//(cnt>0? val/cnt:0);
      //out(r,2)=(cnt_objects>0? objects/cnt_objects:0);
      out(r,2)=(cnt_tbe_mean>0? tbe_mean/cnt_tbe_mean:0);
      out(r,3)=(cnt_tbe_sd>0? tbe_sd/cnt_tbe_sd:0);
      total+=out(r,1);
      //Rcpp::Rcout <<key<<" "<<out(r,1)<< std::endl ;

      r++;
    }
  }
  IntegerVector keys(r);
  NumericVector probs(r);
  //std::map<int,double > objs;
  std::map<int,double > freqs;
  std::map<int,double > tbe_means;
  std::map<int,double > tbe_sds;


  if(total>0){
    for(int i = 0; i < r; i++){
      keys[i]=out(i,0);
      probs[i]=out(i,1)/total;
      //Rcpp::Rcout <<keys[i]<<" "<<out(i,1)<<"/"<<total<<"="<<probs[i]<< std::endl ;
      //objs[i]=out(i,2);
      //objs.insert(std::pair<int,double > (out(i,0),out(i,2)));
      freqs.insert(std::pair<int,double > (out(i,0),probs[i]));
      tbe_means.insert(std::pair<int,double > (out(i,0),out(i,2)));
      tbe_sds.insert(std::pair<int,double > (out(i,0),out(i,3)));
    }

    double rnd = R::runif(0,1);
    double sum=0;
    int z=0;

    while((z<r) && (rnd>sum)){
      sum+=probs[z];
      z++;
    }
    //Rcpp::Rcout <<" r:"<<r<<" probs.length():"<<probs.length()<<" rnd:"<<rnd<<" sum:"<<sum<<" key:"<<keys[z-1]<<"z-1:"<<z-1<< std::endl ;

    //cid = RcppArmadillo::sample(ids, 1, TRUE, freqs)[0] ;
    key = keys[z-1];
    NumericVector ret(2);
    //ret[0]= RcppArmadillo::sample(keys, 1, TRUE, probs)[0] ;
    ret[0]= key;
    ctbe_mean= tbe_means.find(key)->second;
    ctbe_sd= tbe_sds.find(key)->second;
    ctbe= R::rnorm(ctbe_mean,ctbe_sd);
    ret[1]=(ctbe>0? ctbe:0);
    //ret[2]= objs.find(key)->second;
    //ret[2]= freqs.find(key)->second;
    ret[0] = (!NumericVector::is_na(ret[0])? ret[0]:-1);
    ret[1] = (!NumericVector::is_na(ret[1])? ret[1]:0);
    //ret[2] = (!NumericVector::is_na(ret[2])? ret[2]:0);
  //  ret[3] = (!NumericVector::is_na(ret[3])? ret[3]:0);
    return ret;
  }
  else{
    NumericVector ret(3);
    ret[0]=-1;
    ret[1] = 0;
    ret[2] = 0;
    //ret[3] = 0;
    return ret;
  }
}

//Returns <next_id,tbe,#objects,freq>
