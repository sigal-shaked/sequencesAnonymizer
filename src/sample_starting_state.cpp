#include <Rcpp.h>
using namespace Rcpp;



////#define COMPILING_SEQUENCESGENERATOR
#include <iostream>
#include <map>
#include <sstream>
#include <string>
// Returns a starting state and the minimal number of objects that started with this state
// according to models calculated by factor1,factor2 and factor3, or with no factor at all when relevant factors contain no matchings.
// Input parameters: three models m1-m3 by possible factors, model m4 with no factor, current clusterid and factor values, and weights w of the different factors


// [[Rcpp::export]]
NumericVector sample_starting_state(NumericMatrix  m1, NumericMatrix m2, NumericMatrix m3,NumericMatrix m4,NumericVector w=NumericVector::create(0.33,0.33,0.33)){
    std::multimap<int,NumericVector > states;
    int pstate_id=0,  ptotal=1;

    int nrow = m1.nrow(),key=0,r=0;
    double val=0,total=0;
    if((nrow>0) & (w[0]>0)){
      for (int i = 0; i < nrow; i++) {
      //  if (m1(i,pcluster)==cur_cluster_id && m1(i,pfactor)==cur_factor1 && !NumericVector::is_na(m1(i,pstate_id))&& !NumericVector::is_na(m1(i,ptotal))) {
              states.insert(std::pair<int,NumericVector > (m1(i,pstate_id),(NumericVector::create(w[0],m1(i,ptotal)))));
        }
      //}
    }
    nrow = m2.nrow();
    if((nrow>0) & (w[1]>0)){
      for (int i = 0; i < nrow; i++) {
    //    if (m2(i,pcluster)==cur_cluster_id && m2(i,pfactor)==cur_factor2 && !NumericVector::is_na(m2(i,pstate_id))&& !NumericVector::is_na(m2(i,ptotal))) {
              states.insert(std::pair<int,NumericVector > (m2(i,pstate_id),(NumericVector::create(w[1],m2(i,ptotal)))));
     //   }
      }
    }
    nrow = m3.nrow();
    if((nrow>0) & (w[2]>0)){
      for (int i = 0; i < nrow; i++) {
       // if (m3(i,pcluster)==cur_cluster_id && m3(i,pfactor)==cur_factor3 && !NumericVector::is_na(m3(i,pstate_id))&& !NumericVector::is_na(m3(i,ptotal))) {
              states.insert(std::pair<int,NumericVector > (m3(i,pstate_id),(NumericVector::create(w[2],m3(i,ptotal)))));
       // }
      }
    }

    if (states.empty()){
      nrow = m4.nrow();
      if(nrow>0){
        //pcluster=0,  pstate_id=1,  pobjects=3, ptotal=4;
        for (int i = 0; i < nrow; i++) {
        //  if (m4(i,pcluster)==cur_cluster_id && !NumericVector::is_na(m4(i,pstate_id))&& !NumericVector::is_na(m4(i,ptotal))) {
                states.insert(std::pair<int,NumericVector > (m4(i,pstate_id),(NumericVector::create(1,m4(i,ptotal)))));
        //  }
        }
      }
    }


    NumericMatrix out(states.size(),3);

    if (!states.empty()){
      std::multimap<int,NumericVector >::iterator it = states.begin();
      //for (std::multimap<int,NumericVector >::const_iterator it = states.begin(); it != states.end(); ++it) {
      while(it != states.end()){
        key = (*it).first;
        val=((*it).second)[0]*((*it).second)[1];
        //objects= ((*it).second)[0]*((*it).second)[2];
        //cnt_objects= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
        //cnt=((*it).second)[0];
        ++it;
        while (it != states.end() && key == it->first) {
          val += ((*it).second)[0]*((*it).second)[1];
          //objects +=((*it).second)[0]*((*it).second)[2];
          //cnt += ((*it).second)[0];
          //cnt_objects+= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
          ++it;
        }
        out(r,0)=key;
        out(r,1)=val;//(cnt>0? val/cnt:0);
        //objects = objects/cnt;
        //out(r,2)=(cnt_objects>0? objects/cnt_objects:0);
        total+=out(r,1);
        r++;
    }
  }
    IntegerVector keys(r);
    NumericVector probs(r);
    std::map<int,double > objs;
    std::map<int,double > freqs;

    if(total>0){
      for(int i = 0; i < r; i++){
          keys[i]=out(i,0);
          probs[i]=out(i,1)/total;
          //objs[i]=out(i,2);
          //objs.insert(std::pair<int,double > (out(i,0),out(i,2)));
          freqs.insert(std::pair<int,double > (out(i,0),probs[i]));
      }
    }
   double rnd = R::runif(0,1);
   double sum=0;
   int z=0;
   while((z<r) && (rnd>sum)){
     sum+=probs[z];
     z++;
   }
   //cid = RcppArmadillo::sample(ids, 1, TRUE, freqs)[0] ;
   key = keys[z-1];
   //Rcpp::Rcout <<keys[z]<<" "<< freqs.find(keys[z])->second<< std::endl ;

   NumericVector ret(1);
   ret[0]= key;//RcppArmadillo::sample(keys, 1, TRUE, probs)[0] ;
   //ret[1]= objs.find(key)->second;
   //ret[2]= freqs.find(key)->second;
   //ret[1] = (!NumericVector::is_na(ret[1])? ret[1]:0);
   //ret[2] = (!NumericVector::is_na(ret[2])? ret[2]:0);
return ret;
}



