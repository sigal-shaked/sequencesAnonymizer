#include <Rcpp.h>
using namespace Rcpp;



////#define COMPILING_SEQUENCESGENERATOR
#include <iostream>
#include <map>
#include <sstream>
#include <string>



// [[Rcpp::export]]
CharacterVector sample_start_info(NumericMatrix  m1, NumericMatrix m2, NumericMatrix m3,NumericMatrix m4,NumericVector w=NumericVector::create(0.33,0.33,0.33)){
    std::multimap<std::string,NumericVector > states;
    int  phour=0, pduration=1, ptotal=2, pdur_mean=3, pdur_sd=4;

    //int size =m2.nrow()+m2.nrow()+m3.nrow();
    //size = (size > m4.nrow()? size:m4.nrow());
    //SEXP x;
    int nrow = m1.nrow(),r=0,cid=0;
    std::string key;
    double cnt_duration_mean=0,cnt_duration_sd=0,val=0,total=0,duration_mean=0,duration_sd=0,cduration_mean=0,cduration_sd=0,cduration=0;

    if((nrow>0) & (w[0]>0)){
      for (int i = 0; i < nrow; i++) {
       // if (m1(i,pcluster)==cur_cluster_id && m1(i,pfactor)==cur_factor1 && !NumericVector::is_na(m1(i,phour))&& !NumericVector::is_na(m1(i,pduration))&&!NumericVector::is_na(m1(i,ptotal))) {
              std::stringstream ckey;
              ckey<< m1(i,phour)<<"_"<<m1(i,pduration) ;
              states.insert(std::pair<std::string ,NumericVector > (ckey.str(),(NumericVector::create(w[0],m1(i,ptotal),m1(i,pdur_mean),m1(i,pdur_sd)))));
              //Rcpp::Rcout <<" i:"<<i<<" ckey:"<<ckey.str()<< std::endl ;
      //  }
      }
    }
    nrow = m2.nrow();
    if((nrow>0) & (w[1]>0)){
      for (int i = 0; i < nrow; i++) {
        //if (m2(i,pcluster)==cur_cluster_id && m2(i,pfactor)==cur_factor2 && !NumericVector::is_na(m2(i,phour))&& !NumericVector::is_na(m2(i,pduration))&&!NumericVector::is_na(m2(i,ptotal))) {
              std::stringstream ckey ;
              ckey<< m2(i,phour)<<"_"<<m2(i,pduration) ;
              states.insert(std::pair<std::string ,NumericVector > (ckey.str(),(NumericVector::create(w[1],m2(i,ptotal),m2(i,pdur_mean),m2(i,pdur_sd)))));
        //}
      }
    }
    nrow = m3.nrow();
    if((nrow>0) & (w[2]>0)){
      for (int i = 0; i < nrow; i++) {
       // if (m3(i,pcluster)==cur_cluster_id && m3(i,pfactor)==cur_factor3 && !NumericVector::is_na(m3(i,phour))&& !NumericVector::is_na(m3(i,pduration))&&!NumericVector::is_na(m3(i,ptotal))) {
              std::stringstream ckey ;
              ckey<< m3(i,phour)<<"_"<<m3(i,pduration) ;
              states.insert(std::pair<std::string ,NumericVector > (ckey.str(),(NumericVector::create(w[2],m3(i,ptotal),m3(i,pdur_mean),m3(i,pdur_sd)))));
       // }
      }
    }
    if (states.empty()){
      nrow = m4.nrow();
      if(nrow>0){
        //pcluster=0,  phour=1, pduration=2, pobjects=4, ptotal=5, pdur_mean=6, pdur_sd=7;
        for (int i = 0; i < nrow; i++) {
        //  if (m4(i,pcluster)==cur_cluster_id &&  !NumericVector::is_na(m4(i,phour))&& !NumericVector::is_na(m4(i,pduration))&&!NumericVector::is_na(m4(i,ptotal))) {
                std::stringstream ckey ;
                ckey<< m4(i,phour)<<"_"<<m4(i,pduration) ;
                states.insert(std::pair<std::string ,NumericVector > (ckey.str(),(NumericVector::create(1,m4(i,ptotal),m4(i,pdur_mean),m4(i,pdur_sd)))));
        //  }
        }
      }
    }

    CharacterVector out_keys(states.size());
    NumericMatrix out(states.size(),4);

    if (!states.empty()){
      std::multimap<std::string,NumericVector >::iterator it = states.begin();
      //for (std::multimap<int,NumericVector >::const_iterator it = states.begin(); it != states.end(); ++it) {
      while(it != states.end()){
        key = (*it).first;
        val=((*it).second)[0]*((*it).second)[1];
  //      objects= ((*it).second)[0]*((*it).second)[2];
  //      cnt_objects= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
        duration_mean= ((*it).second)[0]*((*it).second)[2];
        cnt_duration_mean= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
        duration_sd= ((*it).second)[0]*((*it).second)[3];
        cnt_duration_sd= (!NumericVector::is_na(((*it).second)[3])? ((*it).second)[0]:0);

        //cnt=((*it).second)[0];
        ++it;
        while (it != states.end() && key == it->first) {
          val += ((*it).second)[0]*((*it).second)[1];
          //objects +=((*it).second)[0]*((*it).second)[2];
          //cnt_objects += (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
          duration_mean= ((*it).second)[0]*((*it).second)[2];
          cnt_duration_mean= (!NumericVector::is_na(((*it).second)[2])? ((*it).second)[0]:0);
          duration_sd= ((*it).second)[0]*((*it).second)[3];
          cnt_duration_sd= (!NumericVector::is_na(((*it).second)[3])? ((*it).second)[0]:0);

          //cnt += ((*it).second)[0];
          ++it;
        }
        out_keys[r]=key;
        out(r,0)=val;//(cnt>0? val/cnt:0);
        //out(r,1)=(cnt_objects>0? objects/cnt_objects:0);
        out(r,1)=(cnt_duration_mean>0? duration_mean/cnt_duration_mean:0);
        out(r,2)=(cnt_duration_sd>0? duration_sd/cnt_duration_sd:0);
        total+=out(r,0);
        r++;
    }
  }
    CharacterVector keys(r);
    IntegerVector ids(r);
    NumericVector freqs(r);
    //std::map<std::string,double > objs;
    std::map<std::string,double > duration_means;
    std::map<std::string,double > duration_sds;

    if(total>0){
      for(int i = 0; i < r; i++){
          ids[i]=i;
          key = Rcpp::as<std::string>(out_keys[i]);
          keys[i]=key;
          freqs[i]=out(i,0)/total;
          //objs[i]=out(i,2);
          //objs[key]=out(i,1);
          duration_means[key]=out(i,1);
          duration_sds[key]=out(i,2);
          //objs.insert(std::pair<std::string,double > (out_keys[i],out(i,1)));
          //duration_means.insert(std::pair<std::string,double > (key,out(i,2)));
          //duration_sds.insert(std::pair<std::string,double > (key,out(i,3)));
      }
    }
   CharacterVector ret(3);
   double rnd = R::runif(0,1);
   double sum=0;
   int z=0;
   while((z<r) && (rnd>sum)){
     sum+=freqs[z];
     z++;
   }
   cid=z-1;
   //cid = RcppArmadillo::sample(ids, 1, TRUE, freqs)[0] ;
   key = Rcpp::as<std::string>(keys[cid]);
   std::size_t pos = (key.length()>0? key.find("_"):-1);
   ret[0] =(((key.length()>0) & (pos>=0)? key.substr(0,pos):"-1"));
   ret[1] =(((key.length()>0) & (pos>=0)? key.substr(pos+1,key.length()):"-1"));
   cduration_mean= duration_means.find(key)->second;
   cduration_sd= duration_sds.find(key)->second;
   cduration= R::rnorm(cduration_mean,cduration_sd);
   ret[2]=(((!NumericVector::is_na(cduration)) & (cduration>0))? cduration:0);
   //ret[3]= objs.find(key)->second;
   //ret[3] = (!NumericVector::is_na(objs.find(key)->second)? objs.find(key)->second:0);
   return ret;
}
