#include "RDClusterization.hh"
#include <iostream>

using namespace std;

RDClusterization::RDClusterization(const vector<double> *wf , int start, 
				   int n_samp, double thr, int sigma,  
				   int pre_cl, int post_cl) 
{
   vector<int> x;
   int step = 100;  // hardcoded FIXME   

   for (int i = start; i < n_samp - sigma; i++) {
      double mw = wf->at(i);
      for (int k = -sigma; k < sigma; k = k + step) { 
         if (wf->at(i + k) < mw) mw = wf->at(i + k);
      } 
       
      if (wf->at(i) <= mw && mw < thr && (wf->at(i) - wf->at(i - pre_cl)) < thr) {
         x.push_back(i);
         i += sigma;  
      }
   }
    
   if (x.size()) {
      for (size_t p = 0; p < x.size(); p++) {
         //stime_.push_back(x[p]); 
         str_.push_back(x[p] - pre_cl);

         if (p == x.size() - 1) {
            if (x[p] + post_cl >= n_samp)  {
               stp_.push_back(n_samp - 1);
               rp_.push_back(0);
            } 
            else { 
               stp_.push_back(x[p] + post_cl);
               rp_.push_back(1); 
            }
         }
         else {
            if (x[p + 1] - x[p] - pre_cl < post_cl) {
               stp_.push_back(x[p + 1] - pre_cl);
               rp_.push_back(0); 
            }
            else {
               stp_.push_back(x[p] + post_cl);
               rp_.push_back(1); 
            } 
         }
      }
    
      // fixed thresghold startime on the smoothed wf
    
      for (size_t p = 0; p < x.size(); p++) {
         for (int k = str_[p]; k < stp_[p]; k++) {
            if (wf->at(k) < thr) { 
               stime_.push_back(k);
               break;
            }    
         }   
      }

   } 
}
