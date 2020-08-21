// Peak finder module
#include <algorithm>

class peakfinder {
   int pb_;
   double *x_, *y_;
   int n_;
   double wmin_, wmax_;
   int imin_, imax_;
   int st_;
   int Nsample;
   double *WF;

   public:
      peakfinder(double *wf, int min, int max, int sigma, double threshold, int peak_buffer);
      double *x_peak() {return x_;}
      double *y_peak() {return y_;}
      double n_peak()  {return n_;}  
      int x_min()      {return imin_;}     
      double y_min()   {return wmin_;}     
      int x_max()      {return imax_;}    
      double y_max()   {return wmax_;}    
      int start_t()    {return st_;}
      void dump();
};

peakfinder::peakfinder(double *wf, int min, int max, int sigma, double threshold, int peak_buffer) {
    int step = 10;
    double cdf = 0.5;

    pb_ = peak_buffer;
    x_ = new double[pb_];
    y_ = new double[pb_];

    wmin_ = wf[min];
    wmax_ = wf[min];
    
    imin_ = min;
    imax_ = min;
    
    Nsample = max;
    
    WF=wf;
    
    //this->dump();
    
    // min-max
    for(int i = min; i < max; i++) {
        
        //cout <<" wf[" << i << "]=" << wf[i] << endl;
        
       if (wf[i] < wmin_) {
          wmin_ = wf[i];
          imin_ = i;
       }
          
       if (wf[i] > wmax_) {
          wmax_ = wf[i];        
          imax_ = i;
       }

    }    
 
    // start time   
    st_ = min;   
    for(int i = min; i < max; i++) {
       if (wf[i] < wmax_ - (wmax_ - wmin_)*cdf) {
          st_ = i; 
          break;
       }
    }     

    // peak_finder
    for (int i = 0; i < pb_; i++) {
       x_[i] = 0;
       y_[i] = 0;   
    }
     
    int np = 0;
     
    for (int i = min + sigma; i < max - sigma; i++) {
       if (np == pb_) break;

       int mw = wf[i];
       for (int k = - sigma; k < sigma; k = k + step) { 
          if (wf[i + k] < mw) mw = wf[i + k];
       } 
       
       if ( wf[i] <= mw && wf[i]  < wmax_ - threshold*(wmax_ - wmin_) ) {
          x_[np] = i;  
          y_[np] = wf[i];
          np++;
          i += sigma;
           
       }
        //cout <<" x_[" << i << "]=" << i <<" y_[" << i << "]=" << wf[i] << endl;
    }
    
    if (np == 0) {
      x_[0] = imin_;
      y_[0] = wmin_;
    }

    n_ = np; 
}

/*
void peakfinder::dump() {
    
    for (int t=0; t<Nsample; t++)
        
        cout <<" WF[" << t << "]=" << WF[t] << endl;
    
}
*/
