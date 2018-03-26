//
// Integer representation of floating point arithmetic suitable for FPGA designs
// 
// Author: Yuri Gershtein 
// Date:   March 2018
//
// Functionality:
//
//  *note* all integers are assumed to be signed
//
//  all variables have units, stored in a map <string,int>, with string a unit (i.e. "phi") and int the power
//                   "2" is always present in the map, and it's int pair is referred to as 'shift'
//                   units are properly combined / propagated through calculations
//                   adding/subtracting variables with different units throws an exception
//                   adding/subtracting variables with different shifts is allowed and is handled correctly
//
// calculate() method re-calculates the variable double and int values based on its operands
//                   returns false in case of overflows and/or mismatches between double and int calculations.
//
// the maximum and minimum values that the variable assumes are stored and updated each time calculate() is called
// if IMATH_ROOT is defined, all values are also stored in a histogram
//
// var_def     (string name, string units, double fmax, double K):
//                   define variable with bit value fval = K*ival, and maximum absolute value fmax.
//                   calculates nbins on its own
//                   one can assign value to it using set_ methods 
//
// var_add     (string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
// var_subtract(string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
//                   add/subtract variables. Bit length increases by 1, but capped at nmax
//                   if range>0 specified, bit length is decreased to drop unnecessary high bits
//
// var_mult    (string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
//                   multiplication. Bit length is a sum of the lengths of the operads, but capped at nmax
//                   if range>0 specified, bit length is decreased to drop unnecessary high bits or post-shift is reduced
//
// var_timesC  (string name, var_base *p1, double cF, int ps = 17):
//                   multiplication by a constant. Bit length stays the same
//                   ps defines number of bits used to represent the constant
//
// var_shift  (string name, var_base *p1, int shift):
//                   shifts the variable right by shift (equivalent to multiplication by pow(2, -shift));
//                   Units stay the same, nbits are adjusted.
//
// var_neg    (string name, var_base *p1):
//                   multiplies the variable by -1
//
// var_inv     (string name, var_base *p1, double offset, int nbits, int n, unsigned int shift, mode m, int nbaddr=-1):
//                   LUT-based inversion, f = 1./(offset + f1) and  i = 2^n / (offsetI + i1)
//                   nbits is the width of the LUT (signed)
//                   m is from enum mode {pos, neg, both} and refers to possible sign values of f
//                            for pos and neg, the most significant bit of p1 (i.e. the sign bit) is ignored
//                   shift is a shift applied in i1<->address conversions (used to reduce size of LUT)
//                   nbaddr: if not specified, it is taken to be equal to p1->get_nbits()
//                           
//
// var_nounits (string name, var_base *p1, int ps = 17):
//                   convert a number with units to a number - needed for trig function expansion (i.e. 1 - 0.5*phi^2)
//                   ps is a number of bits to represent the unit conversion constant
//
// var_adjustK (string name, var_base *p1, double Knew, double epsilon = 1e-5, bool do_assert = false, int nbits = -1)
//                   adjust variable shift so the K is as close to Knew as possible (needed for bit length adjustments) 
//                   if do_assert is true, throw an exeption if Knew/Kold is not a power of two
//                   epsilon is a comparison precision, nbits forces the bit length (possibly discarding MSBs)
//
//
// bool calculate(int debug_level) runs through the entire formula tree recalculating both ineteger and floating point values
//                     returns true if everything is OK, false if obvious problems with the calculation exist, i.e
//                                  -  integer value does not fit into the alotted number of bins
//                                  -  integer value is more then 10% or more then 2 away from fval_/K_ 
//                     debug_level:  0 - no warnings
//                                   1 - limited warning
//                                   2 - as 1, but also include explicit warnings when LUT was used out of its range
//                                   3 - maximum complaints level
//
//
#ifndef IMATH_H
#define IMATH_H

//use root if uncommented
#define IMATH_ROOT

#include <limits>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <assert.h>

#ifdef IMATH_ROOT
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#endif

class var_base {
  
 public:

  var_base(std::string name, var_base *p1, var_base *p2, int l){
    p1_ = p1;
    p2_ = p2;
    name_ = name;
    latency_ = l;
    int step1 = (p1)? p1->get_step()+p1->get_latency() : 0;
    int step2 = (p2)? p2->get_step()+p2->get_latency() : 0;
    step_ = step1>step2? step1 : step2;
    
    minval_ = std::numeric_limits<double>::max();
    maxval_ = -std::numeric_limits<double>::max();
    readytoprint_ = true;
    readytoanalyze_ = true;
    Kmap_.clear();
    Kmap_["2"] = 0; // initially, zero shift
#ifdef IMATH_ROOT
    h_ = 0;
    h_nbins_ = 1024;
    h_precision_ = 0.05;
    if(h_file_ == 0){
      h_file_ = new TFile("imath.root","RECREATE");
      printf("recreating file imath.root\n");
    }
#endif
  }
  ~var_base(){
#ifdef IMATH_ROOT
    if(h_file_) {
      h_file_->ls();
      h_file_->Close();
      h_file_ = 0;
    }
#endif
  }
  
  std::string get_kstring();
  std::string get_name() {return name_;}
  var_base*   get_p1(){return p1_;}
  var_base*   get_p2(){return p2_;}
  double      get_fval(){return fval_;}
  long int    get_ival(){return ival_;}

  double get_minval(){return minval_;}
  double get_maxval(){return maxval_;}
  void   analyze();
#ifdef IMATH_ROOT
  TH2F*  get_h(){return h_;}
#endif
  void reset(){
    minval_ = std::numeric_limits<double>::max();
    maxval_ = -std::numeric_limits<double>::max();
#ifdef IMATH_ROOT
    h_->Clear();
#endif
  }    
  
  int get_nbits(){return nbits_;}
  std::map<std::string, int> get_Kmap(){return Kmap_;}
  double get_range(){return (1<<(nbits_-1))*K_;}// everything is signed
  double get_K()    {return K_;};
  int   get_shift(){return Kmap_["2"];}
  

  void makeready();
  int  get_step(){return step_;}
  int  get_latency(){return latency_;}
  bool calculate(int debug_level);
  bool calculate(){return calculate(1);}
  virtual void local_calculate(){}
  virtual void print(std::ofstream& fs, int l1=0, int l2=0){fs<<"// var_base here. Soemthing is wrong!! "<<l1<<", "<<l2<<"\n";}
  void print_step(int step, std::ofstream& fs);
  void print_all (std::ofstream& fs);
  
  void        dump_cout();
  std::string dump();
  std::string itos(int i);
  
 protected:
  std::string name_;
  var_base *p1_;
  var_base *p2_;
  std::string op_;    // operation
  int latency_;       // number of clock cycles for the operation (for HDL output)
  int step_;          // step number in the calculation (for HDL output)

  double fval_;
  long int ival_;

  int nbits_;
  double K_;
  std::map<std::string, int> Kmap_;  

  bool readytoanalyze_;
  bool readytoprint_;
  static const int MULT_LATENCY = 1;
  static const int LUT_LATENCY = 2;

  double minval_;
  double maxval_;
#ifdef IMATH_ROOT
  void set_hist_pars(int n = 256, double p=0.05){
    h_nbins_     = n;
    h_precision_ = p;
  }
  int    h_nbins_;
  double  h_precision_;
  static TFile* h_file_;
  TH2F *h_;
#endif
  
};

class var_adjustK : public var_base {

 public:

 var_adjustK(std::string name, var_base *p1, double Knew, double epsilon = 1e-5, bool do_assert = false, int nbits = -1):
  var_base(name,p1,0,0){
    op_ = "adjustK";
    K_     = p1->get_K();
    Kmap_  = p1->get_Kmap();
    
    double r = Knew / K_;

    lr_ = (r>1)? log2(r)+epsilon : log2(r);
    K_ = K_ * pow(2,lr_);
    if(do_assert) assert(fabs(Knew/K_ - 1)<epsilon);
    
    if(nbits>0)
      nbits_ = nbits;
    else
      nbits_ = p1->get_nbits()-lr_;

    Kmap_["2"] = Kmap_["2"] + lr_;
  }
  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);

 protected:
  int lr_;
};

class var_def : public var_base {

 public:

  //construct from scratch
 var_def(std::string name, std::string units, double fmax, double K):
  var_base(name,0,0,1){
    op_ = "def";
    K_    = K;
    nbits_ = log2(fmax / K) + 1.999999; //plus one to round up
    if(units!="")
      Kmap_[units] = 1;
    else{
      //defining a constant, K should be a power of two
      int l = log2(K);
      if(fabs(pow(2,l)/K-1)>1e-5){
	printf("defining unitless constant, yet K is not a power of 2! %g, %g\n",K, pow(2,l));
      }
      Kmap_["2"] = l;
    }
  }
  //construct from abother variable (all provenance info is lost!)
 var_def(std::string name, var_base *p):
  var_base(name,0,0,1){
    op_ = "def";
    K_     = p->get_K();
    nbits_ = p->get_nbits();
    Kmap_  = p->get_Kmap();
  }
  void    set_fval(double fval){
    fval_ = fval;
    ival_ = fval / K_;
  }
  void    set_ival(int ival){
    ival_ = ival;
    fval_ = ival * K_;
  }
  void print(std::ofstream& fs, int l1=0, int l2 = 0);
};


class var_add : public var_base {

 public:

 var_add(std::string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
  var_base(name,p1,p2,1){
    op_ = "add";

    std::map<std::string, int> map1 = p1->get_Kmap();
    std::map<std::string, int> map2 = p2->get_Kmap();
    int s1 = map1["2"];
    int s2 = map2["2"];
    
    //first check if the constants are all lined up
    //go over the two maps subtracting the units
    std::map<std::string,int>::iterator it;
    for(it=map2.begin(); it != map2.end(); ++it){
      if(map1.find(it->first) == map1.end())
	map1[it->first] = -it->second;
      else
	map1[it->first] = map1[it->first]-it->second;
    }
  
    //assert if different
    for(it=map1.begin(); it != map1.end(); ++it){
      if(it->second != 0){
	if(it->first != "2"){
	  printf("var_add: bad units! %s^%i for variable %s\n",
		 (it->first).c_str(),it->second,name_.c_str());
	  printf(" *********************************************************\n");
	  p1->dump_cout();
	  printf(" *********************************************************\n");
	  p2->dump_cout();
	  assert(0);
	}
      }
    }
    
    double ki1 = p1->get_K()/pow(2,s1);
    double ki2 = p2->get_K()/pow(2,s2);
    //those should be the same
    if(fabs(ki1/ki2-1.)>1e-6){
      printf("var_add: bad constants! %f %f for variable %s\n",
	     ki1,ki2,name_.c_str());
      printf(" *********************************************************\n");
      p1->dump_cout();
      printf(" *********************************************************\n");
      p2->dump_cout();
      assert(0);
    }
    //everything checks out!

    Kmap_ = p1->get_Kmap();
    
    int s0 = s1<s2? s1 : s2;
    shift1 = s1-s0;
    shift2 = s2-s0;

    int n1 = p1->get_nbits() + shift1;
    int n2 = p2->get_nbits() + shift2;
    int n0 = 1 + (n1>n2?n1:n2);

    //before shifting, check the range
    if(range > 0){
      n0 = log2(range/ki1/pow(2,s0))+1e-9;
      n0 = n0 + 2;
    }
    
    if(n0<=nmax){ //if it fits, we're done
      ps_ = 0;
      Kmap_["2"] = s0;
      nbits_ = n0;
    }
    else{
      ps_ = n0 - nmax;
      Kmap_["2"] = s0+ps_;
      nbits_ = nmax;
    }

    K_ = ki1 * pow(2,Kmap_["2"]);

  }  

  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);

 protected:
  int ps_;
  int shift1;
  int shift2;
  
};


class var_subtract : public var_base {

 public:

 var_subtract(std::string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
  var_base(name,p1,p2,1){
    op_ = "subtract";

    std::map<std::string, int> map1 = p1->get_Kmap();
    std::map<std::string, int> map2 = p2->get_Kmap();
    int s1 = map1["2"];
    int s2 = map2["2"];
    
    //first check if the constants are all lined up
    //go over the two maps subtracting the units
    std::map<std::string,int>::iterator it;
    for(it=map2.begin(); it != map2.end(); ++it){
      if(map1.find(it->first) == map1.end())
	map1[it->first] = -it->second;
      else
	map1[it->first] = map1[it->first]-it->second;
    }
  
    //assert if different
    for(it=map1.begin(); it != map1.end(); ++it){
      if(it->second != 0){
	if(it->first != "2"){
	  printf("var_add: bad units! %s^%i for variable %s\n",
		 (it->first).c_str(),it->second,name_.c_str());
	  printf(" *********************************************************\n");
	  p1->dump_cout();
	  printf(" *********************************************************\n");
	  p2->dump_cout();
	  assert(0);
	}
      }
    }
    
    double ki1 = p1->get_K()/pow(2,s1);
    double ki2 = p2->get_K()/pow(2,s2);
    //those should be the same
    if(fabs(ki1/ki2-1.)>1e-6){
      printf("var_add: bad constants! %f %f for variable %s\n",
	     ki1,ki2,name_.c_str());
      printf(" *********************************************************\n");
      p1->dump_cout();
      printf(" *********************************************************\n");
      p2->dump_cout();
      assert(0);
    }
    //everything checks out!

    Kmap_ = p1->get_Kmap();
    
    int s0 = s1<s2? s1 : s2;
    shift1 = s1-s0;
    shift2 = s2-s0;

    int n1 = p1->get_nbits() + shift1;
    int n2 = p2->get_nbits() + shift2;
    int n0 = 1 + (n1>n2?n1:n2);

    //before shifting, check the range
    if(range > 0){
      n0 = log2(range/ki1/pow(2,s0))+1e-9;
      n0 = n0 + 2;
    }
    
    if(n0<=nmax){ //if it fits, we're done
      ps_ = 0;
      Kmap_["2"] = s0;
      nbits_ = n0;
    }
    else{
      ps_ = n0 - nmax;
      Kmap_["2"] = s0+ps_;
      nbits_ = nmax;
    }

    K_ = ki1 * pow(2,Kmap_["2"]);
  }  

  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);

 protected:
  int ps_;
  int shift1;
  int shift2;
};

class var_nounits : public var_base {

 public:
 var_nounits(std::string name, var_base *p1, int ps = 17):
  var_base(name,p1,0,var_base::MULT_LATENCY){
    op_ = "nounits";
    ps_ = ps;
    nbits_ = p1->get_nbits();

    int s1 = p1->get_shift();
    double ki = p1->get_K()/pow(2,s1);
    int m = log2(ki);

    K_ = pow(2,s1+m);
    Kmap_["2"] = s1+m;
    double c = ki * pow(2,-m);
    cI_ = c * pow(2,ps_);
  }

  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);

 protected:
  int ps_;
  int cI_;
};

class var_shift : public var_base {
 public:
 var_shift(std::string name, var_base *p1, int shift):
  var_base(name,p1,0,0){
    op_    = "shift";
    shift_ = shift;

    nbits_ = p1->get_nbits()-shift;
    Kmap_  = p1->get_Kmap();
    K_     = p1->get_K();
  }
  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);
  
  protected:
    int shift_;
    
};

class var_neg : public var_base {
 public:
 var_neg(std::string name, var_base *p1):
  var_base(name,p1,0,1){
    op_    = "neg";
    nbits_ = p1->get_nbits();
    Kmap_  = p1->get_Kmap();
    K_     = p1->get_K();
  }
  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);
};


class var_timesC : public var_base {

 public:
 var_timesC(std::string name, var_base *p1, double cF, int ps = 17):
  var_base(name,p1,0,var_base::MULT_LATENCY){
    op_ = "timesC";
    cF_ = cF;
    ps_ = ps;

    nbits_ = p1->get_nbits();
    Kmap_  = p1->get_Kmap();
    K_     = p1->get_K();
    
    int   s1 = Kmap_["2"];
    double l  = log2(fabs(cF));
    if(l>0)
      l += 0.999999;
    int m    = l;

    cI_ = cF * pow(2,ps - m);
    K_  = K_ * pow(2, m);
    Kmap_["2"] = s1 + m;
  }

  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);

 protected:
  int ps_;
  int cI_;
  double cF_;
};

class var_mult : public var_base {
  
 public:
 var_mult(std::string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
  var_base(name,p1,p2,var_base::MULT_LATENCY){
    op_ = "mult";
    
    std::map<std::string,int> map1 = p1->get_Kmap();
    std::map<std::string,int> map2 = p2->get_Kmap();
    std::map<std::string,int>::iterator it;
    for(it=map1.begin(); it != map1.end(); ++it){
      if(Kmap_.find(it->first) == Kmap_.end())
	Kmap_[it->first] = it->second;
      else
	Kmap_[it->first] = Kmap_[it->first]+it->second;
    }
    for(it=map2.begin(); it != map2.end(); ++it){
      if(Kmap_.find(it->first) == Kmap_.end())
	Kmap_[it->first] = it->second;
      else
	Kmap_[it->first] = Kmap_[it->first]+it->second;
    }
    K_ = p1->get_K()*p2->get_K();
    int s0 = Kmap_["2"];

    int n0 = p1->get_nbits()+p2->get_nbits();
    if(range>0){
      n0 = log2(range/K_) + 1e-9;
      n0 = n0 + 2;
    }
    if(n0<nmax){
      ps_ = 0;
      nbits_ = n0;
    }
    else{
      ps_ = n0 - nmax;
      nbits_ = nmax;
      Kmap_["2"] = s0 + ps_;
      K_ = K_ * pow(2,ps_);
    }
  }
  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);

 protected:
  int ps_;
};
    
class var_inv : public var_base {

 public:
  enum mode {pos, neg, both};
  
 var_inv(std::string name, var_base *p1, double offset, int nbits, int n, unsigned int shift, mode m, int nbaddr = -1):
  var_base(name,p1,0,var_base::LUT_LATENCY){
    op_ = "inv";
    offset_ = offset;
    nbits_  = nbits;
    n_      = n;
    shift_  = shift;
    m_      = m;
    if(nbaddr <0) nbaddr = p1->get_nbits();
    nbaddr_ = nbaddr - shift;
    if(m_!=mode::both) nbaddr_--;
    Nelements_ = 1<<nbaddr_;
    mask_      = Nelements_-1;
    ashift_    = sizeof(int)*8-nbaddr_;

    std::map<std::string, int> map1 = p1->get_Kmap();
    std::map<std::string,int>::iterator it;
    for(it=map1.begin(); it != map1.end(); ++it)
      Kmap_[it->first] = -it->second;
    Kmap_["2"] = Kmap_["2"] - n;
    K_ = pow(2,-n)/p1->get_K();

    LUT = new int[Nelements_];
    unsigned int ms = sizeof(int)*8-nbits_;
    int offsetI = offset_ / p1_->get_K();
    for(int i=0; i<Nelements_; ++i){
      int i1 = addr_to_ival(i);
      if(offsetI+i1)
	LUT[i] = (round_int((1<<n_)/(offsetI+i1))<<ms)>>ms;
      else
	LUT[i] = 0;
    }
  }
  ~var_inv(){
    if(LUT) delete LUT;
  }

  void initLUT(double offset);  
  
  void local_calculate();
  void print(std::ofstream& fs, int l1=0, int l2 = 0);
  void writeLUT(std::ofstream& fs);

  int ival_to_addr(int ival){ return ((ival>>shift_)&mask_);}
  int addr_to_ival(int addr){
    switch(m_){
    case mode::pos  :  return addr<<shift_;
    case mode::neg  :  return (addr-Nelements_)<<shift_;
    case mode::both :  return (addr<<ashift_)>>(ashift_-shift_);
    }
  }

  int round_int( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
  }

 protected:
  double offset_;
  int n_;
  mode m_;
  unsigned int shift_;
  unsigned int mask_;
  unsigned int ashift_;
  int Nelements_;
  int nbaddr_;
  
  int *LUT;
};




#endif
