//
// Integer representation of floating point arithmetic suitable for FPGA designs
// 
// Author: Yuri Gershtein 
// Date:   March 2018
//

#include "imath.h"

std::string var_base::itos(int i) 
{
    std::ostringstream os;
    os << i;
    return os.str();
}

std::string var_base::get_kstring()
{

  char s[1024];
  std::string t="";
  std::map<std::string,int>::iterator it;
  for(it = Kmap_.begin(); it != Kmap_.end(); ++it){
    sprintf(s,"^(%i)",it->second);
    std::string t0(s);
    t = t + it->first + t0;
  }

  return t;
}

void var_base::analyze()
{
  if(!readytoanalyze_) return;
  
  double u = maxval_;
  if(u < -minval_)
      u = -minval_;

  int iu = log2(get_range()/u);
  if(iu>1){
    printf("analyzing %s: range %g is much larger then %g. suggest cutting by a factor of 2^%i\n",name_.c_str(),get_range(),u,iu);
  }
#ifdef IMATH_ROOT
  if(h_){
    double eff = h_->Integral()/h_->GetEntries();
    if(eff<0.99) {
      printf("analyzing %s: range is too small, contains %f\n",name_.c_str(),eff);
      h_->Print();
    }
    h_file_->cd();
    TCanvas *c = new TCanvas();
    c->cd();
    h_->Draw("colz");
    h_->Write();
  }
  else{
    printf("analyzing %s: no histogram!\n",name_.c_str());
  }
#endif

  if(p1_) p1_->analyze();
  if(p2_) p2_->analyze();

  readytoanalyze_ = false;
}

std::string var_base::dump()
{
  char s[1024];
  std::string u = get_kstring();
  sprintf(s,"Name = %s \t Op = %s \t nbits = %i \n       ival = %li \t fval = %g \t K = %g Range = %f\n       units = %s\n",  
	  name_.c_str(), op_.c_str(), nbits_, ival_, fval_, K_, get_range(), u.c_str());
  std::string t(s);
  return t;
}

void var_base::dump_cout()
{
  char s[2048];
  std::string u = get_kstring();
  sprintf(s,"Name = %s \t Op = %s \t nbits = %i \n       ival = %li \t fval = %g \t K = %g Range = %f\n       units = %s\n       step = %i, latency = %i\n", name_.c_str(), op_.c_str(), nbits_, ival_, fval_, K_, get_range(), u.c_str(), step_, latency_);
  std::string t(s);
  std::cout<<t;
  if(p1_) p1_->dump_cout();
  if(p2_) p2_->dump_cout();
}

bool var_base::calculate(int debug_level)
{
  bool ok1 = true;
  bool ok2 = true;
  
  if(p1_) ok1 = p1_->calculate(debug_level);
  if(p2_) ok2 = p2_->calculate(debug_level);
  
  local_calculate();

  bool all_ok = ok1 && ok2;

  if(debug_level == 3 || all_ok){
    if(fval_ > maxval_) maxval_ = fval_;
    if(fval_ < minval_) minval_ = fval_;
#ifdef IMATH_ROOT
    if(h_==0){
      h_file_->cd();
      std::string hname = "h_"+name_;
      h_ = (TH2F*) h_file_->Get(hname.c_str());
      if(h_ == 0){
	std::string st = name_+";fval;fval-ival*K";
	h_ = new TH2F(hname.c_str(),name_.c_str(),
		      h_nbins_,-get_range(), get_range(),
		      h_nbins_,-get_range()*h_precision_, get_range()*h_precision_);
	if(debug_level==3) std::cout<<" booking histogram "<<hname<<"\n";
      }
    }
    h_->Fill(fval_, K_*ival_-fval_);
#endif
  }
  
  bool todump = false;
  int nmax = sizeof(long int)*8;
  int ns = nmax - nbits_;
  long int itest = ival_;
  itest = itest<<ns;
  itest = itest>>ns;
  if(itest!=ival_){
    if(debug_level == 3 || all_ok){
      std::cout<<"imath: truncated value mismatch!! "<<ival_<<" != "<<itest<<"\n";
      todump = true;
    }
    all_ok = false;
  }
  
  float ftest = ival_ * K_;
  float tolerance = 0.1 * fabs(fval_);
  if(tolerance < 2 * K_) tolerance = 2 * K_;
  if(fabs(ftest-fval_)> tolerance){
    if(debug_level == 3 || (all_ok && (op_!="inv" ||debug_level>=2 ))){
      std::cout<<"imath: **GROSS** value mismatch!! "<<fval_<<" != "<<ftest<<"\n";
      todump = true;
    }
    all_ok = false;
  }
  
  if(todump)
    std::cout<<dump();

  return all_ok;
}

void var_inv::writeLUT(std:: ofstream& fs)
{
  for(int i=0; i<Nelements_; ++i){
    fs<<LUT[i]<<"\n";
  }
}

void var_adjustK::adjust(double Knew, double epsilon, bool do_assert, int nbits)
{
  //WARNING!!!
  //THIS METHID CAN BE USED ONLY FOR THE FINAL ANSWER
  //THE CHANGE IN CONSTANT CAN NOT BE PROPAGATED UP THE CALCULATION TREE
  
    K_     = p1_->get_K();
    Kmap_  = p1_->get_Kmap();
    double r = Knew / K_;

    lr_ = (r>1)? log2(r)+epsilon : log2(r);
    K_ = K_ * pow(2,lr_);
    if(do_assert) assert(fabs(Knew/K_ - 1)<epsilon);
    
    if(nbits>0)
      nbits_ = nbits;
    else
      nbits_ = p1_->get_nbits()-lr_;

    Kmap_["2"] = Kmap_["2"] + lr_;
    
}

//
//  local calculations
//

void var_adjustK::local_calculate()
{
  fval_ = p1_->get_fval();
  ival_ = p1_->get_ival();
  if(lr_>0)
    ival_ = ival_ >> lr_;
  else if(lr_<0)
    ival_ = ival_ <<(-lr_);
}
void var_add::local_calculate()
{
  fval_ = p1_->get_fval() + p2_->get_fval();
  long int i1 = p1_->get_ival();
  long int i2 = p2_->get_ival();
  if(shift1>0) i1 = i1 << shift1;
  if(shift2>0) i2 = i2 << shift2;
  ival_ = i1 + i2;
  if(ps_>0) ival_ = ival_ >> ps_;
}
void var_subtract::local_calculate()
{
  fval_ = p1_->get_fval() - p2_->get_fval();
  long int i1 = p1_->get_ival();
  long int i2 = p2_->get_ival();
  if(shift1>0) i1 = i1 << shift1;
  if(shift2>0) i2 = i2 << shift2;
  ival_ = i1 - i2;
  if(ps_>0) ival_ = ival_ >> ps_;
}
void var_nounits::local_calculate()
{
  fval_ = p1_->get_fval();
  ival_ = (p1_->get_ival() * cI_)>>ps_;
}
void var_timesC::local_calculate()
{
  fval_ = p1_->get_fval() * cF_;
  ival_ = (p1_->get_ival() * cI_)>>ps_;
}
void var_neg::local_calculate()
{
  fval_ = -p1_->get_fval();
  ival_ = -p1_->get_ival();
}
void var_shift::local_calculate()
{
  fval_ = p1_->get_fval() * pow(2,-shift_);
  ival_ = p1_->get_ival();
  if(shift_>0) ival_ = ival_>>shift_;
  if(shift_<0) ival_ = ival_<<(-shift_);
}
void var_mult::local_calculate()
{
  fval_ = p1_->get_fval() * p2_->get_fval();
  ival_ = (p1_->get_ival() * p2_->get_ival())>>ps_;
}
void var_inv::initLUT(double offset)
{
  offset_ = offset;
  unsigned int ms = sizeof(int)*8-nbits_;
  int offsetI = round_int(offset_ / p1_->get_K());
  for(int i=0; i<Nelements_; ++i){
    int i1 = addr_to_ival(i);
    if(offsetI+i1)
      LUT[i] = (round_int((1<<n_)/(offsetI+i1))<<ms)>>ms;
    else
      LUT[i] = 0;
  }
}
void var_inv::local_calculate()
{
  fval_ = 1./(offset_ + p1_->get_fval());
  ival_ = LUT[ival_to_addr(p1_->get_ival())];
}
//
// print functions
//

void var_base::makeready()
{
  readytoprint_ = true;
  if(p1_) p1_->makeready();
  if(p2_) p2_->makeready();
}

void var_base::print_step(int step, std::ofstream& fs){
  if(!readytoprint_) return;
  if(step > step_) return;
  int l1 = 0;
  int l2 = 0;
  if(p1_) {
    p1_->print_step(step, fs);
    l1 = step_ - p1_->get_latency() - p1_->get_step();
  }
  if(p2_) {
    p2_->print_step(step, fs);
    l2 = step_ - p2_->get_latency() - p2_->get_step();
  }
  if(step==step_){
    if(l1<0 || l2<0 || (l1>0&&l2>0) ){
      printf("%s::print_step(%i): something wrong with latencies! %i %i\n",name_.c_str(),step,l1, l2);
      dump_cout();
      assert(0);      
    }
    if(l1>0) fs<<"pipe_delay("+p1_->get_name()+", "+p1_->get_name()+"_delay"+itos(l1)+", "+itos(l1)+");\n";
    if(l2>0) fs<<"pipe_delay("+p2_->get_name()+", "+p2_->get_name()+"_delay"+itos(l2)+", "+itos(l2)+");\n";
  
    print(fs, l1, l2);
    readytoprint_ = false;
  }
}

void var_base::print_all(std::ofstream& fs)
{
  for(int i=0; i<=step_; ++i){
    fs<<"//\n// STEP "<<i<<"\n\n";
    print_step(i,fs);
  }
}

void var_adjustK::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  assert(l2==0);

  std::string shift = "";
  if(lr_>0)
    shift = " >> " + itos(lr_);
  else if(lr_<0)
    shift = " << " + itos(-lr_);

  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string t = name_ + " = " + n1 + shift;

  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n";
  
}

void var_def::print(std::ofstream& fs, int l1, int l2)
{
  assert(l1==0);
  assert(l2==0);
  std::string t = name_ + " = " + itos(ival_);
  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}

void var_add::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  assert(p2_);
  std::string t = p1_->get_name();
  if(l1>0) t = t+"_delay"+itos(l1);
  if(shift1>0) t += "<<"+itos(shift1);
  t += " + " + p2_->get_name();
  if(l2>0) t = t+"_delay"+itos(l2);
  if(shift2>0) t += "<<"+itos(shift2);
  if(ps_>0) t = "("+t+")>>"+itos(ps_);
  t = name_ + " = " + t;
  
  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}

void var_subtract::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  assert(p2_);
  std::string t = p1_->get_name();
  if(l1>0) t = t+"_delay"+itos(l1);
  if(shift1>0) t += "<<"+itos(shift1);
  t += " - " + p2_->get_name();
  if(l2>0) t = t+"_delay"+itos(l2);
  if(shift2>0) t += "<<"+itos(shift2);
  if(ps_>0) t = "("+t+")>>"+itos(ps_);
  t = name_ + " = " + t;

  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}

void var_nounits::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  assert(l2==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string t = name_ + " = (" + n1 + " * " + itos(cI_) + ")";
  if(ps_>0) t = t + ">>" + itos(ps_);
  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}

void var_timesC::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  assert(l2==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string t = name_ + " = (" + n1 + " * " + itos(cI_) + ")";
  if(ps_>0) t = t + ">>" + itos(ps_);
  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}
void var_neg::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  assert(l2==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string t = name_ + " = - " + n1;
  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}
void var_shift::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  assert(l2==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string t = name_ + " = " + n1;
  if(shift_>0) t = t + ">>" + itos(shift_);
  if(shift_<0) t = t + "<<" + itos(-shift_);
  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}
void var_mult::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  assert(p2_);
  std::string n2 = p2_->get_name();
  if(l2>0) n2 = n2 + "_delay"+itos(l2);
  std::string t = name_ + " = (" + n1 + " * " + n2 + ")";
  if(ps_>0) t = t + ">>" + itos(ps_);
  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}
void var_inv::print(std::ofstream& fs, int l1, int l2)
{
  assert(p1_);
  assert(l2==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  //first calculate address
  std::string t1 = "addr_" + name_;
  std::string t = t1 + " = ";
  if(shift_>0)
    t = t + "(" + n1 + ">>"+itos(shift_)+") & "+itos(mask_);
  else
    t = t + n1 + " & "+itos(mask_);
  fs<<t<<"; // address for the LUT\n"; 
  
  std::string t2 = "LUT_" + name_;
  t = name_ + " = " + t2 + "[" + t1 + "]";
  fs<<t<<"; "<<std::string(40-t.size(),' ')<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
}


