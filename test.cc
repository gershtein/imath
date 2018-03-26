// ROOT includes
#include "TMath.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TRandom3.h>

#include <iostream>

#include "imath.h"

// a hack to permit the static file handle
#ifdef IMATH_ROOT
TFile* var_base::h_file_=0;
#endif

using namespace std;

int main()
{

  TRandom3 rr;
  
  static var_def  v1("v1","m",40, 0.02);
  static var_def  v2("v2","m",40, 0.02);
  static var_def  v3("v3","",10,  1./(1<<15));
  static var_inv  v4("v4",&v1, 0, 18, 25, 0, var_inv::mode::both);
  static var_mult v5("v5",&v2,&v4);
  // the analysis of v5 says the range by default is too large
  // so one can put it in explicitely like this:
  // static var_mult v5("v5",&v2,&v4, 3.);
  static var_add  v6("v6",&v3,&v5);
  static var_adjustK v7("v7",&v6,v6.get_K()/4.);
  
  
  for(int i=0; i<100000; ++i){
    float fv1 = 10+rr.Rndm()*27.;
    float fv2 = 20+10*rr.Rndm();
    if(rr.Rndm()>0.5) fv2 = -fv2;
    float fv3 = 6+rr.Rndm()*3.5;
    v1.set_fval(fv1);
    v2.set_fval(fv2);
    v3.set_fval(fv3);
    v7.calculate(1);
  }

  ofstream fs("test.out");
  
  v7.print_all(fs);
  v7.analyze();

  
};
