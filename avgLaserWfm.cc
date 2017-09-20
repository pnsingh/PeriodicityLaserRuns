#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"

#include "TNtuple.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TChain.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TError.h"
#include "TPDF.h"
#include "TPad.h"
#include "TProfile2D.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "TLine.h"
#include "TLatex.h"

#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

#include "fftw3.h"
#include <numeric>      // std::iota
using namespace std;

size_t next_fftw_size (int dw) {
    int bases[4] = { 2, 3, 5, 7 };
    int min_d = 0, min_b = 0;
    for (int b=0;b<4;b++) {
      int k = pow (b, ceil (log (dw) /log (b))) - dw;
      if (b == 2 || k < min_d) {
        min_d = k;
        min_b = b;
      }
    }
    return pow (min_b, ceil (log (dw) /log (min_b)));
  }


int main(int argc,char *argv[]){
  
  TApplication myapp("myapp", 0, 0);
  
  gStyle->SetOptStat("neMRou");
  
  int oneChannel=1;

  
  TFile inF(argv[1]);
//   TFile inF("/home/yguardin/Work/sensei/theCode/procMacros/out/catalog.root");
//   TTree *hitSumm = (TTree*)inF.Get("hitSumm");
  TTree *Events = (TTree*)inF.Get("treeBuilderAvgWfmLaser/Events");
   Events->Print(); //Prints whats in the tree
  TH2D *h2_wfm_6 = (TH2D*)inF.Get("treeBuilderAvgWfmLaser/h2_wfm_6");
  TCanvas *c0;
//   h2_wfm->Print();
  h2_wfm_6->Draw("colz");
  h2_wfm_6->GetYaxis()->SetTitle("Signal [ADC]");
  h2_wfm_6->GetXaxis()->SetTitle("Time [4ns]");
//   Events->ls();
  Events->Print();
//   cout <<hitSumm->GetEntries()<<endl;
  
  //Create a cannonical wfm with period w
  double period=0.501;
  vector<double> sinDiscrete;
//   cout <<TMath::ASin(1)<<endl;
//   cout<<sin(TMath::ASin(1)/2)<<endl;
//   cout<<sin((5.0/5)*TMath::ASin(-1)/2)<<endl;
  for(int i=0;i<1500;i++){
//     cout<<sin(TMath::ASin(1))<<endl;
//     cout <<sin(double(i)*TMath::ASin(1)/period)<<endl;
    sinDiscrete.push_back(1*sin(double(i)*TMath::ASin(1)/period));
  }
  
  
  //vector< vector<double> > *wfm1_mSamp1 =new std::vector< vector<double> >;
  vector< vector<double> > *wfm_mSamp1= new vector< vector<double> >;
  Events->SetBranchAddress("wfm_mSamp1", &wfm_mSamp1);
  vector< vector<double> > *wfm2_mSamp1= new vector< vector<double> >;
  Events->SetBranchAddress("wfm2_mSamp1", &wfm2_mSamp1);
  Events->GetEntry(0);
  cout <<wfm_mSamp1->size()<<" "<<wfm_mSamp1->at(0).size()<<endl;
  TGraphErrors *gch[38];
  TGraphErrors *gch_rms[38];
  for(int i=0; i<wfm_mSamp1->size();i++){
    gch[i] = new TGraphErrors();
    gch_rms[i] = new TGraphErrors();
  }
  double nWfms=100000;
  
  for(int j=0;j<wfm_mSamp1->size();j++){
    for(int i =0;i<wfm_mSamp1->at(j).size();i++){
      double mean=wfm_mSamp1->at(j).at(i)/nWfms;
      double e_mean=sqrt((wfm2_mSamp1->at(j).at(i)/nWfms)-pow(mean,2))/sqrt(nWfms);
      gch[j]->SetPoint(gch[j]->GetN(),i,mean);
      gch[j]->SetPointError(gch[j]->GetN()-1,0,e_mean);
      gch_rms[j]->SetPoint(gch_rms[j]->GetN(),i,sqrt((wfm2_mSamp1->at(j).at(i)/nWfms)-pow(mean,2)));
       gch[j]->SetPoint(gch[j]->GetN(),i,sinDiscrete[i]);
       gch[j]->SetPointError(gch[j]->GetN()-1,0,e_mean);

    }
  }
  
  TCanvas *c1[38];
  
  for(int i=0;i<wfm_mSamp1->size();i++){
    int ch=i;
    if(i==oneChannel){
      c1[i] = new TCanvas(Form("ch%d",ch),"c1",0,8,1300,800);
      gch[i]->Draw("APL");
      gch[i]->GetXaxis()->SetTitle("Time [4ns]");
      gch[i]->GetXaxis()->SetTitle("Amplitude [ADC]");
      gch[i]->SetTitle("Waveform");
      gch[i]->SetMarkerStyle(20);
    }
  }
  
  TCanvas *c2[38];
  
  for(int i=0;i<wfm_mSamp1->size();i++){
    int ch=i;
    if(i==oneChannel){
      c2[i] = new TCanvas(Form("ch_rms%d",ch),"c2",0,8,1300,800);
      gch_rms[i]->Draw("APL");
      gch_rms[i]->GetXaxis()->SetTitle("Time [4ns]");
      gch_rms[i]->GetXaxis()->SetTitle("Amplitude [ADC]");
      gch_rms[i]->SetTitle("RMS");
      gch[i]->SetMarkerStyle(20);
    }
  }
  
  
  ////////////
  int X=5;
  vector<double> signalEveryX(1500,0);
  double x, y;
  for(int i=0;i<wfm_mSamp1->at(oneChannel).size();i++){

    gch[oneChannel]->GetPoint(i,x,y);
//     cout <<i<<" "<<i/30<<" "<<floor(i/30)<<endl;
    //signalEveryX[floor(i/X)]+=y;
    for(int j=0;j<5;j++){
      gch[oneChannel]->GetPoint(i+j,x,y);
      signalEveryX[i]+=y;
    }
  }
  
  TGraphErrors *gAvg= new TGraphErrors();
  for(int i =0;i<signalEveryX.size();i++){
//      gAvg->SetPoint(gAvg->GetN(),i*X+(X-1)*1.0/2,signalEveryX[i]);
    gAvg->SetPoint(gAvg->GetN(),i,signalEveryX[i]);
  }
  TCanvas *c3=new TCanvas("c3","c3",0,8,1300,800);;
  gAvg->Draw("APL");
  gAvg->SetMarkerStyle(20);
  
  size_t fft_n_;
  double fft_fs_, *fft_in_, *fft_out_;
  fftw_plan fft_plan_;
  std::vector<double> frequencies_, power_spectrum_;
  size_t power_spectrum_entries_;
  int bits_;
  
  if (frequencies_.size ()) std::fill (fft_in_, fft_in_ + fft_n_, 0);
  bits_ = 100000;
  fft_fs_= 250000000;//Hz
//   fft_fs_= 1;//Hz
  
  fftw_cleanup();
  
  int dw=1000;
  
  fft_n_ = next_fftw_size (dw);//Me fijo en los primeros 2 us
  
  
  
  fft_in_ = (double *)fftw_malloc(sizeof(double) * fft_n_);
  std::fill (fft_in_, fft_in_ + fft_n_, 0);
  fft_out_ = (double *)fftw_malloc(sizeof(double) * fft_n_);
  fft_plan_ = fftw_plan_r2r_1d (fft_n_, fft_in_, fft_out_, FFTW_R2HC, FFTW_ESTIMATE);
  std::cerr << "DW = " << dw << " N = " << fft_n_ << " FS = " << fft_fs_ << " BITS = " << bits_ << std::endl;
  
  power_spectrum_.resize (fft_n_ / 2 + 1);
  frequencies_.resize (fft_n_ / 2 + 1);
  //std::iota (frequencies_.begin (), frequencies_.end (), 0);
  for(int i=0;i<frequencies_.size();i++)frequencies_[i]=i;
  double fs_step = fft_fs_ / 2. / 1e6 / frequencies_.back ();
  
  cout <<fs_step<<endl;
  //   std::for_each (frequencies_.begin(), frequencies_.end (), [&fs_step](Double_t & d) { d *= fs_step; });
  for(int i=0;i<frequencies_.size();i++)frequencies_[i]=frequencies_[i]*fs_step;
   
  power_spectrum_entries_ = 0; 
  TGraph *graph_;
  graph_ = new TGraph(frequencies_.size ());
  std::copy (frequencies_.begin (), frequencies_.end (), graph_->GetX ());
  graph_->GetXaxis()->SetTitle("Frequency (MHz)");
  graph_->GetYaxis()->SetTitle("Amplitude [ADC]");
  graph_->SetTitle("Spectrum");
  //double mwHz_scale = 2. / (1 << bits_) * std::sqrt (1 / (dw * fft_fs_ / 2 * 50e-3));
  for (size_t i = 0; i < dw; i++) {
//     fft_in_[i] = sinDiscrete[i];
    fft_in_[i] = wfm_mSamp1->at(oneChannel).at(i)/nWfms;
  }
  
   fftw_execute (fft_plan_);
   
  double sum = 0;
  for (size_t i = 1; i < 20; i++) sum += fft_out_[i] * fft_out_[i] + fft_out_[fft_n_ - i] * fft_out_[fft_n_ - i];
  cout <<sum<<endl;
  //if (sum > 1e-6) return 0;

  if (power_spectrum_entries_ > power_spectrum_.size () / 10 && power_spectrum_entries_ > 1000) power_spectrum_entries_ = 0;
  power_spectrum_[0] = (power_spectrum_[0] * power_spectrum_entries_ + fft_out_[0] * fft_out_[0]) / (power_spectrum_entries_ + 1);  // DC component

  for (size_t i = 1; i < (fft_n_ + 1) / 2; i++) power_spectrum_[i] = (power_spectrum_[i] * power_spectrum_entries_ + (fft_out_[i] * fft_out_[i] + fft_out_[fft_n_ - i] * fft_out_[fft_n_ - i])) / (power_spectrum_entries_ + 1);
  cout<<"---------------------";
  cout<<"\n";
  cout<<power_spectrum_entries_<<"\n";
  if (fft_n_ % 2 == 0) power_spectrum_[fft_n_ / 2] = (power_spectrum_[fft_n_ / 2] * power_spectrum_entries_ + fft_out_[fft_n_ / 2] * fft_out_[fft_n_ / 2]) / (power_spectrum_entries_ + 1);  // Nyquist freq
  
//   for (size_t i = 0; i < power_spectrum_.size (); i++) graph_->GetY ()[i] = 10 * log10 (power_spectrum_[i]);
  double max=-1;
  double min=1000;
  for (size_t i = 0; i < power_spectrum_.size (); i++){
    graph_->GetY ()[i] = sqrt(power_spectrum_[i])/(dw/2);
    if(sqrt(power_spectrum_[i])/(dw/2)>max)max=sqrt(power_spectrum_[i])/(dw/2);
    if(sqrt(power_spectrum_[i])/(dw/2)<min)min=sqrt(power_spectrum_[i])/(dw/2);
  }
  
  TCanvas *cFFT=new TCanvas("fft","fft",0,8,1300,800);
//   cFFT->SetLogy();
   graph_->GetYaxis()->SetRangeUser (min-(max-min)*0.1, max+(max-min)*0.1);
  graph_->Draw ("APL");
  
  
  double t_per_sample=0.004;//ns
  int lala=201;
  TLine *lines[lala];
  TLatex *texts[lala];
  for(int i=2;i<lala;i++){
    if(i<7 || i==10  || i==60 || i==20 || i==15){
      lines[i] =new TLine(1.0/(i*t_per_sample),min-(max-min)*0.1,1.0/(i*t_per_sample),max+(max-min)*0.1);
      lines[i]->SetLineColor(kRed);
      lines[i]->Draw("same");

      texts[i] = new TLatex(1.0/(i*t_per_sample),min+(max-min)*0.85,Form("%d samples",i));
      texts[i]->SetTextSize(0.03);
      texts[i]->SetTextAngle(90);
      texts[i]->SetTextColor(kRed);
      texts[i]->Draw("same");
    }

  }
  

  
  graph_->SetLineWidth(3);
  graph_->Draw ("PL");
  
  
//   TH2F *ejes=new TH2F("ejes","",1000,0,10,1000,-10,10);
//   ejes->SetXTitle("");
//   ejes->SetYTitle("");
//   ejes->SetTitleOffset(1.2,"Y");
//   ejes->SetStats(0);
//   ejes->Draw();


  
 
    myapp.Run();
//   outF.Close();  

  return 0;
}