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
#include <numeric>

using namespace std;

double nWfms=10000;

int main(int argc,char *argv[])

{
    ofstream f,f1,f2;
    f.open("log.txt");
    f1.open("log1.txt");
    //f2.open("log2.txt");
    double amp[750];   
    TApplication myapp("myapp", 0, 0);
    gStyle->SetOptStat("neMRou");
    TFile inF(argv[1]);
    TTree *Events = (TTree*)inF.Get("treeBuilderAvgWfmLaser/Events");
    vector< vector<double> > *wfm_mSamp1= new vector< vector<double> >;
    Events->SetBranchAddress("wfm_mSamp1", &wfm_mSamp1);
    vector< vector<double> > *wfm2_mSamp1= new vector< vector<double> >;
    Events->SetBranchAddress("wfm2_mSamp1", &wfm2_mSamp1);
    Events->GetEntry(0);

    double fft_fs_, *fft_in_, *fft_out_;
    int size=1500;
    //fft_n_ = next_fftw_size (dw);
    fft_in_ = (double *)fftw_malloc(sizeof(double) * size);
    //std::fill (fft_in_, fft_in_ + fft_n_, 0);
    fft_out_ = (double *)fftw_malloc(sizeof(double) * size);


    vector<double> sinDiscrete;
    for(int i=0;i<1500;i++)
    {
        sinDiscrete.push_back(sin(double(i)/10.0));
        //f1<<double(i)<<"\t"<<sinDiscrete[i]<<"\n";
    }    
    for (int i = 0; i < size; i++) 
    {
     fft_in_[i] = wfm_mSamp1->at(20).at(i)/nWfms;
     f1<<i<<"\t"<<fft_in_[i]<<"\n";
     //fft_in_[i] = sinDiscrete[i];
    }
    // Variable Declaration
    //double array[] = {0.1, 0.6, 0.1, 0.4, 0.5, 0, 0.8, 0.7, 0.8, 0.6, 0.1,0};
    double *out,*mag,*phase;
    double real,imag;
    int i,j;
    fftw_complex *out_cpx, *mod_cpx;
    fftw_plan fft; 
    fftw_plan ifft; 
    
    //Allocate Memory
    //out_cpx = (fftw_complex*) fftw_malloc(sizeof(fft_n_)*(size/2+1));
    out_cpx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (size));
    mod_cpx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* (size));
    out = (double *) malloc(size*sizeof(double));
    mag = (double *) malloc(size*sizeof(double));
    phase = (double *) malloc(size*sizeof(double));
   
    //fft = fftw_plan_dft_r2c_1d(size, array, out_cpx, FFTW_ESTIMATE);  //Setup fftw plan for fft 1 dimensional, real signal
    //ifft = fftw_plan_dft_c2r_1d(size, mod_cpx, out, FFTW_ESTIMATE);   //Setup fftw plan for ifft 1 dimensional, complex signal

    fft = fftw_plan_dft_r2c_1d(size, fft_in_, out_cpx, FFTW_ESTIMATE);  //Setup fftw plan for fft 1 dimensional, real signal
    ifft = fftw_plan_dft_c2r_1d(size, mod_cpx, out, FFTW_ESTIMATE);   //Setup fftw plan for ifft 1 dimensional, complex signal

    fftw_execute(fft);	//perform fft
	for(j=0;j<size/2+1;j++)
	{
        real = out_cpx[j][0];	//Extract real component
		    imag = out_cpx[j][1];   //Extract imaginary component
		    mag[j] = sqrt((real*real)+(imag*imag));  // Calculate the Magnitude
		    phase[j] = atan2(imag,real); // Calculate the Phase
        amp[j]=(mag[j]/(double(size)/2.0));          
        f<<(double(j)/6.0)<<"\t"<<phase[j]<<"\t"<<amp[j]<<"\n";
        //f<<(double(j)/6.0)<<"\t"<<phase[j]<<"\t"<<sqrt(mag[j])*(double(j)/(double(size)))<<"\n";	
    }
    f.close();
    f1.close();
	cout<<"Done\n";
	/***********MODIFICATION****************************/
	//You can perform frequency domain modification here  	
	/***************************************************/	

	for(j=0;j<size/2+1;j++)
	{
		mod_cpx[j][0] = (mag[j]*cos(phase[j]));  //Construct new real component
		mod_cpx[j][1] = (mag[j]*sin(phase[j]));  //Construct new imaginary  component
        //cout<<j<<"\t"<<mod_cpx[j][0]<<"\t"<<mod_cpx[j][1]<<"\t";
        //cout<<"\n";
	}

    fftw_execute(ifft); //perform ifft
    // Print input and output
    cout<<("Input:    Output:");
    cout<<"\n";
    for(i=0;i<size;i++)
    {
	out[i] = out[i]/size;
    //cout<<i<<"\t"<<out[i]<<"\n";	
	//cout<<i<<"\t"<<out[i]<<"\n";
    //f2<<i<<"\t"<<out[i]<<"\n";
    //printf("%ft%fn",(array[i]),out[i]);
    }

   // f2.close();

   TCanvas *c = new TCanvas("c","Complete FT",0,0,600,400);
   TGraph2D *dt = new TGraph2D("log.txt");
   dt->SetTitle("Fourier Transform");
   dt->GetXaxis()->SetTitle("Frequency (MHz)");
   dt->GetYaxis()->SetTitle("Phase (Radians)");
   dt->GetZaxis()->SetTitle("Amplitude (ADC)");
   dt->Draw("colz");

   TCanvas *c1 = new TCanvas("c1","Waveform",0,0,600,400);
   TGraph *g1 = new TGraph("log1.txt");
   g1->SetTitle("Waveform");
   g1->GetXaxis()->SetTitle("Samples (4ns)");
   g1->GetYaxis()->SetTitle("Amplitude (ADC)");
   g1->Draw();


   TCanvas *c2 = new TCanvas("c2","FFT",0,0,600,400);
   TGraph *g2 = new TGraph("log.txt","%lg %*lg %lg");
   g2->SetTitle("FFT");
   g2->GetXaxis()->SetTitle("Frequency (MHz)");
   g2->GetYaxis()->SetTitle("Amplitude (ADC)");
   g2->Draw();
/*
   TCanvas *c3 = new TCanvas("c3","TestEx",0,0,600,400);
   TGraph *g3 = new TGraph("log2.txt");
   g3->SetTitle("Comparison");
   g3->GetXaxis()->SetTitle("Input");
   g3->GetYaxis()->SetTitle("Output");
   g3->Draw();
*/   



    // Free all memory
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    fftw_free(out_cpx);
    fftw_free(mod_cpx);
    free(out);
    free(mag);
    free(phase);

    myapp.Run();
    return 0;
}