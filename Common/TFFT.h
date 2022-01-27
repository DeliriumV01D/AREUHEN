#pragma once

#include "TVirtualFFT.h"

// - "C2CFORWARD" - a complex input/output discrete Fourier transform (DFT) 
//                  in one or more dimensions, -1 in the exponent
// - "C2CBACKWARD"- a complex input/output discrete Fourier transform (DFT) 
//                  in one or more dimensions, +1 in the exponent
// - "R2C"        - a real-input/complex-output discrete Fourier transform (DFT)
//                  in one or more dimensions,
// - "C2R"        - inverse transforms to "R2C", taking complex input 
//                  (storing the non-redundant half of a logically Hermitian array) 
//                  to real output
// - "R2HC"       - a real-input DFT with output in ?Ehalfcomplex?E format, 
//                  i.e. real and imaginary parts for a transform of size n stored as
//                  r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
// - "HC2R"       - computes the reverse of FFTW_R2HC, above
// - "DHT"        - computes a discrete Hartley transform
//Returns a pointer to the FFT of requested size and type.
//Parameters:
// -ndim : number of transform dimensions
// -n    : sizes of each dimension (an array at least ndim long)
// -option : consists of 3 parts - flag option and an option to create a new TVirtualFFT
//         1) transform type option:
//           Available transform types are:
//           C2CForward, C2CBackward, C2R, R2C, R2HC, HC2R, DHT
//           see class description for details
//         2) flag option: choosing how much time should be spent in planning the transform:
//           Possible options:
//           "ES" (from "estimate") - no time in preparing the transform,
//                                  but probably sub-optimal  performance
//           "M"  (from "measure")  - some time spend in finding the optimal way
//                                  to do the transform
//           "P" (from "patient")   - more time spend in finding the optimal way
//                                  to do the transform
//           "EX" (from "exhaustive") - the most optimal way is found
//           This option should be chosen depending on how many transforms of the
//           same size and type are going to be done.
//           Planning is only done once, for the first transform of this size and type.
//         3) option allowing to choose between the global fgFFT and a new TVirtualFFT object
//           ""  - default, changes and returns the global fgFFT variable
//           "K" (from "keep")- without touching the global fgFFT,
//           creates and returns a new TVirtualFFT*. User is then responsible for deleting it.


////Test:

//TRootApp app("App", &argc, argv);

////No need to delete this canvas pointer. TRootApp will take care about it.
//TCanvas * canvas = app.CreateCanvas("c1", "A Simple Graph Example", 200, 10, 1600, 800);
//canvas->Divide(3);
//canvas->cd(1);

////A function to sample
//TF1 *fsin = new TF1("fsin", "sin(x)+sin(2*x)+sin(0.5*x)+1", 0, 4*TMath::Pi());
//fsin->Draw();

//const Int_t n = 25;
//TH1D * hsin = new TH1D("hsin", "hsin", n, 0, 4*TMath::Pi());

////Fill the histogram with function values
//for (Int_t i = 0; i < n; i++){
//	Double_t x = (Double_t(i - 0.5)/n)*(4*TMath::Pi());
//	hsin->SetBinContent(i, fsin->Eval(x));
//}
//hsin->Draw("same");

//canvas->cd(2);
//std::vector<Double_t> vin(static_cast<size_t>(n)),
//											re_out(static_cast<size_t>(n)),
//											im_out(static_cast<size_t>(n)),
//											back(static_cast<size_t>(n));
//for (size_t i = 0; i < static_cast<size_t>(n); i++)
//{
//	Double_t x = (static_cast<Double_t>(i-0.5)/n)*(4*TMath::Pi());
//	vin[i] = fsin->Eval(x);
//}
//TFFT fft({TFFTOperation::R2C_FORWARD, TFFTNormalization::TWO_SIDED, n});
//fft.operator()(re_out.data(), im_out.data(), vin.data());
//fft.SetFFTProperties({TFFTOperation::C2R_BACKWARD, TFFTNormalization::TWO_SIDED, n});
//fft.operator()(re_out.data(), im_out.data(), back.data());

//TH1D hist("fft", "fft", n, 0, 4*TMath::Pi());
//for (size_t i =0 ; i < static_cast<size_t>(n); i++)
//	 hist.SetBinContent(static_cast<Int_t>(i), re_out[i]);
//hist.Draw();

//canvas->cd(3);
//TH1D hist_back("back", "back", n, 0, 4*TMath::Pi());
//for (size_t i =0 ; i < static_cast<size_t>(n); i++)
//	 hist_back.SetBinContent(static_cast<Int_t>(i), back[i]);
//hist_back.Draw();

//try {
// app.Run(kTRUE);
//} catch (std::exception &e) {
// std::cerr<<e.what();
//};

enum class TFFTOperation {C2C_FORWARD, C2C_BACKWARD, R2C_FORWARD, C2R_BACKWARD, R2R_FORWARD, R2R_BACKWARD};
///1/sqrt(n) в обе стороны либо 1/n в обратную
enum class TFFTNormalization{NONE, BACKWARD, TWO_SIDED};

struct TFFTProperties {
	TFFTOperation FFTOperation;
	TFFTNormalization FFTNormalization;
	int N;
};

///Implementation of the Fourier transform
class TFFT {
private:
	TFFTProperties FFTProperties;
public:
	TFFT(const TFFTProperties &fft_properties);
	virtual ~TFFT(){}

	void SetFFTProperties(const TFFTProperties &fft_properties) noexcept;
	TFFTProperties GetFFTProperties() const noexcept;

	void operator () (const Double_t * re_in, const Double_t * im_in, Double_t * re_out, Double_t * im_out);

	void operator () (Double_t * re_c, Double_t * im_c, Double_t * r); 

  //void operator () (const Double_t * r_in, Double_t * r_out);
};	//TFFT
