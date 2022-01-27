#pragma once

#include <vector>
#include <exception>
#include <memory>

#include "TApplication.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TAxis.h"


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

static const Color_t FILL_COLOR = kWhite;
static const Color_t FRAME_FILL_COLOR = kWhite;

template <class T>
inline void ApplyAxisStyle (T * obj, const char * obj_title = "", const char * xtitle = "", const char * ytitle = "")
{
	obj->SetStats(kFALSE);
  obj->SetTitle(obj_title);
  obj->Draw("ACP");   //AL
  obj->GetXaxis()->SetLabelFont(22);
  obj->GetXaxis()->SetTitleFont(22);
	obj->GetXaxis()->SetTickLength(0.02f);
  obj->GetYaxis()->SetLabelFont(22);
  obj->GetYaxis()->SetTitleFont(22);
	obj->GetYaxis()->SetTickLength(0.02f);
	obj->GetXaxis()->SetLabelSize(0.04f);
	obj->GetXaxis()->SetLabelOffset(0.01f);
	obj->GetXaxis()->SetTitleSize(0.04f);
	obj->GetXaxis()->SetTitleOffset(1.1f);
	obj->GetYaxis()->SetLabelSize(0.04f);
	obj->GetYaxis()->SetLabelOffset(0.01f);
	obj->GetYaxis()->SetTitleSize(0.04f);
  obj->GetYaxis()->SetTitleOffset(1);
  obj->GetXaxis()->SetTitle(xtitle);
  obj->GetYaxis()->SetTitle(ytitle);
}

class TRootApp : public TApplication {
protected:
public:
  std::vector<std::unique_ptr<TCanvas>> Canvases;

  TRootApp(const char * app_class_name, Int_t * argc, char **argv)
    : TApplication(app_class_name, argc, argv){}

	virtual ~TRootApp(){}

  ///Remember that origin of coordinate system in the bottom left
  TCanvas * CreateCanvas(
    const std::string &name = "",
    const std::string &title = "",
    Int_t x = 0,
    Int_t y = 0,
    Int_t width = 1024,
    Int_t height = 768
  ){
    std::string temp = "C" + std::to_string(Canvases.size());
    TCanvas * canvas = new TCanvas(name.empty()?temp.c_str():name.c_str(), title.empty()?temp.c_str():title.c_str(), x, y, width, height);
    canvas->SetFillColor(FILL_COLOR);
    canvas->GetFrame()->SetBorderMode(0);
    canvas->SetFrameFillColor(FRAME_FILL_COLOR);
    canvas->Draw();
    Canvases.emplace_back(std::unique_ptr<TCanvas>(canvas));
    return canvas;
  }

  ///Get canvas by number
	TCanvas * GetCanvas(size_t canvas_number) const
  {
    if (canvas_number < size_t(0) || canvas_number >= Canvases.size())
      throw std::runtime_error("TRootApplication::GetCanvas error: Incorrect canvas number!");
    return Canvases[canvas_number].get();
  }

  ///Draw all the canvases
	void DrawAll() const
  {
    for (auto &it : Canvases)
      it->Draw();
  }
}; //TRootApp
