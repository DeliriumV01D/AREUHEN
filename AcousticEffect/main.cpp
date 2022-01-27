
#include "TCanvas.h"
#include "TMath.h"
#include "TRootApp.h"
#include "TAcousticPressure.h"

#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"

int main(int argc, char ** argv)
{
	TRootApp app("App", &argc, argv);
	auto canvas = app.CreateCanvas("c1", "Energy density distribution", 200, 10, 1600, 800);
	canvas->Divide(2);
	canvas->Draw();
	canvas = app.CreateCanvas("c2", "Acoustic pressure", 200, 10, 1600, 800);
	canvas->Divide(2);
	canvas->Draw();

	const fp	rd = 100000.,
						zd = 350;
	fp fmax = 1.e5;
	const size_t n = 128;
	std::vector<fp> mpwx(n),
									pt(n);
	std::vector<std::complex<fp>> pwx(n);

	TAcousticPressureParams params;
	params.E0 = 1.e20;				//total energy of the hadronic cascade, eV
	params.cs = 145000.;			//speed of sound in a medium, g/cm2/s
	params.Rd = rd;						//detector coordinates, g/cm2
	params.Zd = zd;						//GetLMax(), detector coordinates, g/cm2
	params.LonCutThreshold = 1000.;		//longitudinal cutoff threshold, g/cm2
	params.TrCutThreshold = 20.;			//transverse cutoff threshold, g/cm2
	params.Parameterization = TEDDParameterization::CORSIM;

	//Monte Carlo method for the calculation of the acoustic pressure of the hadronic cascade
	//test of unuran passing as input a multi-dimension distribution object
	TMCICAP acoustic_pressure(1e6);
	acoustic_pressure.SetConditions(params);

	TH3D h1(	"h13D","Energy density distribution",300,
						-0.5*params.TrCutThreshold,
						0.5*params.TrCutThreshold,
						300,
						-0.5*params.TrCutThreshold,
						0.5*params.TrCutThreshold,
						300,
						0,
						params.LonCutThreshold
					);
	TH2D h12(  "h12D","Energy density distribution",
						 100,
						 -0.5*params.TrCutThreshold,
						 0.5*params.TrCutThreshold,
						 150,
						 0,
						 params.LonCutThreshold
					);

	//Acoustic pressure depending on sound frequency (frequency spectrum)
	acoustic_pressure.GetComplexPwSeries(rd, zd, fmax, n, pwx, mpwx, &h1, &h12);

	app.GetCanvas(0)->cd(1);
	ApplyAxisStyle (&h1, "Energy density distribution", "x, g/cm^{2}", "");
	h1.GetZaxis()->SetLabelFont(22);
	h1.GetZaxis()->SetTitleFont(22);
	h1.GetZaxis()->SetTickLength(0.02f);
	h1.GetZaxis()->SetLabelSize(0.04f);
	h1.GetZaxis()->SetLabelOffset(0.01f);
	h1.GetZaxis()->SetTitleSize(0.04f);
	h1.GetZaxis()->SetTitleOffset(1.1f);
	h1.GetZaxis()->SetTitle("z, g/cm^{2}");
	h1.GetYaxis()->SetLabelSize(0);
	h1.Draw();

	app.GetCanvas(0)->cd(2);
	ApplyAxisStyle (&h12, "Energy density distribution", "x, g/cm^{2}", "z, g/cm^{2}");
	h12.Smooth(1);
	h12.SetContour(11);
	h12.Draw("CONT3");

	TH1D	Pw("Pw", "Pw-title", n, 0, fp(n-1)*fmax/n/1000),
				Pt("Pt", "Pt-title", n, 0, fp(n-1)/fmax);

	for (size_t i = 0; i < n; i++ )
		Pw.SetBinContent(static_cast<Int_t>(i), mpwx[i] * abs(/*(*(PRZt->Gamma))(4,0,isBaikal)*/.1/6.25/1e18));	// /6.25/1E18 ->Pa
	Pw.SetBins(n-1, 0, fp(n-1)*fmax/n/1000);		//Hz->kHz

	//Acoustic pressure depending on time
	acoustic_pressure.GetPtSeries(pwx, pt, fmax, n);

	for (size_t i = 0; i < n; i++ )
		Pt.SetBinContent(static_cast<Int_t>(i), pt[i] * /*(*(PRZt->Gamma))(4,0,isBaikal)*/.1/6.25/1e18);// /6.25/1E18 -> Pa
	Pt.SetBins(n-1, 0, fp(n-1)/fmax);
	app.GetCanvas(1)->cd(1);
	ApplyAxisStyle(&Pt,"P(t)", "t,sec", "P,Pa");
	Pt.Draw("L");
	app.GetCanvas(1)->cd(2);
	ApplyAxisStyle(&Pw, "|P(f)|","f,kHz","|P(f)|, Hz#upointkg#upointm^{-1}");
	Pw.Draw("C");

	try {
	 app.Run(kTRUE);
	} catch (std::exception &e) {
	 std::cerr<<e.what();
	};
}
