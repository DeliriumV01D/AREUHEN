#include <iostream>
#include <iomanip>

#include "TCanvas.h"
#include "TMath.h"
#include "TRootApp.h"
#include "TAcousticPressure.h"
#include "TCylindersSuperpositionApproximatedEDD.h"

#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"

class TLongitudinallyIntegratedEDDWrapper {
protected:
	TEnergyDensityDistribution * EDD;
public:
	TLongitudinallyIntegratedEDDWrapper(TEnergyDensityDistribution * edd)
	{
		EDD = edd;
	}

	Double_t operator () (const Double_t * x, const Double_t * p)
	{
		Double_t	phi = static_cast<Double_t>(x[0]),
							r = static_cast<Double_t>(x[1]);
		//Wrapper functor for the subintegral expression of the energy longitudinal distribution density
		auto longitudinally_integrable_red = [=/*&EDD*/](const Double_t * x, const Double_t * p)
		{
			fp	z = static_cast<fp>(x[0]),
					r = static_cast<fp>(p[0]),
					phi = static_cast<fp>(p[1]);
			return static_cast<Double_t>((*EDD)(r, z));	//does not depend on the variable phi
		};
		//name, func, xmin, xmax, nparameters
		TF1 f("func", &longitudinally_integrable_red, 0, 1500, 2);
		f.SetParameters(r, phi);
		return f.Integral(0, 1500);
	}
};		//TLongitudinallyIntegratedRedWrapper

class TTransverselyIntegratedEDDWrapper {
protected:
	TEnergyDensityDistribution * EDD;
public:
	TTransverselyIntegratedEDDWrapper(TEnergyDensityDistribution * edd)
	{
		EDD = edd;
	}

	Double_t operator () (const Double_t * x, const Double_t * p)
	{
		Double_t	z = static_cast<Double_t>(x[1]),
							r = static_cast<Double_t>(x[0]);
		//Wrapper functor for the subintegral expression of the energy transverse distribution density
		auto transversely_integrable_red =[=/*&EDD*/](const Double_t * x, const Double_t * p)
		{
			fp	phi = static_cast<fp>(x[0]),
					z = static_cast<fp>(p[1]),
					r = static_cast<fp>(p[0]);
			return static_cast<Double_t>((*EDD)(r, z));	//does not depend on the variable phi
		};
		//name, func, xmin, xmax, nparameters
		TF1 f("func", &transversely_integrable_red, 0, 2 * TMath::Pi(), 2);
		f.SetParameters(r, z);
		return f.Integral(0, 2 * TMath::Pi());
	}
};	//TTransverselyIntegratedREDWrapper

int main(int argc, char ** argv)
{
	TRootApp app("App", &argc, argv);
	auto canvas = app.CreateCanvas("c1", "Energy distribution", 200, 10, 1600, 800);
	canvas->Divide(4);
	canvas->Draw();

	TEnergyDensityDistributionParams params = {
		1.e18f,						///total energy of the hadronic cascade, eV
		TEDDParameterization::CORSIM,
		20.f,							///transverse cutoff threshold, g/cm2
		1500.f						///longitudinal cutoff threshold, g/cm2
	};
	std::unique_ptr<TEnergyDensityDistribution> edd = std::unique_ptr<TEnergyDensityDistribution>(new TEnergyDensityDistribution(params.E0, params.Parameterization, params.TrCutThreshold, params.LonCutThreshold));
	TLongitudinallyIntegratedEDDWrapper longitudinally_integrated_edd(edd.get());
	TTransverselyIntegratedEDDWrapper transversely_integrated_edd(edd.get());

	//name, func, xmin, xmax, ymin, ymax, nparameters
	TF2 lf("longitudinally integrated edd", &longitudinally_integrated_edd, 0, 2 * TMath::Pi(), 0, 10,  0);
	lf.SetParameters(0., 0.); //
	app.GetCanvas(0)->cd(1);
	lf.Draw("SURF7 POL");//("POLAR SURF3"); //pollego2z

	//name, func, xmin, xmax, ymin, ymax, nparameters
	TF2 trf("transversely integrated edd", &transversely_integrated_edd, 0, 10, 0, 1500, 0);
	trf.SetParameters(0., 0.); //
	app.GetCanvas(0)->cd(2);
	trf.Draw("SURF7");

	std::unique_ptr<TCylindersSuperpositionApproximatedEDD> csaedd = std::unique_ptr<TCylindersSuperpositionApproximatedEDD>(new TCylindersSuperpositionApproximatedEDD(params, 25));
	TLongitudinallyIntegratedEDDWrapper longitudinally_integrated_csaedd(csaedd.get());
	TTransverselyIntegratedEDDWrapper transversely_integrated_csaedd(csaedd.get());

	//name, func, xmin, xmax, ymin, ymax, nparameters
	TF2 lf2("longitudinally integrated cylinders superposition approximated edd", &longitudinally_integrated_csaedd, 0, 2 * TMath::Pi(), 0, 10,  0);
	lf2.SetParameters(0., 0.); //
	app.GetCanvas(0)->cd(3);
	lf2.Draw("SURF7 POL");//("POLAR SURF3"); //pollego2z

	//name, func, xmin, xmax, ymin, ymax, nparameters
	TF2 trf2("transversely integrated cylinders superposition approximated edd", &transversely_integrated_csaedd, 0, 10, 0, 1500, 0);
	trf2.SetParameters(0., 0.); //
	app.GetCanvas(0)->cd(4);
	trf2.Draw("SURF7");

	try {
	 app.Run(kTRUE);
	} catch (std::exception &e) {
	 std::cerr<<e.what();
	};
}
