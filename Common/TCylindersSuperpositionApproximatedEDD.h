#pragma once

#include <memory>

#include "TEnergyDistributionDensity.h"

//Test:
//
//#include <iostream>
//#include <iomanip>

//#include "TCanvas.h"
//#include "TMath.h"
//#include "TRootApp.h"
//#include "TAcousticPressure.h"
//#include "TCylindersSuperpositionApproximatedEDD.h"

//#include "TH2D.h"
//#include "TH3D.h"
//#include "TMath.h"
//#include "TF1.h"
//#include "TF2.h"

//class TLongitudinallyIntegratedEDDWrapper {
//protected:
//	TEnergyDensityDistribution * EDD;
//public:
//	TLongitudinallyIntegratedEDDWrapper(TEnergyDensityDistribution * edd)
//	{
//		EDD = edd;
//	}

//	Double_t operator () (const Double_t * x, const Double_t * p)
//	{
//		Double_t	phi = static_cast<Double_t>(x[0]),
//							r = static_cast<Double_t>(x[1]);
//		//Wrapper functor for the subintegral expression of the energy longitudinal distribution density
//		auto longitudinally_integrable_red = [=/*&EDD*/](const Double_t * x, const Double_t * p)
//		{
//			fp	z = static_cast<fp>(x[0]),
//					r = static_cast<fp>(p[0]),
//					phi = static_cast<fp>(p[1]);
//			return static_cast<Double_t>((*EDD)(r, z));	//does not depend on the variable phi
//		};
//		//name, func, xmin, xmax, nparameters
//		TF1 f("func", &longitudinally_integrable_red, 0, 1500, 2);
//		f.SetParameters(r, phi);
//		return f.Integral(0, 1500);
//	}
//};		//TLongitudinallyIntegratedRedWrapper

//class TTransverselyIntegratedEDDWrapper {
//protected:
//	TEnergyDensityDistribution * EDD;
//public:
//	TTransverselyIntegratedEDDWrapper(TEnergyDensityDistribution * edd)
//	{
//		EDD = edd;
//	}

//	Double_t operator () (const Double_t * x, const Double_t * p)
//	{
//		Double_t	z = static_cast<Double_t>(x[1]),
//							r = static_cast<Double_t>(x[0]);
//		//Wrapper functor for the subintegral expression of the energy transverse distribution density
//		auto transversely_integrable_red =[=/*&EDD*/](const Double_t * x, const Double_t * p)
//		{
//			fp	phi = static_cast<fp>(x[0]),
//					z = static_cast<fp>(p[1]),
//					r = static_cast<fp>(p[0]);
//			return static_cast<Double_t>((*EDD)(r, z));	//does not depend on the variable phi
//		};
//		//name, func, xmin, xmax, nparameters
//		TF1 f("func", &transversely_integrable_red, 0, 2 * TMath::Pi(), 2);
//		f.SetParameters(r, z);
//		return f.Integral(0, 2 * TMath::Pi());
//	}
//};	//TTransverselyIntegratedREDWrapper

//int main(int argc, char ** argv)
//{
//	TRootApp app("App", &argc, argv);
//	auto canvas = app.CreateCanvas("c1", "Energy distribution", 200, 10, 1600, 800);
//	canvas->Divide(4);
//	canvas->Draw();

//	TEnergyDensityDistributionParams params = {
//		1.e18f,						///total energy of the hadronic cascade, eV
//		TEDDParameterization::CORSIM,
//		20.f,							///transverse cutoff threshold, g/cm2
//		1500.f						///longitudinal cutoff threshold, g/cm2
//	};
//	std::unique_ptr<TEnergyDensityDistribution> edd = std::unique_ptr<TEnergyDensityDistribution>(new TEnergyDensityDistribution(params.E0, params.Parameterization, params.TrCutThreshold, params.LonCutThreshold));
//	TLongitudinallyIntegratedEDDWrapper longitudinally_integrated_edd(edd.get());
//	TTransverselyIntegratedEDDWrapper transversely_integrated_edd(edd.get());

//	//name, func, xmin, xmax, ymin, ymax, nparameters
//	TF2 lf("longitudinally integrated edd", &longitudinally_integrated_edd, 0, 2 * TMath::Pi(), 0, 10,  0);
//	lf.SetParameters(0., 0.); //
//	app.GetCanvas(0)->cd(1);
//	lf.Draw("SURF7 POL");//("POLAR SURF3"); //pollego2z

//	//name, func, xmin, xmax, ymin, ymax, nparameters
//	TF2 trf("transversely integrated edd", &transversely_integrated_edd, 0, 10, 0, 1500, 0);
//	trf.SetParameters(0., 0.); //
//	app.GetCanvas(0)->cd(2);
//	trf.Draw("SURF7");

//	std::unique_ptr<TCylindersSuperpositionApproximatedEDD> csaedd = std::unique_ptr<TCylindersSuperpositionApproximatedEDD>(new TCylindersSuperpositionApproximatedEDD(params, 25));
//	TLongitudinallyIntegratedEDDWrapper longitudinally_integrated_csaedd(csaedd.get());
//	TTransverselyIntegratedEDDWrapper transversely_integrated_csaedd(csaedd.get());

//	//name, func, xmin, xmax, ymin, ymax, nparameters
//	TF2 lf2("longitudinally integrated cylinders superposition approximated edd", &longitudinally_integrated_csaedd, 0, 2 * TMath::Pi(), 0, 10,  0);
//	lf2.SetParameters(0., 0.); //
//	app.GetCanvas(0)->cd(3);
//	lf2.Draw("SURF7 POL");//("POLAR SURF3"); //pollego2z

//	//name, func, xmin, xmax, ymin, ymax, nparameters
//	TF2 trf2("transversely integrated cylinders superposition approximated edd", &transversely_integrated_csaedd, 0, 10, 0, 1500, 0);
//	trf2.SetParameters(0., 0.); //
//	app.GetCanvas(0)->cd(4);
//	trf2.Draw("SURF7");

//	try {
//	 app.Run(kTRUE);
//	} catch (std::exception &e) {
//	 std::cerr<<e.what();
//	};
//}

struct TCyl {
	fp	AverageED,
			BorderED,
			Zn,
			Zp,
			R;
	TCyl() : AverageED(0), BorderED(0), Zn(0), Zp(0), R(0){}
};

///Energy density distribution is approximated by a superposition of cylinders with a uniform distribution density
class TCylindersSuperpositionApproximatedEDD : public TEnergyDensityDistribution {
protected:
	unsigned int NCyl;
	std::vector<TCyl> Cylinders;
public:
	TCylindersSuperpositionApproximatedEDD(const TEnergyDensityDistributionParams &params, const unsigned int ncyl)
		: TEnergyDensityDistribution(params.E0, params.Parameterization, params.TrCutThreshold, params.LonCutThreshold)
	{
		NCyl = ncyl;
		Cylinders = std::vector<TCyl>(static_cast<size_t>(NCyl));

		const unsigned int nsteps = 100 * NCyl;
		fp	r0 = static_cast<fp>(0.15),
				lmax = this->EvaluateLMax(),
				max_ed = TEnergyDensityDistribution::operator()(0, lmax),
				h_ed = (max_ed - max_ed * static_cast<fp>(0.001))/(NCyl + 1),
				h_r = params.TrCutThreshold/nsteps,
				h_z_n = lmax/nsteps,
				h_z_p = (params.LonCutThreshold - lmax)/nsteps;

		//Initializing cylinders
		for (size_t i = 0; i < NCyl; i++)
		{
			Cylinders[i].BorderED = max_ed - h_ed * (static_cast<fp>(i) + 1);
			Cylinders[i].AverageED = max_ed - h_ed * (static_cast<fp>(i) + static_cast<fp>(0.5));
			Cylinders[i].R = 0;
			Cylinders[i].Zn = lmax;
			Cylinders[i].Zp = lmax;
		}

		//if we go from outside to inside then we should choose > sign to not overwrite the correct values of the cylinders boundaries
		for (unsigned int step = 0; step < nsteps; step++)
		{
			fp	zn = lmax - h_z_n * step,
					zp = lmax + h_z_p * step,
					r = h_r * step,
					zned = TEnergyDensityDistribution::operator()(r0, zn),
					zped = TEnergyDensityDistribution::operator()(r0, zp),
					red = TEnergyDensityDistribution::operator()(r, lmax);

			//in the descending direction
			for (unsigned int i = 0; i < NCyl; i++)
			{
				if (zned >= Cylinders[i].BorderED && zn < Cylinders[i].Zn)
					Cylinders[i].Zn = zn;
				if (zped >= Cylinders[i].BorderED && zp > Cylinders[i].Zp)
					Cylinders[i].Zp = zp;
				if (red >= Cylinders[i].BorderED && r > Cylinders[i].R)
					Cylinders[i].R = r;
			}
		}

		//Don't forget about negative-ed cylinders to avoid multiple accounting in the superposition
		for (unsigned int i = 0; i < NCyl-1; i++)
			Cylinders[i].AverageED -= Cylinders[i+1].AverageED;
	}

	virtual ~TCylindersSuperpositionApproximatedEDD() override {}

	virtual fp operator () (const fp r, const fp z) override
	{
		fp result = 0;
		//Calculate the superposition of the cylinder energy densities in the given point
		for (unsigned int i = 0; i < NCyl; i++)
			if (r <= Cylinders[i].R && z >= Cylinders[i].Zn && z <= Cylinders[i].Zp)
				result += Cylinders[i].AverageED;
		return result;
	}

	std::vector<TCyl> * GetCylinders(){ return &Cylinders; }
};
