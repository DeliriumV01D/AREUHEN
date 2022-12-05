#pragma once

#include <complex>
#include <memory>

#include "TH2.h"
#include "TH3.h"

#include "TEnergyDistributionDensity.h"
#include "TFilonIntegrator.h"
#include "TCylindersSuperpositionApproximatedEDD.h"

//Simple test:
//
//#include "TCanvas.h"
//#include "TMath.h"
//#include "TRootApp.h"
//#include "TAcousticPressure.h"

//#include "TH2D.h"
//#include "TH3D.h"
//#include "TF1.h"

//int main(int argc, char ** argv)
//{
//	TRootApp app("App", &argc, argv);
//	auto canvas = app.CreateCanvas("c1", "Energy density distribution", 200, 10, 1600, 800);
//	canvas->Divide(2);
//	canvas->Draw();
//	canvas = app.CreateCanvas("c2", "Acoustic pressure", 200, 10, 1600, 800);
//	canvas->Divide(2);
//	canvas->Draw();

//	const fp	rd = 100000,
//						zd = 350;
//	fp fmax = static_cast<fp>(1.e5);
//	const size_t n = 128;
//	std::vector<fp> mpwx(n),
//									pt(n);
//	std::vector<std::complex<fp>> pwx(n);

//	TAcousticPressureParams params;
//	params.E0 = static_cast<fp>(1.e20);				//total energy of the hadronic cascade, eV
//	params.cs = static_cast<fp>(145000);			//speed of sound in a medium, g/cm2/s
//	params.Rd = rd;						//detector coordinates, g/cm2
//	params.Zd = zd;						//GetLMax(), detector coordinates, g/cm2
//	params.LonCutThreshold = static_cast<fp>(1000);		//longitudinal cutoff threshold, g/cm2
//	params.TrCutThreshold = static_cast<fp>(20);			//transverse cutoff threshold, g/cm2
//	params.Parameterization = TEDDParameterization::CORSIM;

//	//Monte Carlo method for the calculation of the acoustic pressure of the hadronic cascade
//	//test of unuran passing as input a multi-dimension distribution object
//	TMCICAP acoustic_pressure(1e6);

//	acoustic_pressure.SetConditions(params);

//	TH3D h1(	"h13D","Energy density distribution",300,
//						Double_t(-0.5 * Double_t(params.TrCutThreshold)),
//						Double_t(0.5 * Double_t(params.TrCutThreshold)),
//						300,
//						Double_t(-0.5 * Double_t(params.TrCutThreshold)),
//						Double_t(0.5 * Double_t(params.TrCutThreshold)),
//						300,
//						0,
//						static_cast<Double_t>(params.LonCutThreshold)
//					);
//	TH2D h12(  "h12D","Energy density distribution",
//						 100,
//						 Double_t(-0.5 * Double_t(params.TrCutThreshold)),
//						 0.5 * Double_t(params.TrCutThreshold),
//						 150,
//						 0,
//						 Double_t(params.LonCutThreshold)
//					);

//	//Acoustic pressure depending on sound frequency (frequency spectrum)
//	acoustic_pressure.GetComplexPwSeries(rd, zd, fmax, n, pwx, mpwx, &h1, &h12);

//	app.GetCanvas(0)->cd(1);
//	ApplyAxisStyle (&h1, "Energy density distribution", "x, g/cm^{2}", "");
//	h1.GetZaxis()->SetLabelFont(22);
//	h1.GetZaxis()->SetTitleFont(22);
//	h1.GetZaxis()->SetTickLength(0.02f);
//	h1.GetZaxis()->SetLabelSize(0.04f);
//	h1.GetZaxis()->SetLabelOffset(0.01f);
//	h1.GetZaxis()->SetTitleSize(0.04f);
//	h1.GetZaxis()->SetTitleOffset(1.1f);
//	h1.GetZaxis()->SetTitle("z, g/cm^{2}");
//	h1.GetYaxis()->SetLabelSize(0);
//	h1.Draw();

//	app.GetCanvas(0)->cd(2);
//	ApplyAxisStyle (&h12, "Energy density distribution", "x, g/cm^{2}", "z, g/cm^{2}");
//	h12.Smooth(1);
//	h12.SetContour(11);
//	h12.Draw("CONT3");

//	TH1D	Pw("Pw", "Pw-title", n, 0, Double_t(n-1) * Double_t(fmax) / n / 1000),
//				Pt("Pt", "Pt-title", n, 0, Double_t(n-1) / Double_t(fmax));

//	for (size_t i = 0; i < n; i++ )
//		Pw.SetBinContent(static_cast<Int_t>(i), mpwx[i] * abs(/*(*(PRZt->Gamma))(4,0,isBaikal)*/.1/6.25/1e18));	// /6.25/1E18 ->Pa
//	Pw.SetBins(n-1, 0, static_cast<Double_t>(n-1) * static_cast<Double_t>(fmax) / n / 1000);		//Hz->kHz

//	//Acoustic pressure depending on time
//	acoustic_pressure.GetPtSeries(pwx, pt, fmax, n);

//	for (size_t i = 0; i < n; i++ )
//		Pt.SetBinContent(static_cast<Int_t>(i), pt[i] * /*(*(PRZt->Gamma))(4,0,isBaikal)*/.1/6.25/1e18);// /6.25/1E18 -> Pa
//	Pt.SetBins(n-1, 0, static_cast<Double_t>(n-1)/static_cast<Double_t>(fmax));
//	app.GetCanvas(1)->cd(1);
//	ApplyAxisStyle(&Pt,"P(t)", "t,sec", "P,Pa");
//	Pt.Draw("L");
//	app.GetCanvas(1)->cd(2);
//	ApplyAxisStyle(&Pw, "|P(f)|","f,kHz","|P(f)|, Hz#upointkg#upointm^{-1}");
//	Pw.Draw("C");

//	try {
//	 app.Run(kTRUE);
//	} catch (std::exception &e) {
//	 std::cerr<<e.what();
//	};
//}

enum class TPartOfComplexNumber{Im, Re};

struct TAcousticPressureParams {
  ///own params
	TPartOfComplexNumber PartOfComplexNumber;
	fp Rd; fp Zd;									///detector coordinates, g/cm2
	fp w;													///cyclic frequency
	fp cs;												/////speed of sound in a medium, g/cm2/s
  ///TEnergyDensityDistribution delegate params
	fp E0;												///total energy of the hadronic cascade, eV
  TEDDParameterization Parameterization;
	fp TrCutThreshold;            ///transverse cutoff threshold, g/cm2
	fp LonCutThreshold;           ///longitudinal cutoff threshold, g/cm2
};

///Wrapper functor for the subintegral expression of the Kirchhoff integral depending on detector coordinates,
///sound speed, frequency (Fourier spectrum) and energy distribution density
class TFwWrapper {
private:
  TAcousticPressureParams Params;
  std::unique_ptr<TEnergyDensityDistribution> EDD;
public:
  TFwWrapper(const TAcousticPressureParams &params);
	virtual ~TFwWrapper(){}

  void SetParams(const TAcousticPressureParams &params);
	TAcousticPressureParams GetParams() const noexcept { return Params; }

  fp operator () (const fp phi, const fp r, const fp z) const;
  ///Root 		r = *arg;	z = *(arg + 1);	phi = *(arg + 2);
  Double_t operator () ( const Double_t * arg);

  TEnergyDensityDistribution * GetEDD() const noexcept;
};      //TFwWrapper


enum class TAPCalculationMethod { ROOT, ACORNE, ASKARYAN };

///Base class for calculating the acoustic pressure generated by the hadronic cascade
///with energy density distribution
class TAcousticPressure {
protected:
  std::unique_ptr<TFwWrapper> FwWrapper;
  fp LMax;
public:
	TAcousticPressure() : LMax(0) {}
  virtual ~TAcousticPressure(){}

  ///Set experiment conditions
  virtual void SetConditions(const TAcousticPressureParams &conditions);

  ///Returns the linear coordinate of the maximum of the cascade
	virtual fp GetLMax() const noexcept { return LMax; }

  ///Акустический эффект (на данной частоте w) от ЯЭК в точке распложения детектора (Rd,Zd)
  virtual std::complex<fp> GetComplexPwValue (const fp &rd, const fp &zd, const fp &w) = 0;

  ///Ряд значений из N точек акустического эффекта (на данной частоте w) развития ЯЭК в точке распложения детектора (Rd,Zd)
  ///Спектр акустического импульса, ограниченный частотой fmax
  virtual void GetComplexPwSeries(	const fp &rd,
                                    const fp &zd,
																		fp &fmax,															///Hz
                                    size_t N,
																		std::vector<std::complex<fp>> &pwx,		///Комплексный спектр
																		std::vector<fp> &mpwx                 ///Модуль спектра
                                  );
	///Acoustic pressure depending on time
	virtual void GetPtSeries(	std::vector<std::complex<fp>> &pwx,
                            std::vector<fp> &pt,
														fp &fmax,																			///Hz
                            size_t N
                          );
};     //TAcousticPressure


///Calculation of the acoustic pressure by the integration method implemented in Root
///By default uses adaptive multi-dimensional integration using the algorithm from Genz Mallik implemented in the class ROOT::Math::AdaptiveIntegratorMultiDim
class TRICAP : public TAcousticPressure {
protected:
public:
	virtual std::complex<fp> GetComplexPwValue (const fp &rd, const fp &zd, const fp &w) override;
};		//TRICAP



///Monte Carlo method
/// S. Bevan, S. Danaher, J. Perkin, S. Ralph, C. Rhodes, L. Thompson, T. Sloan, D. Waters.
/// Simulation of Ultra High Energy Neutrino Interactions in Ice and Water.
/// Astroparticle Physics, Volume 28, Issue 3, p. 366-379. 2007. ELSEVIER.
/// DOI:10.1016/j.astropartphys.2007.08.001. arXiv:0704.1025.
class TMCICAP : public TAcousticPressure {
private:
  int NPoints;
public:
  TMCICAP(const int npoints = 300000);

  virtual std::complex<fp> GetComplexPwValue (const fp &rd, const fp &zd, const fp &w) override;

	virtual void GetComplexPwSeries(	const fp &rd,
																		const fp &zd,
																		fp &fmax,															///Hz
																		size_t N,
																		std::vector<std::complex<fp>> &pwx,		///Комплексный спектр
																		std::vector<fp> &mpwx									///Модуль спектра
																		) override {GetComplexPwSeries(rd, zd, fmax, N, pwx, mpwx, nullptr, nullptr);}

	virtual void GetComplexPwSeries(	const fp &rd,
																		const fp &zd,
																		fp &fmax,															///Hz
																		size_t N,
																		std::vector<std::complex<fp>> &pwx,		///Комплексный спектр
																		std::vector<fp> &mpwx,	              ///Модуль спектра
																		TH3D * edd3d,													///Debug options
																		TH2D * edd2d													///Debug options
																	);
};		//TMCICAP




/******************************************************************************
TAICP - Метод интегрирования Аскарьяна, основанный на формуле Филона
Аскарьян Г.А., Долгошеин Б.А., Акустическая регистрация нейтрино высоких энергий.
//Письма в ЖЭТФ.-1977.-Т.25.-С.232-233.
******************************************************************************/

//static const int N_INTERVALS = 16;
//static const unsigned int N_CYLINDERS = 25;  ///!!!include to class params

/////Метод интегрирования Аскарьяна, основанный на формуле Филона decomposition into cylinders
//class TAICAP : public TAcousticPressure {
//private:
//	TAcousticPressureParams Conditions;
//	TFilonIntegrator <fp> FilonIntegrator;
//	std::unique_ptr<TCylindersSuperpositionApproximatedEDD> CSAEDD;
//public:
//	///Set experiment conditions
//	virtual void SetConditions(const TAcousticPressureParams &conditions) override;
//	virtual std::complex<fp> GetComplexPwValue (const fp &rd, const fp &zd, const fp &w) override;
//};		//TAICAP

//static const int N_INTERVALS = 16;
//static const int N_CYLINDERS = 32;	//Линейно влияет на скорость вычислений.

/////Вектор коэффициентов разложения плотности энерговыделения в суперпозицию цилиндров по оси z
/////Логарифм радиальной компоненты плотноcти внутри каждого цилиндра приближен прямой
//struct TCyl {
//	fp zi;
//	fp Ai;
//	fp lambdai;
//};

////Функция, вычисляющая коэффициенты разложения
//std::vector<TCyl> GetExpCoefficients(TEnergyDensityDistribution * edd, int N = N_CYLINDERS);

//inline fp phi_z(const fp z, const std::vector<TCyl> exp_coefficients, const fp Rd, const fp Zd, const fp w, const fp cs);

////Подынтегральная функция
//class TF {
//public:
//	fp	Rd,
//			Zd,
//			w,
//			cs;
//	bool v;
//	std::vector<TCyl> ExpCoefficients; //expansion coefficients
//	std::complex<fp> operator () (const fp &x);
//};

////Метод интегрирования Аскарьяна, основанный на формуле Филона decomposition into cylinders
//class TAICAP : public TAcousticPressure {
//private:
//	TFilonIntegrator FilonIntegrator;
//	std::vector<TCyl> ExpCoefficients;
//public:
//	///Set experiment conditions
//	virtual void SetConditions(const TAcousticPressureParams &conditions) override;

//	virtual std::complex<fp> GetComplexPwValue (const fp &rd, const fp &zd, const fp &w) override;

//	virtual void GetComplexPwSeries(	const fp &rd,
//																		const fp &zd,
//																		fp &fmax,															///Hz
//																		size_t N,
//																		std::vector<std::complex<fp>> &pwx,		///Комплексный спектр
//																		std::vector<fp> &mpwx                 ///Модуль спектра
//																	) override;
//};		//TAICAP
