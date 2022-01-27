#pragma once

#include <cmath>
#include <TMath.h>
#include "CommonDefinitions.h"

enum class TEDDParameterization {CORSIM/*, SAUND, NIESSBERTIN*/};

///Heat energy density distribution depending on hadronic cascade energy
/**
CORSIM:
S. Bevan, S. Danaher, J. Perkin, S. Ralph, C. Rhodes, L. Thompson, T. Sloan, D. Waters.
Simulation of Ultra High Energy Neutrino Interactions in Ice and Water.
Astroparticle Physics, Volume 28, Issue 3, p. 366-379. 2007. ELSEVIER. DOI:10.1016/j.astropartphys.2007.08.001. arXiv:0704.1025.
*/
class TEnergyDensityDistribution {
protected:
  TEDDParameterization Parameterization;
  fp  E0,
      El,
      TrCutThreshold,
      LonCutThreshold;

	///CORSIM CONSTS
  fp P1L, P2L, P3L, P4L, P5L, P6L, A1, B1, C1, A2, B2, C2;

	fp L_CORSIM (const fp &z) const
  {
		return P1L * pow((z - P2L)/(P3L - P2L), ((P3L - P2L)/(P4L + P5L*z + P6L * pow(z, fp(2.))))) * exp((P3L - z)/(P4L + P5L*z + P6L * pow(z, fp(2.))));
  }

	fp R_CORSIM (const fp &r, const fp &z) const
  {
		const fp P1R = A1 + B1 * z + C1 * pow(z, fp(2.));
		const fp P2R = A2 + B2 * z + C2 * pow(z, fp(2.));

		const fp I = P1R * TMath::Gamma(4.5 - 2 * P2R) * TMath::Gamma(P2R) / TMath::Gamma(4.5 - P2R);
		return fp(1) / I * pow((r / P1R), (P2R - 1)) * pow((1 + r / P1R), (P2R - 4.5));
  }
public:
  TEnergyDensityDistribution(
		const fp E0,                                    ///total energy of the hadronic cascade, eV
    const TEDDParameterization parameterization,
		const fp &tr_cut_threshold,                     ///transverse cutoff threshold, g/cm2
		const fp &lon_cut_threshold                     ///longitudinal cutoff threshold, g/cm2
  ) : TEnergyDensityDistribution(parameterization, tr_cut_threshold, lon_cut_threshold)
  {
    SetE0(E0);
  }

  TEnergyDensityDistribution(
    const TEDDParameterization parameterization,
		const fp &tr_cut_threshold,                     ///transverse cutoff threshold, g/cm2
		const fp &lon_cut_threshold                     ///longitudinal cutoff threshold, g/cm2
  ){
    Parameterization = parameterization;
    TrCutThreshold = tr_cut_threshold;
    LonCutThreshold = lon_cut_threshold;
  }

  virtual ~TEnergyDensityDistribution(){}

  void SetE0(const fp &E0)
  {
    this->E0 = E0;
    El = log10(E0);
    if (Parameterization == TEDDParameterization::CORSIM)
    {
			fp El2 = pow(El, static_cast<fp>(2));
			P1L = E0*(2.76e-3 - 1.974e-4 * El + 7.45e-6 * El2);
			P2L = - 210.9 - 6.968e-3 * El + 0.1551 * El2;
			P3L = - 41.50 + 113.9 * El - 4.103 * El2;
			P4L = 8.012 + 11.44 * El - 0.5434 * pow(El, fp(2.));
			P5L = 0.7999 * 1e-5 - 0.004843 * El + 0.0002552 * El2;
			P6L = 4.563e-5 - 3.504e-6 * El + 1.315e-7 * El2;

			A1 = 0.01287 * El2 - 0.2573 * El + 0.9636;
			B1 = -0.4697 * 1e-4 * El2 + 0.0008072 * El + 0.0005404;
			C1 = 0.7344e-7 * El2 - 1.375e-6 * El + 4.488e-6;

			A2 = -0.8905e-3 * El2 + 0.007727 * El + 1.969;
			B2 = 0.1173e-4 * El2 - 0.0001782 * El - 5.093e-6;
			C2 = -0.1058e-7 * El2 + 0.1524e-6 * El - 0.1069e-8;
    }
  }

	fp GetE0() const noexcept { return E0; }
	fp GetEl() const noexcept { return El; }
	TEDDParameterization GetParametrization() const noexcept { return Parameterization; }

	///eps = 1/(E0*2*Pi()*r) * d^2E/drdz //Relative energy density
	///Cylindrical coordinate system
  fp operator () (const fp r, const fp z)
  {
    fp r_int = r;
    if (Parameterization == TEDDParameterization::CORSIM && r_int < 10e-6)
			r_int = 10e-6;
		fp result = L(z) * R(r_int, z) / (GetE0() * 2 * TMath::Pi() * r_int);
    if (!isfinite(result))
      return 0.;
		else
      return result;
  }

	///eps = 1/(E0*2*Pi()*r) * d^2E/drdz //Relative energy density
	///Spherical coordinate system, format required by ROOT double r = sqrt(pow(x[0],2)+pow(x[1],2)); double z = x[2];
  fp operator () (const double * x, const double * p)
  {
    const double  r = sqrt(pow(x[0], 2) + pow(x[1], 2)),
                  z = x[2];
    double result;
    if (z >= 0 && z < LonCutThreshold && r >= 0 && r < TrCutThreshold)
      result = (*this)(r,z);
    else
      result = 0.;
    return result;
  }

  ///Longitudinal distribution of the energy density
  fp L (const fp &z)
  {
    switch ( Parameterization )
    {
      case TEDDParameterization::CORSIM:
        return L_CORSIM(z);
        break;
        //case SAUND:
        //  return L_SAUND(z);
        //  break;
        //case NIESSBERTIN:
        //  return L_NIESSBERTIN(z);
        //  break;
    };
  }

  ///Transverse distribution of the energy density
  fp R (const fp r, const fp z)
  {
    switch ( Parameterization )
    {
      case TEDDParameterization::CORSIM:
        return R_CORSIM(r, z);
        break;
        //case SAUND:
        //  return R_SAUND(r, z);
        //  break;
        //case NIESSBERTIN:
        //  return R_NIESSBERTIN(r, z);
        //  break;
    };
		return 0;
  }
};
