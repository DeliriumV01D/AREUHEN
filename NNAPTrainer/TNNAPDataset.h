#pragma once

#include "TorchHeader.h"

#include <iostream>
#include <tuple>
#include <vector>

#include "CommonDefinitions.h"
#include "TRootApp.h"
#include "TAcousticPressure.h"
#include "TCylindersSuperpositionApproximatedEDD.h"
#include "TRandomInt.h"
#include "TRandomDouble.h"


struct TNNAPDatasetParams {
	fp	fmax0 = static_cast<fp>(1.e5),
			fmax1 = static_cast<fp>(2.e5),
			E0 = static_cast<fp>(1.e16),				//total energy of the hadronic cascade, eV
			E1 = static_cast<fp>(1.e20 - 1),		//total energy of the hadronic cascade, eV
			cs0 = 130000,												//speed of sound in a medium, g/cm2/s
			cs1 = 180000,												//speed of sound in a medium, g/cm2/s
			Rd0 = 5000,								//detector coordinates, g/cm2
			Rd1 = 1000000,						//detector coordinates, g/cm2
			Zd0 = -3000,									//GetLMax(), detector coordinates, g/cm2
			Zd1 = 3000,									//GetLMax(), detector coordinates, g/cm2
			LonCutThreshold = static_cast<fp>(1000),		//longitudinal cutoff threshold, g/cm2
			TrCutThreshold = static_cast<fp>(20);				//transverse cutoff threshold, g/cm2
	int N0 = 128,		//Calculate from fmax?
			N1 = 256,
			DataSize = 100000;
	TEDDParameterization Parameterization = TEDDParameterization::CORSIM;
};

class TNNAPDataset : public torch::data::datasets::Dataset<TNNAPDataset> {
	using Example = torch::data::Example<>;
protected:
	torch::Tensor Data;
	torch::Tensor Labels;
	TNNAPDatasetParams Params;

	std::tuple<torch::Tensor, torch::Tensor> GenerateData(const long sz);
public:
	TNNAPDataset(const TNNAPDatasetParams &params);
	virtual ~TNNAPDataset() override {}
	Example get(size_t index) override;
	torch::optional<size_t> size() const override;
};
