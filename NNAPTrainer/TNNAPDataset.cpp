#include "TNNAPDataset.h"

std::tuple<torch::Tensor, torch::Tensor> TNNAPDataset :: GenerateData(const long sz)
{
	const int data_dim = 5,		//E, cs, Rd, Zd, f/w
						label_dim = 2;	//RePw, ImPw
	torch::Tensor data = torch::zeros({ sz, data_dim }/*, device*/),
								labels = torch::zeros({ sz, label_dim }/*, device*/);
	if (sz <= 0)
		return std::make_tuple(data, labels);

	//Monte Carlo method for the calculation of the acoustic pressure of the hadronic cascade
	TMCICAP acoustic_pressure(1e6);

	for (long pos{0}; pos < sz;)
	{
		//В источнике данных можно замешивать точки нескольких спектров для того чтобы подряд не отдавать все точки одного.
		//Еще как-то двигать бины частот например, двигая fmax и меняя число отчетов
		//Энергия распределена равномерно в логарифмическом масштабе
		TAcousticPressureParams params;
		//total energy of the hadronic cascade, eV
		params.E0 = static_cast<fp>(exp(log(Params.E0) + static_cast<float>(RandomDouble()) * (log(Params.E1) - log(Params.E0))));
		//speed of sound in a medium, g/cm2/s
		params.cs = static_cast<fp>(Params.cs0 + static_cast<float>(RandomDouble()) * (Params.cs1 - Params.cs0));
		//detector coordinates, g/cm2
		params.Rd = static_cast<fp>(Params.Rd0 + static_cast<float>(RandomDouble()) * (Params.Rd1 - Params.Rd0));
		//GetLMax(), detector coordinates, g/cm2
		params.Zd = static_cast<fp>(Params.Zd0 + static_cast<float>(RandomDouble()) * (Params.Zd1 - Params.Zd0));
		params.LonCutThreshold = Params.LonCutThreshold;		//longitudinal cutoff threshold, g/cm2
		params.TrCutThreshold = Params.TrCutThreshold;			//transverse cutoff threshold, g/cm2
		params.Parameterization = Params.Parameterization;
		fp fmax = static_cast<fp>(Params.fmax0 + static_cast<float>(RandomDouble()) * (Params.fmax1 - Params.fmax0));
		size_t n = static_cast<size_t>(Params.N0 + static_cast<float>(RandomDouble()) * (Params.N1 - Params.N0));

		acoustic_pressure.SetConditions(params);

		std::vector<fp> mpwx(n);
		std::vector<std::complex<fp>> pwx(n);
		//Acoustic pressure depending on sound frequency (frequency spectrum)
		std::cout<<"acoustic_pressure.GetComplexPwSeries"<<" "<<params.Rd<<"  "<< params.Zd<<" "<< fmax<<" "<< n<<" "<< pwx.size()<<" "<< mpwx.size() <<" "<< params.E0 <<std::endl;
		acoustic_pressure.GetComplexPwSeries(params.Rd, params.Zd, fmax, n, pwx, mpwx);

		//Однообразно случайно перемешаем массивы
		std::vector<int> indices(n);
		auto seed = RandomInt();
		for (long i = 0; i < static_cast<long>(n); i++)
			indices[i] = i;
		shuffle(pwx.begin(), pwx.end(),	std::default_random_engine(seed));
		shuffle(mpwx.begin(), mpwx.end(),	std::default_random_engine(seed));
		shuffle(indices.begin(), indices.end(),	std::default_random_engine(seed));

		n = (pwx.size() <= sz - pos) ? pwx.size() : sz - pos;
		std::cout<<pos<<" "<<sz<<" "<<n<<std::endl;
		for (long i = 0; i < static_cast<long>(n); i++)
		{
			data[pos + i][0] = log10(params.E0);
			data[pos + i][1] = params.cs;
			data[pos + i][2] = params.Rd;
			data[pos + i][3] = params.Zd;
			data[pos + i][4] = indices[static_cast<size_t>(i)] * fmax / n;
			labels[pos + i][0] = static_cast<float>(log10(static_cast<double>(pwx[static_cast<size_t>(i)].real()) + Params.E1));
			labels[pos + i][1] = static_cast<float>(log10(static_cast<double>(pwx[static_cast<size_t>(i)].imag()) + Params.E1));

			if (!finitef(labels[pos + i][0].item<float>()))
				labels[pos + i][0] = 0.f;
			if (!finitef(labels[pos + i][1].item<float>()))
				labels[pos + i][1] = 0.f;
		}
		pos += n;
	}
	return std::make_tuple(data, labels);
}

TNNAPDataset :: TNNAPDataset(const TNNAPDatasetParams &params) : Params(params)
{
	std::tie(Data, Labels) = GenerateData(params.DataSize);
}

TNNAPDataset::Example TNNAPDataset :: get(size_t index)
{
	auto d = Data[index];
	auto l = Labels[index];
	return {d, l};
}

torch::optional<size_t> TNNAPDataset :: size() const
{
	return Data.size(0);
}
