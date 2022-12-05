#pragma once

#include <vector>

#include "TorchHeader.h"
#include "Trainable.h"

#include "TAcousticPressure.h"

///E, cs, Rd, Zd, f/w -> RePw, ImPw
struct NNAPImpl : public Trainable
{
	torch::nn::Linear fc1,
		fc2,
		//fc3,
		fc4;
	torch::nn::LeakyReLU leaky_relu1,
		leaky_relu2;
		//leaky_relu3;
	torch::nn::BatchNorm1d batch_norm1,
		batch_norm2;
		//batch_norm3;
	//torch::nn::Sigmoid sigmoid;
	//torch::nn::Tanh tanh;

	NNAPImpl(const int input_dim = 5,const int output_dim = 2, const int hidden_dim = 64, const std::string module_name = "nnap")
		: Trainable(module_name),
		fc1(torch::nn::LinearOptions(input_dim, hidden_dim).bias(true)),
		leaky_relu1(torch::nn::LeakyReLUOptions().negative_slope(0.2)),
		batch_norm1(hidden_dim),
		fc2(torch::nn::LinearOptions(hidden_dim, hidden_dim).bias(true)),
		leaky_relu2(torch::nn::LeakyReLUOptions().negative_slope(0.2)),
		batch_norm2(hidden_dim),
		//fc3(torch::nn::LinearOptions(hidden_dim, hidden_dim).bias(true)),
		//leaky_relu3(torch::nn::LeakyReLUOptions().negative_slope(0.2)),
		//batch_norm3(hidden_dim),
		fc4(torch::nn::LinearOptions(hidden_dim, output_dim).bias(true))
		//leaky_relu4(torch::nn::LeakyReLUOptions().negative_slope(0.2))
		////batch_nor4(output_dim)
	{
		// register_module() is needed if we want to use the parameters() method later on
		//register_module(module_name + "_input", input);
		register_module(module_name + "_fc1", fc1);
		register_module(module_name + "_fc2", fc2);
		//register_module(module_name + "_fc3", fc3);
		register_module(module_name + "_fc4", fc4);
		register_module(module_name + "_leaky_relu1", leaky_relu1);
		register_module(module_name + "_leaky_relu2", leaky_relu2);
		//register_module(module_name + "_leaky_relu3", leaky_relu3);
		register_module(module_name + "_batch_norm1", batch_norm1);
		register_module(module_name + "_batch_norm2", batch_norm2);
		//register_module(module_name + "_batch_norm3", batch_norm3);
	}

	virtual ~NNAPImpl() override {}

	virtual torch::Tensor forward(torch::Tensor x) override
	{
		x = batch_norm1(leaky_relu1(fc1(x)));
		x = batch_norm2(leaky_relu2(fc2(x)) + x);
		//x = batch_norm3(leaky_relu3(fc3(x)) + x);
		x = fc4(x);
		return x;
	}
};

TORCH_MODULE(NNAP);



///Neural net approximation
class TNNAP : public TAcousticPressure {
protected:
	NNAP Net;
	torch::Device Device;
	float E1;
public:
	///NNAP is wrapper over std::shared_ptr
	TNNAP(NNAP net, const torch::Device &device, const float &e1) : Net(net), Device(device), E1(e1){}

	///Re или Im часть акустического эффекта (на данной частоте w) развития ЯЭК в точке распложения детектора (Rd,Zd)
	///Can generate whole spectrum at once. Use of GetComplexPwSeries is more preferably.
	virtual std::complex<fp> GetComplexPwValue (const fp &rd, const fp &zd, const fp &w) override;

	///Ряд значений из N точек акустического эффекта (на данной частоте w) развития ЯЭК в точке распложения детектора (Rd,Zd)
	///Спектр акустического импульса, ограниченный частотой fmax
	///Can generate whole spectrum at once
	virtual void GetComplexPwSeries(	const fp &rd,
																		const fp &zd,
																		fp &fmax,
																		size_t N,
																		std::vector<std::complex<fp>> &pwx,		//Комплексный спектр
																		std::vector<fp> &mpwx                 //Модуль спектра
																	) override;
};		//TNNAP



///******************************************************************************
//TNNAP
//******************************************************************************/

///Re или Im часть акустического эффекта (на данной частоте w) развития ЯЭК в точке распложения детектора (Rd,Zd)
///Can generate whole spectrum at once. Use of GetComplexPwSeries is more preferably.
std::complex<fp> TNNAP :: GetComplexPwValue (const fp &rd, const fp &zd, const fp &w)
{
	std::complex<fp> result;

	Net->eval();
	torch::NoGradGuard no_grad;

	auto params = FwWrapper->GetParams();
	params.Rd = rd;
	params.Zd = zd;
	params.w = w;

	const int sz = 1;					//spectrum size == 1 for single value
	const int data_dim = 5,		//E, cs, Rd, Zd, f/w
						label_dim = 2;	//RePw, ImPw
	torch::Tensor data = torch::zeros({ sz, data_dim }/*, device*/),
								labels = torch::zeros({ sz, label_dim }/*, device*/);
	for (int i = 0; i < sz; i++)			//sz == 1
	{
		data[i][0] = log10(params.E0);
		data[i][1] = params.cs;
		data[i][2] = params.Rd;
		data[i][3] = params.Zd;
		data[i][4] = static_cast<float>(w / 2. / TMath::Pi());
	}

	torch::Tensor output = Net->forward(data.to(Device));

	for (int i = 0; i < sz; i++)	//sz == 1
	{
		result.real(static_cast<fp>(pow(10., static_cast<double>(output[i][0].item<float>())) - static_cast<double>(E1)));
		result.imag(static_cast<fp>(pow(10., static_cast<double>(output[i][1].item<float>())) - static_cast<double>(E1)));
	}

	return result;
}

//Ряд значений из N точек акустического эффекта (на данной частоте w) развития ЯЭК в точке распложения детектора (Rd,Zd)
//Спектр акустического импульса, ограниченный частотой fmax
//Can generate whole spectrum at once
void TNNAP :: GetComplexPwSeries(	const fp &rd,
																	const fp &zd,
																	fp &fmax,
																	size_t N,
																	std::vector<std::complex<fp>> &pwx,		//Комплексный спектр
																	std::vector<fp> &mpwx                 //Модуль спектра
																)
{
	if (pwx.size() != N)
		pwx.resize(N);

	if (mpwx.size() != N)
		mpwx.resize(N);

	Net->eval();
	torch::NoGradGuard no_grad;

	auto params = FwWrapper->GetParams();
	params.Rd = rd;
	params.Zd = zd;

	const int data_dim = 5,		//E, cs, Rd, Zd, f/w
						label_dim = 2;	//RePw, ImPw
	torch::Tensor data = torch::zeros({ static_cast<long>(N), data_dim }/*, device*/),
								labels = torch::zeros({ static_cast<long>(N), label_dim }/*, device*/);

	for (int64_t i = 0; i < static_cast<int64_t>(N); i++)
	{
		float f = static_cast<float>(i) * static_cast<float>(fmax) / N;
		//pwx[i] = this->GetComplexPwValue(rd, zd, w);
		data[i][0] = log10(params.E0);
		data[i][1] = params.cs;
		data[i][2] = params.Rd;
		data[i][3] = params.Zd;
		data[i][4] = f;
	};

	torch::Tensor output = Net->forward(data.to(Device));

	for (int64_t i = 0; i < static_cast<int64_t>(N); i++)
	{
		pwx[static_cast<unsigned long>(i)].real(static_cast<fp>(pow(10., static_cast<double>(output[i][0].item<float>())) - static_cast<double>(E1)));
		pwx[static_cast<unsigned long>(i)].imag(static_cast<fp>(pow(10., static_cast<double>(output[i][1].item<float>())) - static_cast<double>(E1)));
		mpwx[static_cast<unsigned long>(i)] = static_cast<fp>(sqrt(pow(pwx[i].real(), 2) + pow(pwx[i].imag(), 2)));
	};
}
