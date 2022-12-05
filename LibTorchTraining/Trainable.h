#pragma once

#include "TorchHeader.h"

struct Trainable : public torch::nn::Module
{
	template <typename T>				//torch::nn::ModuleHolder<ImplType>
	static int ParamsCount(T module);

	Trainable (const std::string module_name) : torch::nn::Module(module_name){}
	virtual ~Trainable(){}
	virtual void Initialize();
	virtual torch::Tensor forward(torch::Tensor x) = 0;
};

template <typename T>				//torch::nn::ModuleHolder<ImplType>
int Trainable :: ParamsCount(T module)
{
	int result = 0;
	for (auto p : module->parameters())
	{
		int ss = 1;
		for (auto s : p.sizes())
			ss *= s;
		result += ss;
	}
	return result;
}
