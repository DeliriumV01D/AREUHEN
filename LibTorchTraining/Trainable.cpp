#include "Trainable.h"

void Trainable :: Initialize()
{
	for (auto &p : this->named_parameters())
	{
		if (p.key().find("norm") != p.key().npos && p.key().find(".weight") != p.key().npos)
		{
			this->named_parameters()[p.key()] = torch::nn::init::constant_(p.value(), 1.);
		} else if (p.key().find(".weight") != p.key().npos)
		{
			this->named_parameters()[p.key()] = torch::nn::init::xavier_normal_(p.value(), 0.1);
		}

		if (p.key().find(".bias") != p.key().npos)
		{
			this->named_parameters()[p.key()] = torch::nn::init::constant_(p.value(), 0.);
		}
	}
}
