#include <cmath>
#include <cstdio>
#include <iostream>
#include <filesystem>
#include <string>
#include <tuple>

#include <iostream>
#include <iomanip>

#include "TCanvas.h"
#include "TMath.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"

#include "TNNAPDataset.h"
#include "TNNAP.h"
#include "Trainer.h"
#include "TensorBoard.h"


int main(int argc, char ** argv)
{
	//setlocale(LC_ALL, "Russian");
	if (argc < 3)
	{
		std::cout << "You must choose action -train or -inference and specify the path to the checkpoint file"<<std::endl;
		return -1;
	}

	std::string CheckpointFile = std::string(argv[2]);
	TRandomInt::Instance().Initialize(static_cast<unsigned long>(std::time(nullptr)));
	TRandomDouble::Instance().Initialize(RandomInt());
	torch::manual_seed(42);
	TRootApp app("App", &argc, argv);


	// Create the device we pass around based on whether CUDA is available.
	torch::Device device(torch::kCPU);
	if (torch::cuda::is_available())
	{
		std::cout << "CUDA is available! Training on GPU." << std::endl;
		device = torch::Device(torch::kCUDA);
	} else	{
		std::cout << "CUDA is not available! Training on CPU." << std::endl;
	}

	if (std::string(argv[1]) == "-inference")
	{
		std::cout << "Inference" << std::endl;
		auto canvas = app.CreateCanvas("NNAP Inference", "Acoustic pressure", 200, 10, 1600, 800);
		canvas->Divide(2);
		canvas->Draw();

		const fp	rd = 100000,
							zd = 0;
		fp fmax = static_cast<fp>(1.e5);
		fp E1 = static_cast<fp>(1.e20-1);
		const size_t n = 128;
		std::vector<fp> mpwx(n),
										pt(n);
		std::vector<std::complex<fp>> pwx(n);

		TAcousticPressureParams params;
		params.E0 = static_cast<fp>(1.e19-1);				//total energy of the hadronic cascade, eV
		params.cs = static_cast<fp>(145000);			//speed of sound in a medium, g/cm2/s
		params.Rd = rd;						//detector coordinates, g/cm2
		params.Zd = zd;						//GetLMax(), detector coordinates, g/cm2
		params.LonCutThreshold = static_cast<fp>(1000);		//longitudinal cutoff threshold, g/cm2
		params.TrCutThreshold = static_cast<fp>(20);			//transverse cutoff threshold, g/cm2
		params.Parameterization = TEDDParameterization::CORSIM;

		//Monte Carlo method for the calculation of the acoustic pressure of the hadronic cascade
		//test of unuran passing as input a multi-dimension distribution object
		NNAP nnap(5, 2, 256, "nnap");
		nnap->to(device);

		std::cout << "nnap:\n" << std::flush;
		for (auto k : nnap->named_parameters())
			std::cout << k.key() << std::endl;

		std::cout << "nnap parameters count: " << Trainable::ParamsCount(nnap) << std::endl;
		torch::load(nnap, CheckpointFile);
		TNNAP acoustic_pressure(nnap, device, E1);

		acoustic_pressure.SetConditions(params);

		//Acoustic pressure depending on sound frequency (frequency spectrum)
		acoustic_pressure.GetComplexPwSeries(rd, zd, fmax, n, pwx, mpwx);

		app.GetCanvas(0)->cd(1);

		TH1D	Pw("Pw", "Pw-title", n, 0, Double_t(n-1) * Double_t(fmax) / n / 1000),
					Pt("Pt", "Pt-title", n, 0, Double_t(n-1) / Double_t(fmax));

		for (size_t i = 0; i < n; i++ )
			Pw.SetBinContent(static_cast<Int_t>(i), static_cast<Double_t>(mpwx[i]) * abs(/*(*(PRZt->Gamma))(4,0,isBaikal)*/.1/6.25/1e18));	// /6.25/1E18 ->Pa
		Pw.SetBins(n-1, 0, static_cast<Double_t>(n-1) * static_cast<Double_t>(fmax) / n / 1000);		//Hz->kHz

		//Acoustic pressure depending on time
		acoustic_pressure.GetPtSeries(pwx, pt, fmax, n);

		for (size_t i = 0; i < n; i++ )
			Pt.SetBinContent(static_cast<Int_t>(i), static_cast<Double_t>(pt[i]) * /*(*(PRZt->Gamma))(4,0,isBaikal)*/.1/6.25/1e18);// /6.25/1E18 -> Pa
		Pt.SetBins(n-1, 0, static_cast<Double_t>(n-1)/static_cast<Double_t>(fmax));
		app.GetCanvas(0)->cd(1);
		ApplyAxisStyle(&Pt,"P(t)", "t,sec", "P,Pa");
		Pt.Draw("L");
		app.GetCanvas(0)->cd(2);
		ApplyAxisStyle(&Pw, "|P(f)|","f,kHz","|P(f)|, Hz#upointkg#upointm^{-1}");
		Pw.Draw("C");

		app.GetCanvas(0)->Update();
		app.GetCanvas(0)->Draw();
		gSystem->ProcessEvents();

		try {		//scope
			app.Run(kFALSE);
		} catch (std::exception &e) {
			std::cerr<<e.what();
		};
	}

	if (std::string(argv[1]) == "-train")
	{
		std::cout << "Training" << std::endl;
		const size_t	BATCH_SIZE = 1024;
		bool RESTORE_FROM_CHECKPOINT = false;

		NNAP nnap(5, 2, 256, "nnap");
		nnap->to(device);

		std::cout << "nnap:\n" << std::flush;
		for (auto k : nnap->named_parameters())
			std::cout << k.key() << std::endl;

		std::cout << "nnap parameters count: " << Trainable::ParamsCount(nnap) << std::endl;

		TNNAPDatasetParams nnap_dataset_params;
		nnap_dataset_params.DataSize = 10000000;
		auto train_dataset = TNNAPDataset(nnap_dataset_params).map(torch::data::transforms::Stack<>());
		std::cout << "prepared dataset of size: " << train_dataset.size().value() << std::endl;

		nnap_dataset_params.DataSize = 100000;
		auto val_dataset = TNNAPDataset(nnap_dataset_params).map(torch::data::transforms::Stack<>());
		std::cout << "prepared dataset of size: " << val_dataset.size().value() << std::endl;

		auto train_data_loader = torch::data::make_data_loader(
			std::move(train_dataset),
			torch::data::DataLoaderOptions().batch_size(BATCH_SIZE).workers(1)
		);

		auto val_data_loader = torch::data::make_data_loader(
			std::move(val_dataset),
			torch::data::DataLoaderOptions().batch_size(BATCH_SIZE).workers(1)
		);

		torch::optim::Adam optimizer(nnap->parameters(), torch::optim::AdamOptions(1e-3).weight_decay(0.001).betas(std::make_tuple (0.9, 0.999)));
		torch::optim::StepLR sheduler(optimizer, 30, 0.1);		//Проверить как sheduler взаимодействует с загрузкой optimizer'а из чекпоинта
		auto loss_function = torch::nn::MSELoss();		//MSE
		auto acc_function = [](const torch::Tensor &t1, const torch::Tensor &t2){return t1.sub(t2).abs().sum();};

		if (RESTORE_FROM_CHECKPOINT)
		{
			std::cout << "restoring parameters from checkpoint..." << std::endl;
			torch::load(nnap, CheckpointFile);
		}	else {
			std::cout << "initializing parameters with xavier..." << std::endl;
			nnap->Initialize();
		}

		TensorBoardRoot tensor_board(&app);

		TrainerParams trainer_params;
		trainer_params.NumberOfEpochs = 100;
		trainer_params.LogInterval = 30;
		trainer_params.CheckpointEvery = 100;
		trainer_params.NetCheckpointPath = CheckpointFile;
		//trainer_params.OptimizerCheckpointPath  = "/DATA/Delirium/PROJECTS/AREUHEN/NNAPTraine/nnap-optimizer-checkpoint.pt";

		Trainer trainer(trainer_params);
		trainer.Train(nnap, device,	train_data_loader, val_data_loader, optimizer,	sheduler,	loss_function, acc_function, &tensor_board);

		std::cout << "Training complete!" << std::endl;

		try {		//scope
			app.Run(kFALSE);
		} catch (std::exception &e) {
			std::cerr<<e.what();
		};
	}


}
